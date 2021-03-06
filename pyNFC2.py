#pyNFC2.py
"""
implementation file for pyNFC2 module

"""

import numpy as np
import pyLattice_nb as pl_nb

class pyNFC_LBM(object):
    def __init__(self,Nx,Ny,Nz,rho_lbm,u_bc,omega,Cs,lattice_type='D3Q15',):
        """
          Nx, Ny, Nz so the partition has info about the overall
                     lattice structure
          lattice_type - ['D3Q15' | 'D3Q19' | 'D3Q27']
          rho_lbm - scaled density for outlet boundary condition
          u_bc - scaled velocity for inlet boundary condition
          omega - relaxation constant for LBM collisions
          Cs - parameter for turbulence model
        """
        
        self.Nx = Nx; self.Ny = Ny; self.Nz = Nz; # global domain structure
        self.nnodes = self.Nx*self.Ny*self.Nz;

        # LBM simulation parameters
        self.rho_lbm = rho_lbm; self.u_bc = u_bc; self.omega = omega; self.Cs = Cs
        
        
        if lattice_type == 'D3Q15':
            self.lattice = pl_nb.D3Q15Lattice(self.Nx, self.Ny, self.Nz)
        elif lattice_type == 'D3Q19':
            self.lattice = pl_nb.D3Q19Lattice(self.Nx, self.Ny, self.Nz)
        else:
            self.lattice = pl_nb.D3Q27Lattice(self.Nx, self.Ny, self.Nz)

        self.numSpd = self.lattice.get_numSpd()
        self.ex = np.array(self.lattice.get_ex(),dtype=np.int32);
        self.ey = np.array(self.lattice.get_ey(),dtype=np.int32);
        self.ez = np.array(self.lattice.get_ez(),dtype=np.int32);
        self.bb_Spd = np.array(self.lattice.get_bbSpd(),dtype=np.int32);
        self.w = np.array(self.lattice.get_w(),dtype=np.float32);

        self.gen_adjacency() 
        self.allocate_data_arrays()
        self.initialize_node_lists()
        self.initialize_lattice_points()

        self.vtk_dump_num = 0;
        self.vtk_ux_stub = 'ux'; self.vtk_uy_stub = 'uy'; self.vtk_uz_stub = 'uz'
        self.vtk_rho_stub = 'density'
        self.vtk_suffix = '.b_dat'

    def initialize_lattice_points(self):
        """
        set density distribution values to zero-speed lattice everywhere
        """
        for nd in range(self.nnodes):
            for spd in range(self.numSpd):
                self.fEven[nd,spd] = self.rho_lbm*self.w[spd]
                self.fOdd[nd,spd]=self.rho_lbm*self.w[spd]
    
    def get_XYZ_index(self,g_nd): # this will depend upon a global geometry structure like a brick
        """
         get X,Y, and Z index for a given
        """
        z = g_nd/(self.Nx*self.Ny)
        y = (g_nd - z*self.Nx*self.Ny)/self.Nx
        x = (g_nd - z*self.Nx*self.Ny - y*self.Nx)
        return (x,y,z)

    def initialize_node_lists(self):
        """
         initialize inl, onl, snl
        """
        inl_filename = "inl.lbm"
        inl_f = open(inl_filename,'r')
        numINL = int(inl_f.readline());
        for i in range(numINL):
            gIN = int(inl_f.readline());
            self.inl[gIN] = 1
        inl_f.close()

        onl_filename = "onl.lbm"
        onl_f = open(onl_filename,'r')
        numONL = int(onl_f.readline())
        for k in range(numONL):
            gOUT = int(onl_f.readline());
            self.onl[gOUT] = 1
        onl_f.close()
        
        snl_filename = "snl.lbm"
        snl_f = open(snl_filename,'r');
        numSNL = int(snl_f.readline());
        for i in range(numSNL):
            gSNL=int(snl_f.readline());
            self.snl[gSNL] = 1
        snl_f.close()

        self.all_nodes = range(self.nnodes) 

    def allocate_data_arrays(self):
        """
         allocate arrays for LBM simulation
        """
        # some thought/testing should be done regarding the shape of this data array.
        self.fEven = np.empty([self.nnodes , self.numSpd],dtype=np.float32)
        self.fOdd = np.empty_like(self.fEven)
        self.snl = np.zeros([self.nnodes],dtype=np.int32);
        self.inl = np.zeros([self.nnodes],dtype=np.int32);
        self.onl = np.zeros([self.nnodes],dtype=np.int32);

    def get_gInd_XYZ(self,x,y,z): 
        """
         generate global node number based on input x y and z
        """
        return x+y*self.Nx + z*self.Nx*self.Ny

    def gen_adjacency(self):
        """
            generate the adjacency list indicating who is neighbor to whom
        """ 
        self.adjacency = np.empty((self.nnodes,self.numSpd),dtype=np.int32) # create the array
        
        # do this the bone-headed way:
        for k in range(self.nnodes):
            (x,y,z) = self.get_XYZ_index(k)
            for spd in range(self.numSpd):
                x_t = x + self.ex[spd]; x_t = x_t % self.Nx;
                y_t = y + self.ey[spd]; y_t = y_t % self.Ny;
                z_t = z + self.ez[spd]; z_t = z_t % self.Nz;
                self.adjacency[k,spd] = self.get_gInd_XYZ(x_t,y_t,z_t)

    def take_LBM_timestep(self,isEven):
        """
            carry out the LBM process for a time step

        """
        self.process_lattice_points(isEven,self.all_nodes)


    def process_lattice_points(self,isEven,lp_list):
        """
          carry out the LBM process for a list of lattice points
          isEven - boolean to indicate if this is an even time step or odd time step
          lp_list - list of lattice points to be processed 
        """
        if isEven:
            fIn = self.fEven; fOut = self.fOdd;
        else:
            fIn = self.fOdd; fOut = self.fEven;

        for lp in lp_list:
            f = fIn[lp,:]
            ndType = 0
            if self.inl[lp]==1:
                ndType = 2
            if self.onl[lp]==1:
                ndType = 3
            if self.snl[lp]==1:
                ndType = 1
           

            f_o = self.lattice.compute_fOut(f,ndType,self.omega,self.Cs,self.u_bc,self.rho_lbm)
            self.stream(fOut,f_o,lp)

    def stream(self,fOut,f,lp):
        """
            stream collided particle density distributions to neighbor lattice points
        """

        for spd in range(self.numSpd):
            tgtNd = self.adjacency[lp,spd]
            fOut[tgtNd,spd] = f[spd]

    def write_data(self,isEven):
        """
          generate binary data files: ux[dump #].b_dat, uy[dump #].b_dat, 
          uz[dump #].b_dat, and density[dump #].b_dat
        """
        ux, uy, uz, rho = self.compute_local_data(isEven);

        # create file names
        ux_fn = self.vtk_ux_stub + str(self.vtk_dump_num) + self.vtk_suffix
        uy_fn = self.vtk_uy_stub + str(self.vtk_dump_num) + self.vtk_suffix
        uz_fn = self.vtk_uz_stub + str(self.vtk_dump_num) + self.vtk_suffix
        rho_fn = self.vtk_rho_stub + str(self.vtk_dump_num) + self.vtk_suffix
        
        ux.astype('float32').tofile(ux_fn)
        uy.astype('float32').tofile(uy_fn)
        uz.astype('float32').tofile(uz_fn)
        rho.astype('float32').tofile(rho_fn)
        self.vtk_dump_num += 1

    def compute_local_data(self,isEven):
        """
         compute macroscopic data from density distribution data for all lattice points

        """

        if isEven:
            f = self.fEven;
        else:
            f = self.fOdd;

        ux = np.zeros([self.nnodes],dtype=np.float32)
        uy = np.zeros([self.nnodes],dtype=np.float32)
        uz = np.zeros([self.nnodes],dtype=np.float32)
        rho = np.zeros([self.nnodes],dtype=np.float32)

        rho = np.sum(f,axis=1)
        ux = np.dot(f,self.ex)/rho
        uy = np.dot(f,self.ey)/rho
        uz = np.dot(f,self.ez)/rho
#        for lp in self.all_nodes:
#            for spd in range(self.numSpd):
#                rho[lp]+=f[lp,spd]
#                ux[lp]+=self.ex[spd]*f[lp,spd]
#                uy[lp]+=self.ey[spd]*f[lp,spd]
#                uz[lp]+=self.ez[spd]*f[lp,spd]
#            ux[lp]/=rho[lp]; uy[lp]/=rho[lp]; uz[lp]/=rho[lp]
        ux[np.where(self.snl[:self.nnodes]==1)]=0.;
        uy[np.where(self.snl[:self.nnodes]==1)]=0.;
        uz[np.where(self.snl[:self.nnodes]==1)]=0.;
        uz[np.where(self.inl[:self.nnodes]==1)]=self.u_bc;
        return ux, uy, uz, rho
