#FluidChannel.py
"""
Class implementation file for the Python class FluidChannel
Depends on vtkHelper module for geometry visualization functionality
"""
import math
import argparse
import numpy as np
from vtkHelper import saveStructuredPointsVTK_ascii as writeVTK
import scipy.io

class EmptyChannel:  
    """
     a channel with nothing in it
    """
    def __init__(self,Lo):
        """
         constructor
        """
        self.Lo = Lo


    def get_Lo(self):
        """
         set Lo if need be ?
        """
        return self.Lo

    def get_obstList(self,X,Y,Z):
        """
          for an empty channel - no obstacles 
        """
        return []


class SphereObstruction(EmptyChannel):
    """
     a channel with a sphere obstruction
    """

    def __init__(self,r,x_c,y_c,z_c):
        """
          just need to define the radius and position of the center of the obstacle.
          it is up to the caller to verify that the object will fit within the intended
          channel.  If it does not fit, the obstacle will effectively be
          truncated at the channel boundaries
         
        """
        self.r = r
        self.x_c = x_c
        self.y_c = y_c
        self.z_c = z_c
 
    def get_Lo(self):
        return self.r*2.

    def get_obstList(self,X,Y,Z):
        """
            return a list of all indices all indices within boundary of sphere 
        """

        x = np.array(X); y = np.array(Y); z = np.array(Z);
        dist = (x - self.x_c)**2 + (y - self.y_c)**2 + (z - self.z_c)**2
       
        return list(np.where(dist < self.r**2))
       

class GolfBall(EmptyChannel):
    """
     a channel with a golf ball obstacle
    """

    def __init__(self,SO,d_dimp,rd_dimp,N_e,N_a):
        """
           SO - pass in a sphericle obstacle as one of the arguments
           d_dimp = diameter of the dimples on the golf ball
           rd_dimp = radial distance of the center of the dimple from the center
                     of the golf ball
           N_e = number of dimples along all [0,pi] elevation angles 
           N_e = number of dimples along all [0,2pi] azimuthal angles
        """
        self.sphere = SO;
        self.d_dimp = d_dimp;
        self.rd_dimp = rd_dimp;
        self.N_e = N_e;
        self.N_a = N_a;

    def get_Lo(self):
        return self.sphere.get_Lo()


    def get_obstList(self,X,Y,Z):
        """
           return the obst list for the golf ball
        """
        obst_list1 = self.sphere.get_obstList(X,Y,Z)
        el_angles = np.linspace(0.,np.pi,self.N_e)
        
        x = np.array(X); y = np.array(Y); z = np.array(Z);
        print "removing the dimples"
        # start removing dimples
        iel = 0;
        for el in el_angles:
            iel+=1
        # for each elevation, we will get a different number of dimples
            N_az_el = np.floor(self.N_a*np.sin(el))+1;
            if N_az_el == 1:
                N_az_el+=1
            
            az_angles = np.linspace(0.,2.*np.pi, N_az_el, endpoint = False)
            print "removing dimples in elevation %g of %g" % (iel, len(el_angles))
            iaz = 0;
            for az in az_angles:
              iaz+=1
              print "removing dimple %g of %g on this elevation" % (iaz,len(az_angles))
              # get coordinates of the center of the spherical dimple
              y_c_d = self.sphere.y_c + self.rd_dimp*np.cos(el);
              z_c_d = self.sphere.z_c + self.rd_dimp*np.sin(az)*np.sin(el);
              x_c_d = self.sphere.x_c + self.rd_dimp*np.cos(az)*np.sin(el);
 
              dist = (x - x_c_d)**2 + (y - y_c_d)**2 + (z - z_c_d)**2
              dimples = np.where(dist <= ((self.d_dimp/2.))**2)
              obst_list1 = np.setxor1d(obst_list1[:],
                  np.intersect1d(obst_list1[:],dimples[:]))
             

        return obst_list1[:] 
        


# class EllipticalScourPit(EmptyChannel):
#     """
#      a channel with an elliptical scour pit with prescribed properties
#      corresponds to case 3 of Bryan's geometry_desc.m
#     """
# 
#     def __init__(self,x_c,z_c,cyl_rad):
#         """
#           constructor giving the x and z coordinates of the scour pit along with
#           the radius of the cylindrical piling
#         """
#         self.x_c = x_c
#         self.z_c = z_c
#         self.cyl_rad = cyl_rad
# 
#     def get_Lo(self):
#         return self.cyl_rad*2.
# 
#     def get_obstList(self,X,Y,Z):
#         """
#          return a list of all indices of lattice points within the boundaries of the
#          scour pit obstacle
#         """
#        
#         ellip_a = 2.*2.*self.cyl_rad
#         ellip_b = 2.*self.cyl_rad
#         ellip_c = 8.*self.cyl_rad
#         ellip_x = self.x_c
#         ellip_z = self.z_c + self.cyl_rad
#         ellip_y = ellip_b 
# 
#         floor_part = np.array(np.where(Y < ellip_b)).flatten()
# 
#         dist = (X - self.x_c)**2 + (Z - self.z_c)**2;
#         cyl_part = list(np.array(np.where( dist < self.cyl_rad**2)).flatten())
# 
#         scour_pit = np.array(np.where( (X - ellip_x)**2/(ellip_a**2) + 
#                         (Y - ellip_y)**2/(ellip_b**2) +
#                         (Z - ellip_z)**2/(ellip_c**2) <= 1.)).flatten()
# 
#         # remove the scour pit from the floor
#         obst_list = np.setxor1d(floor_part[:], 
#                         np.intersect1d(floor_part[:],scour_pit[:]))
# 
# 
#         # then add the cylinder
#         obst_list = np.union1d(obst_list[:],cyl_part[:])
#         
#         return list(obst_list[:])

# class ConeScourPit(EmptyChannel):
#     """
#     a channel with a conical scour pit determined by the angle of repose of the soil particles (assumed to be river sand, phi=30 deg).  
#     """
# 
#     def __init__(self,x_c,z_c,cyl_rad):
#         """
#           constructor giving the x and z coordinates of the scour pit along with the radius of the cylindrical piling
#         """
#         self.x_c = x_c
#         self.z_c = z_c
#         self.cyl_rad = cyl_rad
# 
#     def get_Lo(self):
#         return self.cyl_rad*2.
# 
#     def get_obstList(self,X,Y,Z):
#         """
#          return a list of all indices of lattice points within the boundaries of the conical scour pit obstacle.  x_s is defined in 'Scour at marine structures' by Richard Whitehouse, 1998.  Assumes river sand with phi (angle of repose) equal to 30 degrees.  h_cone is equal to rad_cone*tan(30) = rad_cone*0.57735
#         """
#         x_c_cone = self.x_c
# 	    z_c_cone = self.z_c
#         y_c_cone = 0.
#         x_s = 2.25*2*self.cyl_rad
#         rad_cone = x_s + self.cyl_rad
# 	    h_cone = rad_cone*0.57735
# 
#         floor_part = np.array(np.where(Y < h_cone)).flatten()
# 
#         dist = (X - self.x_c)**2 + (Z - self.z_c)**2;
#         cyl_part = list(np.array(np.where( dist < self.cyl_rad**2)).flatten())
# 
#         scour_pit = np.array(np.where( (X - x_c_cone)**2 + (Z - z_c_cone)**2 <= ((self.cyl_rad/cone)/(h_cone))**2*(Y - y_c_cone)**2))
# 
#         # remove the scour pit from the floor
#         obst_list = np.setxor1d(floor_part[:], 
#                         np.intersect1d(floor_part[:],scour_pit[:]))
# 
# 
#         # then add the cylinder
#         obst_list = np.union1d(obst_list[:],cyl_part[:])
#         
#         return list(obst_list[:])
# 
# class SinglePile(EmptyChannel):
#     """
#     a channel with a single pile, no scour.  Used for comparison to both elliptical and conical scour pits.
#     """
# 
#     def __init__(self,x_c,z_c,cyl_rad):
#         """
#           constructor giving the x and z coordinates of the piling center along with the radius of the cylindrical piling
#         """
#         self.x_c = x_c
#         self.z_c = z_c
#         self.cyl_rad = cyl_rad
# 
#     def get_Lo(self):
#         return self.cyl_rad*2.
# 
#     def get_obstList(self,X,Y,Z):
#         """
#          return a list of all indices of lattice points within the boundaries of the bed Bed thickness is equal to the diameter of the piling (2x radius)
#         """
#        
#     	#Bed
#         floor_part = np.array(np.where(Y < 2*self.cyl_rad)).flatten()
# 	
# 	#Piling
#         dist = (X - self.x_c)**2 + (Z - self.z_c)**2;
#         cyl_part = list(np.array(np.where( dist < self.cyl_rad**2)).flatten())
# 
# 
#         # then add the cylinder
#         obst_list = np.union1d(floor_part[:],cyl_part[:])
#         
#         return list(obst_list[:])
# 
# class WavyBed(EmptyChannel):
#     """
#     a channel with a single pile, Sin-wave bottom.  
#     """
# 
#     def __init__(self,x_c,z_c,cyl_rad):
#         """
#           constructor giving the x and z coordinates of the piling center along with the radius of the cylindrical piling
#         """
#         self.x_c = x_c
#         self.z_c = z_c
#         self.cyl_rad = cyl_rad
# 
#     def get_Lo(self):
#         return self.cyl_rad*2.
# 
#     def get_obstList(self,X,Y,Z):
#         """
# waveh and wavel are used to characterize the sine wave for the bed.  shallower sin waves do better in remaining stable throughout the simulation at low Reynolds numbers.
# 
#         """
#        
#     	#Bed
#         waveh = 0.125
#         wavel = 5        
#         floor_part = np.array(np.where(Y < (waveh*np.sin(wavel*Z) + 2*self.cyl_rad))).flatten()
#         
# 	#Piling
#         dist = (X - self.x_c)**2 + (Z - self.z_c)**2;
#         cyl_part = list(np.array(np.where( dist < self.cyl_rad**2)).flatten())
# 
# 
#         # then add the cylinder
#         obst_list = np.union1d(floor_part[:],cyl_part[:])
#         
#         return list(obst_list[:])
# 
# 
# class PipeContract(EmptyChannel):
#     """
#     a single smooth pipe with diameter in, diam_in, through a contraction and leaving at diameter out, diam_out.  Contraction assumed to be 45 degrees.  Channel assumed to be 2 x 2 x 8.  Lo = diam_out (smaller diameter).  Contraction begins at z = 4.  For a clean pipe, diam_in = 1.8 and diam_out = 0.8.
#     """
# 
#     def __init__(self,diam_in,diam_out):
#         """
#  constructor identifying diameters into and out of contraction.  Recommend diam_in = 1.8 and diam_out = 0.8
#         """
#         self.diam_in = diam_in
#         self.diam_out = diam_out
# 
#     def get_Lo(self):
#         return self.diam_out
# 
#     def get_obstList(self,X,Y,Z):
#         """
#    Define areas external to pipe.
#         """
#        #Pipe in - find all points exterior of large pipe
#         pipe_in = np.array(np.where((X - 1)**2 + (Y - 1)**2 > (self.diam_in/2)**2)).flatten()
#         pipe_in_stop = np.array(np.where(Z <= 4)).flatten()
#         pipe_in = np.intersect1d(pipe_in[:],pipe_in_stop[:])
#     
#         #Contraction - find all points exterior of contraction
#         r_cone = self.diam_out
#         h_cone = self.diam_out	
#         contraction = np.array(np.where((X - 1)**2 + (Y - 1)**2 > (r_cone/h_cone)**2*(Z - (4 + h_cone))**2)).flatten()
#         contraction_start = np.array(np.where(Z >= 4)).flatten()
#         contraction_stop = np.array(np.where(Z <= 4 + .5*self.diam_out)).flatten()
#         contraction = np.intersect1d(contraction[:],contraction_start[:])
#         contraction = np.intersect1d(contraction[:],contraction_stop[:])
#     
#         #Pipe out - final all points exterior of smaller pipe
#         pipe_out = np.array(np.where((X - 1)**2 + (Y - 1)**2 > (self.diam_out/2)**2)).flatten()
#         pipe_out_start = np.array(np.where(Z >= 4 + .5*self.diam_out)).flatten()
#         pipe_out = np.intersect1d(pipe_out[:],pipe_out_start[:])
#     
#     
#         #Put the pieces together
#     
#         #pipe = pipe_in[:]
#         pipe = np.union1d(contraction[:],pipe_in[:])
#         pipe = np.union1d(pipe[:],pipe_out[:])
#     
#         obst_list = pipe[:]
# 
#        
#         return list(obst_list[:])
# 
# class PipeExpand(EmptyChannel):
#     """
#     opposite of pipe contraction.  a single smooth pipe with diameter in, diam_in, through an expansion and leaving at diameter out, diam_out.  Expansion assumed to be 45 degrees.  Channel assumed to be 2 x 2 x 8.  Lo = diam_in (smaller diameter).  Expansion begins at z = 4.  Best works when diam_in = 0.8 and diam_out = 1.8
#     """
# 
#     def __init__(self,diam_in,diam_out):
#         """
#           constructor identifying pipe diameters into and out of expansion.  Recommend diam_in = 0.8 and diam_out = 1.8
#         """
#         self.diam_in = diam_in
# 	    self.diam_out = diam_out
# 
#     def get_Lo(self):
#         return self.diam_in
# 
#     def get_obstList(self,X,Y,Z):
#         """
#    Define areas external to pipe.
#         """
#        #Pipe in - find all points exterior of small
#         pipe_in = np.array(np.where((X - 1)**2 + (Y - 1)**2 > (self.diam_in/2)**2)).flatten()
#         pipe_in_stop = np.array(np.where(Z <= 1.5 + 0.5*(self.diam_out - self.diam_in))).flatten()
#         pipe_in = np.intersect1d(pipe_in[:],pipe_in_stop[:])
#     
#         #Expansion - find all points exterior of expansion
#         r_cone = self.diam_in
#         h_cone = self.diam_in	
#         expansion = np.array(np.where((X - 1)**2 + (Y - 1)**2 > (r_cone/h_cone)**2*(Z - 1.5)**2)).flatten()
#         expansion_start = np.array(np.where(Z >= 1.5 + 0.5*(self.diam_out - self.diam_in)))
#         #expansion_stop = np.array(np.where(Z <= 4)).flatten()
#         expansion = np.intersect1d(expansion[:],expansion_start[:])
#         #expansion = np.intersect1d(expansion[:],expansion_stop[:])
#     
#         #Pipe out - final all points exterior of smaller pipe
#         pipe_out = np.array(np.where((X - 1)**2 + (Y - 1)**2 > (self.diam_out/2)**2)).flatten()
#         pipe_out_start = np.array(np.where(Z >= 1.5 + 0.5*(self.diam_in - self.diam_out))).flatten()
#         pipe_out = np.intersect1d(pipe_out[:],pipe_out_start[:])
#     
#     
#         #Put the pieces together
#     
#         pipe = expansion[:]
#         pipe = np.union1d(expansion[:],pipe_in[:])
#         pipe = np.union1d(pipe[:],pipe_out[:])
#     
#         obst_list = pipe[:]
# 
#        
#         return list(obst_list[:])
# 
# class PipeTurn(EmptyChannel):
#     """
#   Provides an s-shaped pipe of constant radius with two 180-degree turns constructed out of constant-radius tori.  Diameter needs to be 0.5 for 
#     """
# 
#     def __init__(self,diam_in):
#         """
#           constructor providing pipe diameter for use in Lo.  Use 0.5.
#         """
#         self.diam_in = diam_in
# 	
#     def get_Lo(self):
#         return self.diam_in
# 
#     def get_obstList(self,X,Y,Z):
#         """
#    Define areas external to pipe.
#         """
#        #Pipe_1
#         pipe_1 = np.array(np.where((X - 1)**2 + (Y - 4)**2 >= 0.5**2)).flatten()
#         pipe_1_stop_z = np.array(np.where(Z <= 3.0)).flatten()
#         pipe_1_stop_y = np.array(np.where(Y >= 3.25)).flatten()
#         pipe_1_stop = np.intersect1d(pipe_1_stop_z[:],pipe_1_stop_y[:])
#         pipe_1 = np.intersect1d(pipe_1[:],pipe_1_stop[:])
#     
#         #Turn_1
#         turn_1 = np.array(np.where((0.75 - np.sqrt((Y - 3.25)**2 + (Z -3)**2))**2 + (X - 1)**2 >= 0.5**2)).flatten()
#         turn_1_stop_z = np.array(np.where(Z >= 3.0)).flatten()
#         turn_1_stop_y = np.array(np.where(Y>= 1.75)).flatten()
#         turn_1_stop = np.intersect1d(turn_1_stop_z[:],turn_1_stop_y[:])
#         turn_1 = np.intersect1d(turn_1[:],turn_1_stop[:])
#     
#         #Pipe_2
#         pipe_2 = np.array(np.where((X - 1)**2 + (Y - 2.5)**2 >= 0.5**2)).flatten()
#         pipe_2_start_z = np.array(np.where(Z >= 1.5)).flatten()
#         pipe_2_start_y_up = np.array(np.where(Y <= 3.25)).flatten()
#         pipe_2_start_y_down = np.array(np.where(Y >= 1.75)).flatten()
#         pipe_2_start_y = np.intersect1d(pipe_2_start_y_up[:],pipe_2_start_y_down[:])	
#         pipe_2_start = np.intersect1d(pipe_2_start_z[:],pipe_2_start_y[:])
#         pipe_2 = np.intersect1d(pipe_2[:],pipe_2_start[:])
#         pipe_2_stop_z = np.array(np.where(Z <= 3.0)).flatten()
#         pipe_2_stop_y = np.array(np.where(Y <= 3.25)).flatten()
#         pipe_2_stop = np.intersect1d(pipe_2_stop_z[:],pipe_2_stop_y[:])
#         pipe_2 = np.intersect1d(pipe_2[:],pipe_2_stop[:])
#     
#         #Turn_2
#         turn_2 = np.array(np.where((0.75 - np.sqrt((Y - 1.75)**2 + (Z -1.5)**2))**2 + (X - 1)**2 >= 0.5**2)).flatten()
#         turn_2_stop_z = np.array(np.where(Z <= 1.5)).flatten()
#         turn_2_stop_y = np.array(np.where(Y <= 3.25)).flatten()
#         turn_2_stop = np.intersect1d(turn_2_stop_z[:],turn_2_stop_y[:])
#         turn_2 = np.intersect1d(turn_2[:],turn_2_stop[:])
#         
#         #Pipe_3
#         pipe_3 = np.array(np.where((X - 1)**2 + (Y - 1.0)**2 >= 0.5**2)).flatten()
#         pipe_3_start_z = np.array(np.where(Z >= 1.5)).flatten()
#         pipe_3_start_y = np.array(np.where(Y <= 1.75)).flatten()
#         pipe_3_start = np.intersect1d(pipe_3_start_z[:],pipe_3_start_y[:])
#         pipe_3 = np.intersect1d(pipe_3[:],pipe_3_start[:])	
#     
#         #Put the pieces together
#     
#         pipe = np.union1d(pipe_1[:],turn_1[:])
#         pipe = np.union1d(pipe[:],pipe_2[:])
#         pipe = np.union1d(pipe[:],turn_2[:])	
#         pipe = np.union1d(pipe[:],pipe_3[:])
#     
#         obst_list = pipe[:]
# 
#         return list(obst_list[:])
# 
# class PipeOut(EmptyChannel):
#     """
#   Class consisting of a single pipe of diam_in and length length_in exiting a wall into an open space.  
#     """
# 
#     def __init__(self,diam_in,length_in):
#         """
#         defines the diameter and length (z axis) of pipe leading to open area
#         """
#         self.diam_in = diam_in
#         self.length_in = length_in
# 
#     def get_Lo(self):
#         return self.diam_in
# 
#     def get_obstList(self,X,Y,Z):
#         """
#    Define solid areas around pipe.  Everything else will be open.  Ensure coordinates for center of circle match center of Lx-Ly.
#         """
#        #Pipe In
#         pipe_in = np.array(np.where((X - 0.5*(4))**2 + (Y - 0.5*(4))**2 >= (0.5*self.diam_in)**2)).flatten()
#         pipe_in_stop = np.array(np.where(Z <= self.length_in)).flatten()
#         pipe_in = np.intersect1d(pipe_in[:],pipe_in_stop[:])
#     
#     
#         obst_list = pipe_in[:]
# 
#         return list(obst_list[:])
# 
# class Butterfly(EmptyChannel):
#     """
#   A geometry class that defines a fully open butterfly valve within a pipe of diam=1.0.  
#     """
# 
#     def __init__(self,diam):
#         """
#           constructor identifying pipe diameter.  Must be a 1 diam pipe inside a 1.2 x 1.2 x 8 channel.  Valve center at z = 3.
#         """
#         self.diam = diam
# 	
# 
#     def get_Lo(self):
#         return self.diam
# 
#     def get_obstList(self,X,Y,Z):
#         """
#    Define solid areas
# 
#         """
# 
#         #Pipe
#         pipe = np.array(np.where((X - 0.6)**2 + (Y - 0.6)**2 >= 0.5**2)).flatten()
#     
#         #Seat
#         seat = np.array(np.where((X - 0.6)**2 + (Y - 0.6)**2 >= 0.42**2)).flatten()
#         seat_start = np.array(np.where(Z >= 2.975)).flatten()
#         seat_stop = np.array(np.where(Z <= 3.025)).flatten()
#         seat = np.intersect1d(seat[:],seat_start[:])
#         seat = np.intersect1d(seat[:],seat_stop[:])
#     
#         #Pivot
#         pivot = np.array(np.where((X - 0.6)**2 + (Z - 3)**2 <= 0.075**2)).flatten()
#     
#         #Front Disc
#         front_disc = np.array(np.where((Y - 0.6)**2 + (Z - 3)**2 <= 0.5**2)).flatten()
#         front_disc_stop = np.array(np.where(Z <= 3.0)).flatten()
#         front_disc_x_min = np.array(np.where(X >= 0.525)).flatten()
#         front_disc_x_max = np.array(np.where(X <= 0.575)).flatten()
#     
#         front_disc = np.intersect1d(front_disc[:],front_disc_stop[:])
#         front_disc = np.intersect1d(front_disc[:],front_disc_x_min[:])
#         front_disc = np.intersect1d(front_disc[:],front_disc_x_max[:])
#     
#         #Back Disc
#         back_disc = np.array(np.where((Y - 0.6)**2 + (Z - 3)**2 <= 0.5**2)).flatten()
#         back_disc_start = np.array(np.where(Z >= 3.0)).flatten()
#         back_disc_x_min = np.array(np.where(X >= 0.625)).flatten()
#         back_disc_x_max = np.array(np.where(X <= 0.675)).flatten()
#     
#         back_disc = np.intersect1d(back_disc[:],back_disc_start[:])
#         back_disc = np.intersect1d(back_disc[:],back_disc_x_min[:])
#         back_disc = np.intersect1d(back_disc[:],back_disc_x_max[:])
#     
#         #Put the pieces together
#     
#         valve = np.union1d(pipe[:],seat[:])
#         valve = np.union1d(valve[:],pivot[:])
#         valve = np.union1d(valve[:],front_disc[:])
#         valve = np.union1d(valve[:],back_disc[:])
#         
#         obst_list = valve[:]
# 
#         return list(obst_list[:])


class Tee(EmptyChannel):
    """
  establishes a single large pipe with a "tee" into a smaller pipe that loops up and around before rejoining the main line.  The Main line undergoes a contraction after the tee but before the rejoining secondary line.  diam_2 should be smaller than diam_1.
    """

    def __init__(self,diam_1,diam_2):
        """
        Constructor identifying the diameters of the two pipes.  Pipe 1 runs straight through from Z_min to Z_max.  Pipe 2 tees off and runs parallel to Pipe 1.  Pipe 1 enters/exits z planes at y = 1.  Pipe 2 runs at y = 3.  Assumes dimensions of space (X,Y,Z) is (2,4,8).
        """
        self.diam_1 = diam_1
	self.diam_2 = diam_2
	

    def get_Lo(self):
        return self.diam_1

    def get_obstList(self,X,Y,Z):
        """
   Define solid areas

        """

        #Pipe 1
        pipe_1a = np.array(np.where((X - 1)**2 + (Y - 1)**2 <= (self.diam_1/2)**2)).flatten()
        pipe_1a_stop = np.array(np.where(Z<=4.)).flatten()
        pipe_1a = np.intersect1d(pipe_1a[:],pipe_1a_stop[:])	
        pipe_1b = np.array(np.where((X - 1)**2 + (Y - 1)**2 <= (self.diam_1/4)**2)).flatten()
        pipe_1 = np.union1d(pipe_1a[:],pipe_1b[:])
    
        #Pipe 2 Tee Off
        tee_1 = np.array(np.where((X - 1)**2 + (Z - 1.5)**2 <= (self.diam_2/2)**2)).flatten()
        tee_1_start = np.array(np.where(Y >= 1)).flatten()
        tee_1_end = np.array(np.where(Y <= 3 - 0.5*self.diam_2)).flatten()
        tee_1 = np.intersect1d(tee_1[:],tee_1_start[:])
        tee_1 = np.intersect1d(tee_1[:],tee_1_end[:])
    
        #Pipe 2 Elbow 1
        elbow_1 = np.array(np.where((self.diam_2/2 - np.sqrt((Y - (3 - self.diam_2/2))**2 + (Z -(1.5 + self.diam_2/2))**2))**2 + (X - 1)**2 <= (self.diam_2/2)**2)).flatten()
        elbow_1_start = np.array(np.where(Y >= 3- 0.5*self.diam_2)).flatten()
        elbow_1_stop = np.array(np.where(Z <= 1.5 + self.diam_2/2)).flatten()
        elbow_1 = np.intersect1d(elbow_1[:],elbow_1_start[:])
        elbow_1 = np.intersect1d(elbow_1[:],elbow_1_stop[:])
    
    
        #Pipe 2
        pipe_2 = np.array(np.where((X - 1)**2 + (Y - 3)**2 <= (self.diam_2/2)**2)).flatten()
        pipe_2_start = np.array(np.where(Z >= 1.5 + self.diam_2/2)).flatten()
        pipe_2_stop = np.array(np.where(Z <= 5 - self.diam_2/2)).flatten()
        pipe_2 = np.intersect1d(pipe_2[:],pipe_2_start[:])
        pipe_2 = np.intersect1d(pipe_2[:],pipe_2_stop[:])
    
        #Pipe 2 Elbow 2
        elbow_2 = np.array(np.where((self.diam_2/2 - np.sqrt((Y - (3 - self.diam_2/2))**2 + (Z -(5- self.diam_2/2))**2))**2 + (X - 1)**2 <= (self.diam_2/2)**2)).flatten()
        elbow_2_start = np.array(np.where(Y >= 3- 0.5*self.diam_2)).flatten()
        elbow_2_stop = np.array(np.where(Z >= 5- self.diam_2/2)).flatten()
        elbow_2 = np.intersect1d(elbow_2[:],elbow_2_start[:])
        elbow_2 = np.intersect1d(elbow_2[:],elbow_2_stop[:])
    
        #Pipe 2 Tee In
        tee_2 = np.array(np.where((X - 1)**2 + (Z - 5)**2 <= (self.diam_2/2)**2)).flatten()
        tee_2_start = np.array(np.where(Y >= 1)).flatten()
        tee_2_end = np.array(np.where(Y <= 3 - 0.5*self.diam_2)).flatten()
        tee_2 = np.intersect1d(tee_2[:],tee_2_start[:])
        tee_2 = np.intersect1d(tee_2[:],tee_2_end[:])
    
        empty = np.array(np.where(Y>=0.)).flatten() 
    
        #Put the pieces together
        pipe = np.union1d(pipe_1[:],tee_1[:])
        pipe = np.union1d(pipe[:],elbow_1[:])
        pipe = np.union1d(pipe[:],pipe_2[:])
        pipe = np.union1d(pipe[:],elbow_2[:])
        pipe = np.union1d(pipe[:],tee_2[:])
        pipe = np.setxor1d(pipe[:], empty[:])	
     
        obst_list = pipe[:]

        return list(obst_list[:])

def fluid_properties(fluid_str):  
   """
   Return the physical density and kinematic viscosity for the prescribed
   fluid.
   
   """
   fluid_lib = {'water':(1000., 1.0e-6), 
                'glycol':(965.3,6.216e-4),
                'glycerin':(1260,1.18e-3)}
   if fluid_str in fluid_lib.keys():
     return fluid_lib[fluid_str]
   else:
     print 'valid fluids are:'
     for keys in fluid_lib:
       print " '%s' " % keys
     raise KeyError('invalid fluid specified')

class FluidChannel:
    def __init__(self,Lx_p=1.,
        Ly_p=1.,
        Lz_p=6.,
        fluid='water', 
        obst=EmptyChannel(1.0),
        N_divs = 5,
        wallList=['left','right','top','bottom']):
        """
         class constructor
        """
        self.Lx_p = Lx_p
        self.Ly_p = Ly_p
        self.Lz_p = Lz_p
        self.N_divs = N_divs
        self.fluid = fluid
        self.obst = obst
        

        # generate the geometry

        Lo = obst.get_Lo()

        self.Ny = math.ceil((N_divs-1)*(Ly_p/Lo))+1
        self.Nx = math.ceil((N_divs-1)*(Lx_p/Lo))+1
        self.Nz = math.ceil((N_divs-1)*(Lz_p/Lo))+1
        self.nnodes = self.Nx*self.Ny*self.Nz
        print "Creating channel with %g lattice points." % self.nnodes
        x = np.linspace(0.,Lx_p,self.Nx).astype(np.float32);
        y = np.linspace(0.,Ly_p,self.Ny).astype(np.float32);
        z = np.linspace(0.,Lz_p,self.Nz).astype(np.float32);
        self.l_tol = (x[1]-x[0])/2.;
   
        Y,Z,X = np.meshgrid(y,z,x);
    
        self.x = np.reshape(X,self.nnodes)
        self.y = np.reshape(Y,self.nnodes)
        self.z = np.reshape(Z,self.nnodes)

        # get fluid properties from the included fluid library
        self.rho_p, self.nu_p = fluid_properties(fluid)

        # identify inlet and outlet nodes - 
        # require the user to set solid boundaries separately
        #self.inlet_list = np.where(self.z==0)
        self.inlet_list = np.where(np.isclose(self.z,0.,atol=self.l_tol))
        print "maximum value of self.z = %g, Lz_p = %g"%(np.amax(self.z),self.Lz_p)
        #self.outlet_list = np.where(self.z>=self.Lz_p)
        self.outlet_list = np.where(np.isclose(self.z,self.Lz_p,atol=self.l_tol))
        
        print "Getting obstacle list"
        # get obstacle list
        self.obst_list = self.obst.get_obstList(self.x[:],self.y[:],self.z[:])
        

        print "Generating channel solid boundaries"
        # set channel walls
        self.set_channel_walls(wallList)

        # now eliminate overlap between node lists

        self.inlet_list = np.setxor1d(self.inlet_list[:],
            np.intersect1d(self.inlet_list[:],self.solid_list[:]))
        self.inlet_list = np.setxor1d(self.inlet_list[:],
            np.intersect1d(self.inlet_list[:],self.obst_list[:]))
        
        self.outlet_list = np.setxor1d(self.outlet_list[:],
            np.intersect1d(self.outlet_list[:],self.solid_list[:]))
        self.outlet_list = np.setxor1d(self.outlet_list[:],
            np.intersect1d(self.outlet_list[:],self.obst_list[:]))

        self.obst_list = np.setxor1d(self.obst_list[:],
            np.intersect1d(self.obst_list[:],self.solid_list[:]))
       
    def write_mat_file(self, geom_filename):
        """
          generate the mat file to interface with genInput.py.  Needs to save
          Lx_p, Ly_p, Lz_p, Lo, Ny_divs, rho_p, nu_p, snl, inl and onl.
          note that the snl and obst_list need to be combined into one list 
        """
        mat_dict = {}
        mat_dict['Lx_p'] = self.Lx_p
        mat_dict['Ly_p'] = self.Ly_p
        mat_dict['Lz_p'] = self.Lz_p
        mat_dict['Lo'] = self.obst.get_Lo()
        mat_dict['Ny_divs'] = self.N_divs
        mat_dict['rho_p'] = self.rho_p
        mat_dict['nu_p'] = self.nu_p
        mat_dict['snl'] = list(np.union1d(self.obst_list[:],self.solid_list[:]))
        mat_dict['inl'] = list(self.inlet_list[:])
        mat_dict['onl'] = list(self.outlet_list[:])

        scipy.io.savemat(geom_filename,mat_dict)


    
    def write_bc_vtk(self):
        """
         write node lists to properly formatted VTK files
        """
        print "Creating boundary condition arrays"
        obst_array = np.zeros(self.nnodes)
        obst_array[list(self.obst_list)] = 100.

        #print type(self.inlet_list)
        inlet_array = np.zeros(self.nnodes)
        inlet_array[list(self.inlet_list)] = 200.

        outlet_array = np.zeros(self.nnodes)
        outlet_array[list(self.outlet_list)] = 300.

        solid_array = np.zeros(self.nnodes)
        solid_array[list(self.solid_list)] = 500.
        
        dims = [int(self.Nx), int(self.Ny), int(self.Nz)]
        origin = [0., 0., 0.]
        dx = self.x[1] - self.x[0]
        spacing = [dx, dx, dx] #uniform lattice
        
        print "Writing boundary conditions to VTK files"
        writeVTK(inlet_array,'inlet','inlet.vtk',dims,origin,spacing)
        writeVTK(outlet_array,'outlet','outlet.vtk',dims,origin,spacing)
        writeVTK(obst_array,'obst','obst.vtk',dims,origin,spacing)
        writeVTK(solid_array,'solid','solid.vtk',dims,origin,spacing)


     # must have geometry set first
    def set_channel_walls(self,walls=['left','right','top','bottom']): 
        """
         set up to 4 walls as solid walls for the simulation
        """
        solid_list_a = np.empty(0).flatten()
        solid_list_b = np.empty(0).flatten()
        solid_list_c = np.empty(0).flatten()
        solid_list_d = np.empty(0).flatten()

        for w in walls:
            if w=='right':
                #solid_list_a = np.array(np.where((self.x==0.))).flatten()
		solid_list_a = np.array(np.where(np.isclose(self.x,0.,atol=self.l_tol))).flatten()
            elif w=='left':
		
                #solid_list_b = np.array(np.where((self.x >= self.Lx_p))).flatten()
		solid_list_b = np.array(np.where(np.isclose(self.x,self.Lx_p,atol=self.l_tol))).flatten()
		#print "num nodes on left wall = %d" % len(solid_list_b)
            elif w=='top':
                #solid_list_d = np.array(np.where((self.y >= self.Ly_p))).flatten()
		solid_list_d = np.array(np.where(np.isclose(self.y,self.Ly_p,atol=self.l_tol))).flatten()
            elif w=='bottom':
                #solid_list_c = np.array(np.where((self.y == 0.))).flatten()
		solid_list_c = np.array(np.where(np.isclose(self.y,0.,atol=self.l_tol))).flatten()

        solid_list = np.array(np.union1d(solid_list_a,solid_list_b)); 
        solid_list = np.array(np.union1d(solid_list,solid_list_c))
        self.solid_list = np.array(np.union1d(solid_list,solid_list_d))

