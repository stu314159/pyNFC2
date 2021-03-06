#pyNFC_preprocess.py
"""
pre-processor file for pyNFC2.  Define geometry and discretization parameters

"""

import FluidChannel as fc


Lx_p = 1.0
Ly_p = 1.0
Lz_p = 5.0

s_x = 0.5
s_y = 0.5
s_z = 2.0

r_s = 0.25
N_divs = 7

geom_file_stub = 'sphere'
sphereObst = fc.SphereObstruction(r_s,s_x,s_y,s_z)
channel = fc.FluidChannel(Lx_p = Lx_p,Ly_p = Ly_p,Lz_p = Lz_p,
                          N_divs = N_divs,obst = sphereObst, fluid = 'water')
channel.write_bc_vtk()
channel.write_mat_file(geom_file_stub)
