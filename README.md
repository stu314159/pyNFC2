# pyNFC2
second version of pyNFC - really just a do-over.
Starting over with pyNFC; smaller steps

Workflow: 
a) define geometry with script using FluidChannel.py as illustrated in pyNFC_preprocess.py
b) set problem-specific parameters (params.lbm) by invoking genInput.py either from command-line or script
c) run pyNFC_test.py from command line
d) once run is complete, run pyNFC_postprocess.py using mpi4py with your local MPI execution script
e) visualize with a tool that reads legacy VTK data files (like ParaView)
