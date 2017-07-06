#!/bin/bash

# pre-process
python ./pyNFC_preprocess.py

# input
python ./genInput.py sphere.mat D3Q15 100 1 0 10 10 1e-2 0 0

# run the code
python ./pyNFC2_test.py


# post-process the output
mpirun -np 4 ./pyNFC_postprocess.py

