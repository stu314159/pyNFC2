#!/bin/bash

# input
python ./genInput.py sphere.mat D3Q15 25 1 0 1 10 1e-2 0 0

# run the code
python ./pyNFC2_test.py

