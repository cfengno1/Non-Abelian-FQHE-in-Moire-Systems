#!bin sh
ifort -qopenmp constant.f90 basiscon.f90 hamilton.f90 zndrv1.f90 test.f90 -mcmodel=large  -g -CB -traceback  -mkl -L ~/ -larpack
#nohup ./a.out&
