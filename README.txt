This repository holds minimac.

To build optimized (objects in obj/ and bin in bin/):
  make
or
  make opt

To build with openmp (objects in obj/omp/ and bin in bin/omp/):
  make openmp

To build debug (objects in obj/debug/ and bin in bin/debug/):
  make debug

To build for profile (objects in obj/profile/ and bin in bin/profile/):
  make profile

To build everything (optimized, openmp, debug, and profile):
  make all

To clean (removes objects and binaries for all types):
  make clean

