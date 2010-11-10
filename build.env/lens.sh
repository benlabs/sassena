#!/bin/bash

module () 
{ 
    eval `/usr/bin/modulecmd bash $*`
}
    
module unload PE-pgi
module load PE-gnu/4.3.2
module unload ompi
module load ompi/1.4.2-gnu4.3.2  
module load git
module load cmake
module list

    
PREFIX=$HOME/sw/lens

# requirements:
# boost
# fftw
# hdf5
# libxml2
# 

cmake -D CMAKE_EXE_LINKER_FLAGS="-Wl,-rpath,$LD_LIBRARY_PATH" -D CMAKE_PREFIX_PATH=$PREFIX ../..


