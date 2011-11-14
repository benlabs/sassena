#!/bin/bash

# get rid of interference
export LD_LIBRARY_PATH=

PREFIX=$HOME/sw/moldyn

export FFTW3_DIR=$HOME/sw/moldyn/opt/fftw-3.3
export BOOST_ROOT=$HOME/sw/moldyn/opt/boost-1.47-shared-mt
export HDF5_HOME=$HOME/sw/moldyn/opt/hdf5-1.87

export PATH=/shared/apps/gnu4.4.4/openmpi/1.4.3/bin:$PATH

export LD_LIBRARY_PATH=${PREFIX}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib64:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=/opt/ofed/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/shared/apps/gnu4.4.4/openmpi/1.4.3/lib:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=${FFTW3_DIR}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${HDF5_HOME}/lib:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=${BOOST_ROOT}/lib:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH


cmake -D CMAKE_C_COMPILER=gcc44 -D CMAKE_CXX_COMPILER=g++44 -D CMAKE_EXE_LINKER_FLAGS="-Wl,-rpath,$LD_LIBRARY_PATH" ../..