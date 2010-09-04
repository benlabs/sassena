#!/bin/bash

module () 
{ 
 eval `/opt/modules/3.1.6/bin/modulecmd bash $*`
}
    
PREFIX=$HOME/sw/kraken-xt5

# requirements:
# boost
# fftw
# hdf5
# libxml2
# 

module swap PE-pgi PE-gnu
module load hdf5
module load fftw

export BOOST_ROOT=$PREFIX

cmake -D STATIC=true -D CRAY=true ../..



