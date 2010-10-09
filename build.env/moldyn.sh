#!/bin/bash

PREFIX=$HOME/sw/moldyn

# requirements:
# boost
# fftw
# hdf5
# libxml2
# 

export BOOST_ROOT=$PREFIX

cmake ../..
