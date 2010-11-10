#!/bin/bash

PREFIX=$HOME/sw/viper

cmake -D CMAKE_EXE_LINKER_FLAGS="-Wl,-rpath,$LD_LIBRARY_PATH" -D CMAKE_PREFIX_PATH=$PREFIX ../..
