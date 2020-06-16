#!/bin/bash
#

./mkmf -t Makefile.template -p ../lib/libflowsampler.a ../src/*/*.F90

make all

make examples

#make install

