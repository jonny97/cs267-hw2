#!/bin/csh -f
#
# this is a comment
#
make
mpirun -n 1 -c 1 ./serial -n 500 -s serial.txt