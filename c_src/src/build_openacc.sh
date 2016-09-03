#!/bin/bash

gcc-6 -c tasks.c -O3 -Wall -fopenacc -fopenacc-dim=4:4
gcc-6 -c main.c -O3 -Wall -fopenacc -fopenacc-dim=4:4

gcc-6 -o main main.o tasks.o -lgsl -lgslcblas -lm -g -pthread -O3 -fopenacc -fopenacc-dim=4:4


