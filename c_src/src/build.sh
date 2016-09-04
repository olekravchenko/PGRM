#!/bin/bash

gcc -c tasks.c -O3 -Wall -fopenmp
gcc -c main.c -O3 -Wall -fopenmp

gcc -o main main.o tasks.o -lgsl -lgslcblas -lm -g -pthread -O3 -fopenmp


