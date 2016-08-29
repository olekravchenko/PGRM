#!/bin/bash

gcc -c tasks.c -O3
gcc -c main.c -O3

gcc -o main main.o tasks.o -lgsl -lgslcblas -lm -g -pthread -O3


