#!/bin/bash

gcc -c tasks.c
gcc -c main.c

gcc -o main main.o tasks.o -lgsl -lgslcblas -lm -g -pthread -O3


