#!/bin/bash

gcc -c tasks.c -O3 -Wall
gcc -c main.c -O3 -Wall

gcc -o main main.o tasks.o -lgsl -lgslcblas -lm -g -pthread -O3


