#!/bin/bash
g++ -c Random_routines.cpp $(root-config --cflags --libs)
g++ -c rho_test.cpp $(root-config --cflags --libs)
g++ -o rho_test Random_routines.o rho_test.o $(root-config --cflags --libs)
