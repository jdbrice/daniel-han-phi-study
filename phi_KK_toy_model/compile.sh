#!/bin/bash
g++ -c Random_routines.cpp $(root-config --cflags --libs)
g++ -c analysis.cpp $(root-config --cflags --libs)
g++ -o toy_model Random_routines.o analysis.o $(root-config --cflags --libs)
