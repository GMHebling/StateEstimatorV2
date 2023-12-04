CC=clang++
CFLAGS=-std=gnu++11 -Wall -g -I/Users/gustavohebling/Code/Projects/SolverEngine/
LDFLAGS= -lcholmod -lspqr -lsuitesparseconfig -lm -lstdc++ "-Wl,-rpath,/usr/local/lib/"

testEstimator: main.o data_structures.h electrical_utils.h estimator.h input_data.h mat_utils.h optimization_utils.h output.h topology_utils.h wls_utils.h sparseSystem.h engine.h
	$(CC) $(CFLAGS) -o testEstimator main.o $(LDFLAGS)

main.o: main.cc 
	$(CC) $(CFLAGS) -c main.cc $(LDFLAGS)