CC=g++
CFLAGS=-std=gnu++11 -Wall -g

testEstimator: main.o data_structures.h electrical_utils.h estimator.h input_data.h mat_utils.h optimization_utils.h output.h topology_utils.h wls_utils.h
	$(CC) $(CFLAGS) -o testEstimator main.o $(LDFLAGS)

main.o: main.cc 
	$(CC) $(CFLAGS) -c main.cc 