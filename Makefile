# Makefile for g++

CC = g++
SRC = main.cpp diagnose.cpp poisson_solver.cpp cubic_spline_1d_tdma.cpp simulation.cpp
OBJ := diagnose.o poisson_solver.o cubic_spline_1d_tdma.o 
OBJ_D := main.o simulation.o
CFLAGS = -O3 -fopenmp -march=native -std=c++11
TARGET = vlasov
NTHREADS = 64

$(TARGET): $(OBJ) $(OBJ_D)
	$(CC) $(CFLAGS) $(OBJ) $(OBJ_D) -o $(TARGET)
$(OBJ): %.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJ_D): %.o: %.cpp input.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	-rm $(OBJ) $(OBJ_D) $(TARGET)
run: $(TARGET)
	-OMP_NUM_THREADS=$(NTHREADS) ./$(TARGET)
