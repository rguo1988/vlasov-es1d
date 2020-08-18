# Makefile for g++

CC = g++
SRC = main.cpp diagnose.cpp poisson_solver.cpp cubic_spline_1d_tdma.cpp simulation.cpp
OBJ := $(SRC:.cpp=.o)
CFLAGS = -O3 -std=c++11 -march=native -fopenmp
TARGET = vlasov

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET)
$(OBJ): %.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	-rm $(OBJ) $(TARGET)
run: $(TARGET)
	-./$(TARGET)
