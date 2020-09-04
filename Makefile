# Makefile for g++

CC = g++
SRC = main.cpp diagnose.cpp poisson_solver.cpp cubic_spline_1d_tdma.cpp simulation.cpp
OBJ := $(SRC:.cpp=.o)
CFLAGS = -O3 -fopenmp -march=native -std=c++11
TARGET = vlasov

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET)
$(OBJ): %.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	-rm $(OBJ) $(TARGET)
<<<<<<< HEAD

run:$(TARGET)
	./$(TARGET)
=======
run: $(TARGET)
	-./$(TARGET)
>>>>>>> b723b0f9293571755d5b4acaa7e2aa76cb9b6edd
