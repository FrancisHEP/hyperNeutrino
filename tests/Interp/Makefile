FLAGS=-std=c++11 -Wall -O2 $(shell root-config --cflags)
LIBS=$(shell root-config --libs) -lMinuit -lMathCore -lMathMore

all: main

main: main.o
	g++ $(FLAGS) main.o $(LIBS) -o main

main.o: main.cpp
	g++ $(FLAGS) -c main.cpp -o main.o


