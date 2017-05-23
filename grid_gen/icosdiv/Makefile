CC = g++

all: icosdiv

icosdiv: icosdiv.o Point.o Triangle.o
	$(CC) -o icosdiv icosdiv.o Point.o Triangle.o

icosdiv.o: Point.o Triangle.o icosdiv.cc
	$(CC) -c icosdiv.cc

Point.o: Point.cc Point.h
	$(CC) -c Point.cc

Triangle.o: Triangle.cc Triangle.h
	$(CC) -c Triangle.cc

clean:
	rm -f icosdiv.o Point.o Triangle.o icosdiv
