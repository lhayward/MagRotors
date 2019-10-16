CC = c++ 
CFLAGS = -std=c++0x -Wall -funroll-loops -O3 #-pg #-xHOST

PROG = mag
HDRS = MersenneTwister.h FileReading.h SimParams.h
SRCS = Main.cpp FileReading.cpp SimParams.cpp
OBJS = Main.o FileReading.o SimParams.o

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(PROG)
	
FileReading.o: FileReading.cpp $(HDRS)
	$(CC) $(CFLAGS) -c FileReading.cpp -o FileReading.o

SimParams.o: SimParams.cpp $(HDRS)
	$(CC) $(CFLAGS) -c SimParams.cpp -o SimParams.o
	
Main.o: Main.cpp $(HDRS)
	$(CC) $(CFLAGS) -c Main.cpp -o Main.o 
	
clean:
	rm -f $(PROG) $(OBJS) 
