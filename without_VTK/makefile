
CC = g++
CFLAGS = -g 


DEFAULT: FluidElement.o Grid.o Main.o
	$(CC) $(CFLAGS) -o Main.exe FluidElement.o Grid.o Main.o

FluidElement.o: ../sources/FluidElement.h ../sources/FluidElement.cpp 
	$(CC) $(CFLAGS) -o FluidElement.o -c ../sources/FluidElement.cpp

Grid.o: ../sources/Grid.h ../sources/Grid.cpp 
	$(CC) $(CFLAGS) -o Grid.o -c ../sources/Grid.cpp

Main.o: ../sources/Main.cpp 
	$(CC) $(CFLAGS) -o Main.o -c ../sources/Main.cpp




clean: 
	$(RM) *.o *~ *.exe

