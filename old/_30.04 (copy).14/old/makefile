
CC = g++
CFLAGS = -g -Wall 


VTKLIBS = /usr/local/lib/
VTKINCLUDE = /usr/local/include/vtk-6.1/
INCLUDES = -I. -I /usr/X11R6/include/ -I $(VTKINCLUDE)
LDFLAGS = -L $(VTKLIBS) -L /usr/X11R6/lib64/ -L /usr/X11R6/lib/


DEFAULT: FluidElement.o Grid.o Main.o
	$(CC) $(CFLAGS) -o Main.exe FluidElement.o Grid.o Main.o

FluidElement.o: ../sources/FluidElement.h ../sources/FluidElement.cpp 
	$(CC) $(CFLAGS) -o FluidElement.o -c ../sources/FluidElement.cpp

Grid.o: ../sources/Grid.h ../sources/Grid.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFlAGS) -o Grid.o -c ../sources/Grid.cpp

Main.o: ../sources/Main.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFlAGS) -o Main.o -c ../sources/Main.cpp




clean: 
	$(RM) *.o *~ *.exe

