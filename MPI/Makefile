CC=mpiifort
F90_FILES=jacobimpi.f90
F90_OBJ= jacobimpi.o
APP=jacobimpi
OPTIONS=-heap-arrays 100 -check bounds -warn all -traceback
all:
	$(CC) $(OPTIONS) -c $(F90_FILES)
	$(CC) $(OPTIONS) $(F90_OBJ) -o $(APP)
clean:
	rm -f *.o
	rm -f $(APP)
	rm -f *.mod