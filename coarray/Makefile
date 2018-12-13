CC=ifort

F90_FILES= jacobicoarray.f90 
F90_OBJ=jacobicoarray.o
APP=jacobicoarray
OPTIONS=-coarray=shared -warn all 
all:
	$(CC) $(OPTIONS)  -c $(F90_FILES)
	$(CC) $(OPTIONS)  $(F90_OBJ) -o $(APP)
clean:
	rm -f *.o
	rm -f $(APP)
	rm -f *.mod
	rm -f *.out
	rm -f *.error


