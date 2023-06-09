main: main.f90 const.o io.o simulation.o
	gfortran -o main main.f90 const.o io.o simulation.o -O5
	
const.o: const.f90
	gfortran -c const.f90
	
io.o: io.f90
	gfortran -c io.f90
	
simulation.o: simulation.f90
	gfortran -c simulation.f90 -O5
	
clean:
	rm -f *.mod *.o
