all: main2.o Orbital.o Func.o xc_potential.o hartree.o Poisson.o Integrand.o gauss.o scf.o
	g++ main2.o scf.o Orbital.o Func.o xc_potential.o hartree.o Poisson.o Integrand.o gauss.o -pg  -o   atomdft -L/opt/etsf/lib -lgsl -lopenblas -lm -lxc -std=c++11

main2.o: main2.cpp Orbital.h
	g++  main2.cpp -c -I. -I/opt/etsf/include -std=c++11

gauss.o:gauss.cpp gauss.h
	g++ gauss.cpp -c -I. -std=c++11

Integrand.o:Integrand.cpp Integrand.h
	g++ Integrand.cpp -c -I. -std=c++11

scf.o:scf.cpp scf.h
	g++ scf.cpp -c -I/opt/etsf/include -std=c++11 

Orbital.o: Orbital.cpp Orbital.h 
	g++  Orbital.cpp -c -O3 -I.  -std=c++11

Func.o:Func.cpp Func.h 
	g++ Func.cpp -c -O3 -I. -std=c++11

xc_potential.o: xc_potential.cpp xc_potential.h
	g++ xc_potential.cpp   -c -O3 -I/opt/etsf/include 

hartree.o: hartree.cpp hartree.h Poisson.cpp Poisson.h
	g++ hartree.cpp -c -O3 -I/opt/boost_1_62_0/

Poisson.o: Poisson.cpp Poisson.h
	g++ Poisson.cpp -c -O3

clean:
	rm -rf *.o
	rm scf
