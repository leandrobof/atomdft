/*
 * main.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: leandro
 */

#include "Orbital.h"
#include "gauss.h"
#include "SKtable.h"
#include "scf.h"

using namespace std;


//double Hubbard(Scf &a,vector<Orbital*> &b,const int &c );


int main(int argc,char *argv[]){



double tinf=50;
double t0=-8;
double h=0.008;




Scf scf(t0,tinf,h);



scf.initialize(argv[1]);
scf.run(0,1,1,0.2);

//scf.run(Atom,atof(argv[2]),atof(argv[3]),atof(argv[4]),0.4);

double e[3]={0.,0.,0.};
double ocupation[3]={0.,0.,0.};
double U[3]={0.,0.,0.};

scf.energy(e,ocupation);


//scf.run(Atom,atof(argv[2]),atof(argv[3]),atof(argv[4]),0.2);


//potencial de referencia V0



//*****

scf.run(atof(argv[2]),atof(argv[3]),atof(argv[4]),0.4);

//******
//Array  splines de los orbitales.
Orbital_spline *C[3]={NULL,NULL,NULL};
scf.orbital(C);


Potential_spline Veff;
Potential_spline vconf;

scf.Vconf(vconf);
scf.Veff_noconf(Veff);

string simbolo(argv[1]);

SKtable sk;

sk.create(C,C);

sk.run(C,C,Veff,vconf,e,U,ocupation,simbolo+"-"+simbolo+".skf");

for(int i=0;i<3;i++){
    delete C[i];
}



//*******************
return 0;
}


/*
double Hubbard(Scf &a,vector<Orbital*> &b,const int &c ){
	double h=0.001;
	double n=b[c]->ocup();
	double nf=n+0.5*h;
	double nb=n-0.5*h;
	double e=b[c]->energy();
	double ef;
	double eb;
	b[c]->set_ocup(nf);
	a.run(b,0,1,1,0.3);
	ef=b[c]->energy();
	b[c]->set_ocup(nb);
	a.run(b,0,1,1,0.3);
	eb=b[c]->energy();
	b[c]->set_ocup(n);
	a.run(b,0,1,1,0.3);
	return (ef-eb)/h;
}
*/

