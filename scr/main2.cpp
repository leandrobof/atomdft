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

int N=scf.Nt();
double *r=scf.rgrid();
double *veff=scf.Veff_noconf();
double *vconf=scf.Vconf();
double veff0[N];
double veff2[N];

//potencial de referencia V0
for(int i=0;i<N;i++){
	veff0[i]=veff[i]-vconf[i];
}

Potential_spline Veff0(veff0,r,N);
//*****

scf.run(atof(argv[2]),atof(argv[3]),atof(argv[4]),0.4);

//******
//Array  splines de los orbitales.
Orbital_spline *C[3]={NULL,NULL,NULL};
scf.orbital(C);

for(int i=0;i<N;i++){
	veff2[i]=veff[i]-vconf[i];
}

Potential_spline Veff(veff2,r,N);
Potential_spline Vconf(vconf,r,N);

string simbolo(argv[1]);

SKtable sk;

sk.create(C,C);

sk.run(C,C,Veff,Veff0,Vconf,e,U,ocupation,simbolo+"-"+simbolo+".skf");

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

