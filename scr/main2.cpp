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

double e[3]={0.,0.,0.};
double ocupation[3]={0.,0.,0.};
double U[3]={0.,0.,0.};

double tinf=50;
double t0=-8;
double h=0.008;




Scf atomo1(t0,tinf,h);
Scf atomo2(t0,tinf,h);


atomo1.initialize(argv[1]);
atomo1.run(0,1,1,0.2);

atomo2.initialize(argv[2]);
atomo2.run(0,1,1,0.2);

//scf.run(Atom,atof(argv[2]),atof(argv[3]),atof(argv[4]),0.4);



atomo1.energy(e,ocupation);
atomo2.energy(e,ocupation);

//scf.run(Atom,atof(argv[2]),atof(argv[3]),atof(argv[4]),0.2);


//potencial de referencia V0



//*****

atomo1.run(atof(argv[3]),atof(argv[4]),atof(argv[5]),0.4);
atomo2.run(atof(argv[6]),atof(argv[7]),atof(argv[8]),0.4);

//******
//Array  splines de los orbitales.
Orbital_spline *A[3]={NULL,NULL,NULL};
Orbital_spline *B[3]={NULL,NULL,NULL};

atomo1.orbital(A);
atomo2.orbital(B);

Potential_spline Veff_A;
Potential_spline vconf_A;

Potential_spline Veff_B;
Potential_spline vconf_B;

atomo1.Vconf(vconf_A);
atomo1.Veff_noconf(Veff_A);

atomo2.Vconf(vconf_B);
atomo2.Veff_noconf(Veff_B);


string simboloA(argv[1]);
string simboloB(argv[2]);

SKtable skaa;
SKtable skab;
SKtable skba;
SKtable skbb;

skab.create(A,A);
skab.create(A,B);
skba.create(B,A);
skab.create(B,B);

skab.run(A,A,Veff_A,vconf_A,e,U,ocupation,simboloA+"-"+simboloA+"2.skf");
skab.run(A,B,Veff_A,vconf_B,e,U,ocupation,simboloA+"-"+simboloB+"2.skf");
skba.run(B,A,Veff_B,vconf_A,e,U,ocupation,simboloB+"-"+simboloA+"2.skf");
skab.run(B,B,Veff_B,vconf_B,e,U,ocupation,simboloB+"-"+simboloB+"2.skf");

for(int i=0;i<3;i++){
    delete A[i];
    delete B[i];
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

