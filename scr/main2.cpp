/*
 * main.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: leandro
 */

#include "Orbital.h"
#include "gauss.h"

#include "scf.h"

using namespace std;

void create( Orbital_spline **a,Orbital_spline **b,vector<Integrand*> &c,vector<Integrand*> &d,Potential_spline &va,Potential_spline &vb,Potential_spline &vc);
void borrar(vector<Integrand*> &,vector<Integrand*> &);
double Hubbard(Scf &a,vector<Orbital*> &b,const int &c );
void SK_table(Orbital_spline **A,Orbital_spline **B,Potential_spline &Veff,Potential_spline &Veff0,Potential_spline &Vconf,double *e,double *U,double *ocupation,string archivo);

int main(int argc,char *argv[]){


vector <Orbital*> Atom;



double tinf=50;
double t0=-8;
double h=0.001;




Scf scf(t0,tinf,h);



scf.initialize(Atom,argv[1]);

scf.run(Atom,0,1,1,0.3);

//scf.run(Atom,atof(argv[2]),atof(argv[3]),atof(argv[4]),0.4);

double e[3]={0.,0.,0.};
double ocupation[3]={0.,0.,0.};
double U[3]={0.,0.,0.};

scf.energy(Atom,e,ocupation);


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

scf.run(Atom,atof(argv[2]),atof(argv[3]),atof(argv[4]),0.2);

//*****
//Array  splines de los orbitales.
Orbital_spline *C[3]={NULL,NULL,NULL};
scf.orbital(Atom,C);

for(int i=0;i<N;i++){
	veff2[i]=veff[i]-vconf[i];
}

Potential_spline Veff(veff2,r,N);
Potential_spline Vconf(vconf,r,N);

string simbolo(argv[1]);

SK_table(C,C,Veff,Veff0,Vconf,e,U,ocupation,simbolo+"-"+simbolo+".skf");


delete C[0];
delete C[1];
delete C[2];

for (int i=0;i<Atom.size();i++){
	delete Atom[i];
}

//*******************
return 0;
}


void create( Orbital_spline **a,Orbital_spline **b,vector<Integrand*> &c,vector<Integrand*> &d,Potential_spline &va,Potential_spline &vb,Potential_spline &vc){
    if(a[0]!=NULL and b[0]!=NULL){
    	c[0]=new S_ss(a[0],b[0]);
    	d[0]=new V_ss(a[0],b[0],va,vb,vc);
    }

    if(a[0]!=NULL and b[1]!=NULL){
    	c[1]=new S_sp(a[0],b[1]);
    	d[1]=new V_sp(a[0],b[1],va,vb,vc);
    }

    if(a[0]!=NULL and b[2]!=NULL){
    	c[2]=new S_sd(a[0],b[2]);
    	d[2]=new V_sd(a[0],b[2],va,vb,vc);
    }

    if(a[1]!=NULL and b[1]!=NULL){
    	c[3]=new S_pp_pi(a[1],b[1]);
    	d[3]=new V_pp_pi(a[1],b[1],va,vb,vc);
    	c[4]=new S_pp_sig(a[1],b[1]);
    	d[4]=new V_pp_sig(a[1],b[1],va,vb,vc);
    }

    if(a[1]!=NULL and b[2]!=NULL){
    	c[5]=new S_pd_pi(a[1],b[2]);
    	d[5]=new V_pd_pi(a[1],b[2],va,vb,vc);
    	c[6]=new S_pd_sig(a[1],b[2]);
    	d[6]=new V_pd_sig(a[1],b[2],va,vb,vc);
    }

    if(a[2]!=NULL and b[2]!=NULL){
    	c[7]=new S_dd_del(a[2],b[2]);
    	d[7]=new V_dd_del(a[2],b[2],va,vb,vc);
    	c[8]=new S_dd_pi(a[2],b[2]);
    	d[8]=new V_dd_pi(a[2],b[2],va,vb,vc);
    	c[9]=new S_dd_sig(a[2],b[2]);
    	d[9]=new V_dd_sig(a[2],b[2],va,vb,vc);
    }


};



void borrar(vector<Integrand*> &c,vector<Integrand*> &d){
	for(int i=0;i<c.size();i++){
		delete c[i];
		delete d[i];
	}
}
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

void SK_table(Orbital_spline **A,Orbital_spline **B,Potential_spline &Veff,Potential_spline &Veff0,Potential_spline &Vconf,double *e,double *U,double *ocupation,string archivo){
vector<Integrand*> S(10,NULL);
vector<Integrand*> V(10,NULL);

double d=0.2;
create(A,B,S,V,Veff,Veff0,Vconf);


gauss g(50,d);

ofstream salida2(archivo);

salida2<<"0.02"<<" "<<"599"<<endl;
salida2<<e[2]<<" "<<e[1]<<"  "<<e[0]<<" 0.0  "<<U[2]<<"  "<<U[1]<<"  "<<U[0]<<"  "<<ocupation[2]<<"  "<<ocupation[1]<<"  "<<ocupation[0]<<endl;
salida2<<"12.01, 19*0.0"<<endl;

for(int i=0;i<19;i++){
	salida2<<"20*1.0,"<<endl;
}

double overlap[10];
double H[10];
double s[10];
double v[10];

for(int i=0;i<10;i++){
		H[i]=0.00;
		overlap[i]=0.00;
		s[i]=0.0;
		v[i]=0.0;
}

while(2*d<12){
	g.integrate2d(0,1,0,pi,S,V,s,v);

    for(int i=0;i<10;i++){
    	if(V[i]!=NULL){
    		overlap[i]=s[i];
    		H[i]=V[i]->energy()*overlap[i]+v[i];
    	}
	}


	for(int i=9;i>=0;i--){
		salida2<<H[i]<<"  ";
	}

	for(int i=9;i>=0;i--){
			salida2<<overlap[i]<<"  ";
	}
    for(int i=0;i<10;i++){
    	s[i]=0.;
    	v[i]=0.;

    }


	salida2<<endl;
    d=d+0.01;
    g.update_a(d);


}


borrar(S,V);
salida2.close();
}

