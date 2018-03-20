/*
 * SKtable.cpp
 *
 *  Created on: Feb 27, 2018
 *      Author: leandro
 */

#include "SKtable.h"


SKtable::SKtable() {
	// TODO Auto-generated constructor stub
	step=0.02;
	rmin=0.4;
	rmax=12;
	ngrid=50;
	S=vector<Integrand*>(10,NULL);
	V=vector<Integrand*>(10,NULL);
	homonuclear=true;
}

SKtable::~SKtable() {
	// TODO Auto-generated destructor stub

for(int i=0;i<S.size();i++){
			delete S[i];
			delete V[i];
		}
}

void SKtable::create( Orbital_spline **a,Orbital_spline **b,bool type){
    double e;
    homonuclear=type;
	if(a[0]!=NULL and b[0]!=NULL){
    	e=b[0]->energy();
		S[0]=new S_ss();
    	V[0]=new V_ss(e);
    }

    if(a[0]!=NULL and b[1]!=NULL){
    	e=b[1]->energy();
    	S[1]=new S_sp();
    	V[1]=new V_sp(e);
    }

    if(a[0]!=NULL and b[2]!=NULL){
    	e=b[2]->energy();
    	S[2]=new S_sd();
    	V[2]=new V_sd(e);
    }

    if(a[1]!=NULL and b[1]!=NULL){
    	e=b[1]->energy();
    	S[3]=new S_pp_pi();
    	V[3]=new V_pp_pi(e);
    	S[4]=new S_pp_sig();
    	V[4]=new V_pp_sig(e);
    }

    if(a[1]!=NULL and b[2]!=NULL){
    	e=b[2]->energy();
    	S[5]=new S_pd_pi();
    	V[5]=new V_pd_pi(e);
    	S[6]=new S_pd_sig();
    	V[6]=new V_pd_sig(e);
    }

    if(a[2]!=NULL and b[2]!=NULL){
    	e=b[2]->energy();
    	S[7]=new S_dd_del();
    	V[7]=new V_dd_del(e);
    	S[8]=new S_dd_pi();
    	V[8]=new V_dd_pi(e);
    	S[9]=new S_dd_sig();
    	V[9]=new V_dd_sig(e);
    }


};

void SKtable::run( Orbital_spline **A, Orbital_spline **B,Potential_spline &Veffa,Potential_spline &Veff0b,Potential_spline &Veffb,Potential_spline &Vconf,double *e,double *U,double *ocupation,string archivo){
	double d=rmin/2;

	gauss g(ngrid,d);

	ofstream salida2(archivo.c_str());

	salida2<<step<<" "<<int((rmax-step)/step)<<endl;

	if(homonuclear==true){
		salida2<<e[2]<<" "<<e[1]<<"  "<<e[0]<<" 0.0  "<<"Udd"<<"  "<<"Upp"<<"  "<<"Uss"<<"  "<<ocupation[2]<<"  "<<ocupation[1]<<"  "<<ocupation[0]<<endl;
	}
	salida2<<"12.01, 19*0.0"<<endl;

	int j=0;
	while(step*(1+j)<rmin){
		salida2<<"20*1.0,"<<endl;
	    j++;
	}

	double overlap[10];
	double H[10];
	double s[10];
	double v[10];
	double econf[10];

	for(int i=0;i<10;i++){
			H[i]=0.00;
			overlap[i]=0.00;
			s[i]=0.0;
			v[i]=0.0;
			if(V[i]!=NULL){
				econf[i]=V[i]->energy();
			}
			else {econf[i]=0.;}
	}

	while(2*d<rmax){
		g.integrate2d(Veffa,Veffb, Vconf,Veff0b,A,B,S,V,s,v);

	    for(int i=0;i<10;i++){
	    	if(V[i]!=NULL){
	    		overlap[i]=s[i];
	    		H[i]=econf[i]*overlap[i]+v[i];
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
	    d=d+step/2;
	    g.update_a(d);


	}



	salida2.close();

};

void SKtable::set_grid(int a){ngrid=a;}
void SKtable::set_rmax(double a){rmax=a;}
void SKtable::set_rmin(double c){rmin=c;}
void SKtable::set_step(double z){step=z;}
