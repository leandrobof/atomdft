/*
 * scf.cpp
 *
 *  Created on: 24/05/2017
 *      Author: leandro
 */
#include "scf.h"

void dens_inicial(double *r,double *rho,int N,int Z){
	double x;
	double a=0.7280642371;
	double b=-0.5430794693;
	double gamma=0.3612163121;
	double Zeff;
	double v;
	for(int i=0;i<N;i++){
		x=r[i]*pow((128*Z/9*pi*pi),1/3.);
		Zeff=Z*(1+a*sqrt(x)+b*x*exp(-gamma*sqrt(x)))*(1+a*sqrt(x)+b*x*exp(-gamma*sqrt(x)))*exp(-2*a*sqrt(x));
        v=Zeff/r[i];
        rho[i]=-1/(3*pi*pi)*pow(-2*v,3/2.);
	}
}

void new_rho(double *rho, double *grad,vector<Orbital*> &Atom,double alfa,int N){
	double nrho;
	double ngrad;
	for(int i=0;i<N;i++){
		nrho=0;
		ngrad=0;
		for (int j=0;j<Atom.size();j++){
	         nrho=nrho+(Atom[j])->operator[](i);
             ngrad=ngrad+(Atom[j]->grd(i));
	    }
		rho[i]=(1-alfa)*rho[i]+alfa*nrho;
		grad[i]=(1-alfa)*grad[i]+alfa*ngrad;

	}
};
Scf::Scf(int z,double tmin,double tmax,double step,bool rel,bool gga){
Z=z;
h=step;
N=(tmax-tmin)/h;

t=new double [N];
r=new double [N];
veff=new double [N];
vn=new double [N];
vconf=new double [N];
rho=new double [N];

grad=new double [N];
alfa=0.3;
W=0;
a=1;
r0=1;
Relativistic=rel;

t[0]=tmin;
if(gga==false){
		vxc=new Xc_lda(N,Relativistic);
	}
	else{
		vxc=new Xc_gga(N);
	}
};
void Scf::initialize(){


/***Genera malla uniforme t, y potencial*****************/
	for (int i=0;i<N;i++){
		t[i]=t[0]+i*h;
		r[i]=exp(t[i])/Z;

		vn[i]=-Z/r[i];
	    rho[i]=0;
	    vconf[i]=W/(1+exp(-a*(r[i]-r0)));
	    veff[i]=vn[i]+vconf[i];
	}

}
Scf::~Scf(){
	delete [] t;
	delete [] r;
	delete [] veff;
	delete [] vn;
	delete [] vconf;
	delete [] rho;

    delete [] grad;
    delete vxc;
}

void Scf::run(vector<Orbital*> &Atom,double w,double al,double ro,double alf){
	/**************************************************************/
	alfa=alf;
	W=w;
	a=al;
	r0=ro;
	for(int i=0;i<N;i++){
		vconf[i]=W/(1+exp(-a*(r[i]-r0)));


	}

	/*Se crea spline y Orbital.Se resuelve el atomo hidrogenoide veff=-Z/r */



	for(int i=0;i<Atom.size();i++){
		Atom[i]->resolver(veff,r,W);

	}


	/*************************************************************/
	//Se obtiene rho inicial a partir del atomo hidrogenoide

	new_rho(rho,grad,Atom, alfa, N);
	//Se crea e inicializa los potenciales de hartree e intercambio y correlacion.y se crea veff
	Hartree vh(r,rho,t,h,Z,N);

    vxc->update(r,rho,grad,N);
	for(int i=0;i<N;i++){
	        veff[i]=vn[i]+vh[i]+vxc->operator[](i)+vconf[i];
	     }



	for(int i=0;i<Atom.size();i++){
		Atom[i]->resolver(veff,r,W);

	}

	double *e1;
	double *e2;
	e1=new double [Atom.size()];
	e2=new double [Atom.size()];



	for(int i=0;i<Atom.size();i++){
		e1[i]=0.;
		e2[i]=Atom[i]->energy();

	}

	double error=0.00000001;
	int iteraciones=0;
	while(fabs(e2[Atom.size()-1]-e1[Atom.size()-1])> error ){
	    new_rho(rho,grad,Atom, alfa, N);
		vh.update(r,rho,t,h,Z,N);
		vxc->update(r,rho,grad,N);
	    for(int i=0;i<N;i++){
	        veff[i]=vn[i]+vh[i]+vxc->operator[](i)+vconf[i];
	     }
	     e1[Atom.size()-1]=e2[Atom.size()-1];


	     for(int i=0;i<Atom.size();i++){
	    	 Atom[i]->resolver(veff,r,W);
	    	 e2[i]=Atom[i]->energy();

	     }
	     iteraciones++;
	     if(fabs(e2[Atom.size()-1]-e1[Atom.size()-1])<0.0001){
	    	 alfa=0.5;
	     }
	}
	cout<<"iteraciones: "<<iteraciones<<endl;
	cout<<"Eh"<<"   "<<vh.energy()<<endl;
	cout<<"Exc"<<"  "<<vxc->energy()<<endl;

	for (int i=0;i<Atom.size();i++){
		cout<<Atom[i]->prin()<<angular[Atom[i]->ang()]<<"  "<<Atom[i]->energy()<<endl;

	}
	//***Se imprime  Rln de todos los orbitales en archivo radial.txt

	ofstream salida("radial2.txt");
	for(int i=0;i<N;i++){
		salida<<exp(t[i])/Z<<"   ";
		for(int j=0;j<Atom.size();j++){
			salida<<Atom[j]->operator()(i)<<"   ";
		}

		salida<<endl;
	}
	salida.close();

	/******Se libera spline y acc*********************/


	delete [] e1;
	delete [] e2;

}
double* Scf::Veff_noconf(){
	/*for(int i=0;i<N;i++){
		veff[i]=veff[i]-vconf[i];
	}*/
return veff;};
