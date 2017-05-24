/*
 * Orbital.cpp
 *
 *  Created on: Nov 14, 2016
 *      Author: leandro
 */
#include "Orbital.h"
#include "Func.h"
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;


Orbital::Orbital(int principal,int angular,double num_ocupacion,int z,double e_orbital,double *r,int N) {


n=principal;
l=angular;
e=e_orbital;
Z=z;
Nt=N;
t=r;
h=t[1]-t[0];
tinf=t[N];
t0=t[0];
tk=0;
max=Nt;
a=1;
noc=num_ocupacion;
Rl=new double [Nt];
rho=new double [Nt];

}

Orbital::~Orbital() {
	// TODO Auto-generated destructor stub
	delete [] Rl;
    delete [] rho;

}

void Orbital::inward(gsl_spline *spline,gsl_interp_accel *acc){


      Func Sch(e,l,Z,spline,acc);
	  state_type  y(2);
	  y[0]=0;
	  y[1]=0;
	  int i=Nt-1;
      while(y[0]<1.0e-15){
    	  final(exp (t[i])/Z,y);
          i--;
      }
      max=i;

	  
	adams_bashforth_moulton< 5 , state_type > abm;
         Rl[i]=y[0];
	 while (t[i] > tk)
	    {

	      abm.do_step( Sch, y , t[i] , -h );



	    i--;
	    Rl[i]=y[0];
	    }

	  ci=y[0];dRi=y[1];
      for(int k=max;k<Nt;k++){
    	  Rl[k]=0.;
      }

};



void Orbital::outward(gsl_spline *spline,gsl_interp_accel *acc){


	nodos=0;

	  Func Sch(e,l,Z,spline,acc);
	  state_type y(2);
	  double ya;
      inicial(exp(t[0])/Z,y);
	  int i=0;
      double x0,x1;
      adams_bashforth_moulton< 5 , state_type > abm;
      Rl[i]=y[0];
      x0=gsl_spline_eval (spline, t[i], acc)-e;
      x1=x0;
      while (x0*x1>0 )
	    {
	      x0=x1;
    	  ya=y[0];

	        abm.do_step( Sch , y , t[i] , h );

          if (ya*y[0]<0){
        	  nodos++;
          }
          i++;
          x1=gsl_spline_eval (spline, t[i], acc)-e;

        Rl[i]=y[0];
	    }

	  tk=t[i];
	  co=y[0];dRo=y[1];


};


void Orbital::final(double r,state_type &R){
		/* Asumiendo V->0 cuando r->0 ,segun
		  O Čertík ,J E. Pask ,J Vackář Computer Physics Communications 184 (2013) 1777–1791*/

        double k=sqrt(-2*e);

	    R[0]=a*exp(-k*r);
        R[1]=	-k*R[0]*r;               //(df/dt)=(df/dr)*r
};

void Orbital::inicial(double r,state_type &R){
	 /*O Čertík ,J E. Pask ,J Vackář Computer Physics Communications 184 (2013) 1777–1791*/

	 R[0]=pow(r,l+1);
     R[1]=(l+1)*pow(r,l+1);          //(df/dt)=(df/dr)*r por eso l+1.

}




double Orbital::correct_e(){
	//correccion utilizando teoria de perturbaciones.
	/*O Čertík ,J E. Pask ,J Vackář Computer Physics Communications 184 (2013) 1777–1791*/

	double de=co*(dRo-a*dRi)*Z/exp(tk);   // Si lo divido por /exp(tk) converge muy lentamente. Ver cual es la forma correcta.


    int i=0;
    while(t[i]<tk){
    		i++;
    }

    double P1[i];
    double t1[i];
    double P2[max-i];
    double t2[max-i];

    for(int j=0;j<i;j++){
    	P1[j]=Rl[j]*Rl[j]*exp(t[j])/Z;
        t1[j]=t[j];
    }
    for(int k=i;k<max;k++){
        	P2[k-i]=Rl[k]*Rl[k]*exp(t[k])/Z;
            t2[k-i]=t[k];
    }

    gsl_spline *p1 = gsl_spline_alloc (gsl_interp_cspline, i);
    gsl_spline *p2 = gsl_spline_alloc (gsl_interp_cspline, max-i);
    gsl_interp_accel *acc1 = gsl_interp_accel_alloc ();
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();


    gsl_spline_init (p1, t1, P1, i);
    gsl_spline_init (p2, t2, P2, max-i);

    double norm=gsl_spline_eval_integ (p1, t1[0], t[i-1], acc1)+gsl_spline_eval_integ (p2, t2[0], t2[max-1-i], acc2);

    gsl_spline_free (p1);
    gsl_spline_free (p2);
    gsl_interp_accel_free (acc1);
    gsl_interp_accel_free (acc2);

    a=1;
    return de/(2*norm);
};


void Orbital::resolver(gsl_spline *spline,gsl_interp_accel *acc){
    double e_old=0.;
    double error=0.000000001;
    double de=0;
    double inf=-Z*Z/2.;
    double sup=0.;
    outward(spline,acc);
    int nod=0;
    int perturbaciones=0;
    while(nodos !=n-l-1  ){
    	      if (nodos> n-l-1){
            	  sup=e;
            	  e=(sup+inf)/2;

              }
              else if(nodos<n-l-1){
            	  inf=e;
            	  e=(sup+inf)/2;

              }
    outward(spline,acc);
    nod++;
    }

   while(fabs(e-e_old)>error ){

    e_old=e;
	inward(spline,acc);
    estimate_a();
    radial();
    de=correct_e();
    perturbaciones++;
    e=e+de;
    outward(spline,acc);
    while(nodos !=n-l-1){
    	de=de/2;
    	e=e_old+de;
        outward(spline,acc);

    };
    }



cout<<e<<endl;



cout<<"evaluciones nodo: "<<nod<<"    "<<"evaluciones pertubaciones: "<<perturbaciones<<endl;
dens();
}


void Orbital::radial(){

      int i=Nt-1;

      while (t[i] >= tk){
    	  Rl[i]=a*Rl[i];
          i--;

      }


};

void Orbital::estimate_a(){a=co/ci;};





void Orbital::print(){
	ofstream archivo("radial.txt");


	for (int i=0;i<Nt;i++){
		archivo<<exp(t[i])/Z<<"   "<<Rl[i]<<"   "<<rho[i]<<endl;

	}
	archivo.close();
}
void Orbital::dens(){
	for (int i=0;i<Nt;i++){
		rho[i]=Rl[i]*Rl[i]*exp(t[i])/Z;
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, Nt);

	    gsl_spline_init (spline, t, rho, Nt);


	    double norm=1/sqrt(gsl_spline_eval_integ (spline, t[0], t[Nt-1], acc));

        for(int i=0;i<Nt;i++){
        	Rl[i]=norm*Rl[i];
        	rho[i]=Rl[i]*Rl[i]/(4*pi*exp(t[i])*exp(t[i]));

        };

        gsl_spline_free (spline);
	    gsl_interp_accel_free (acc);

}



double Orbital::operator[](int i){
	return rho[i];
};
double Orbital::operator()(int i){
	return Rl[i]*Rl[i];
};
Orbital_rel::Orbital_rel(int principal,int angular,double num_ocupacion,int s,int z,double energy,double *r,int N) : Orbital(principal,angular,num_ocupacion,z,energy,r,N){
	if (l == 0){k=-1;}
	else if(s == 1 ){k=-l-1;}
	else{k=l;}
    Ql=new double [Nt];
} ;

Orbital_rel::~Orbital_rel(){
	delete [] Ql;
};

void Orbital_rel::inicial(double r,state_type &R){
	double b=sqrt(k*k-(Z/c)*(Z/c));
	R[0]=pow(r,b);
	R[1]=R[0]*c*(b+k)/Z;
}
void Orbital_rel::final(double r,state_type &R){
	double lamb=sqrt(-2*e-e*e/(c*c));
    R[0]=exp(-lamb*r);
    R[1]=-sqrt(-(e/(e+2*c*c)))*R[0];
}
void Orbital_rel::outward(gsl_spline *spline,gsl_interp_accel *acc){


	nodos=0;

	  Dirac Ec(e,l,k,Z,spline,acc);
	  state_type y(2);
	  double ya;
	  int i=0;
	  inicial(exp (t[i])/Z,y);
	  while(y[0]<1.0e-16){
		      Rl[i]=0; Ql[i]=0;
		      inicial(exp (t[i])/Z,y);
	      	  i++;

	  }
	  double x0,x1;
      adams_bashforth_moulton< 5 , state_type > abm;
      Rl[i]=y[0]; Ql[i]=y[1];
      x0=gsl_spline_eval (spline, t[i], acc)-e; //Ver cual es la condicion correcta.
      x1=x0;
      while (x0*x1>0 )
	    {
	      x0=x1;
    	  ya=y[0];

	        abm.do_step( Ec , y , t[i] , h );

          if (ya*y[0]<0){
        	  nodos++;
          }
          i++;
          x1=gsl_spline_eval (spline, t[i], acc)-e;

        Rl[i]=y[0]; Ql[i]=y[1];
	    }

	  tk=t[i];
	  co=y[0];dRo=y[1];


};
void Orbital_rel::inward(gsl_spline *spline,gsl_interp_accel *acc){


	  Dirac Ec(e,l,k,Z,spline,acc);
	  state_type  y(2);
	  y[0]=0;
	  y[1]=0;
	  int i=Nt-1;
      while(y[0]<1.0e-12){
    	  final(exp (t[i])/Z,y);
          i--;
      }
      max=i;


	adams_bashforth_moulton< 5 , state_type > abm;
         Rl[i]=y[0]; Ql[i]=y[1];
	 while (t[i] > tk)
	    {

	      abm.do_step( Ec, y , t[i] , -h );



	    i--;
	    Rl[i]=y[0]; Ql[i]=y[1];
	    }

	  ci=y[0];dRi=y[1];
      for(int k=max;k<Nt;k++){
    	  Rl[k]=0.; Ql[k]=0;
      }

};
double Orbital_rel::correct_e(){
	{
		//correccion utilizando teoria de perturbaciones.
		/*O Čertík ,J E. Pask ,J Vackář Computer Physics Communications 184 (2013) 1777–1791*/

		double de=co*(dRo-a*dRi);   // En la ecuacion Dirac no va con exp(tk)/Z.


	    int i=0;
	    while(t[i]<tk){
	    		i++;
	    }

	    double P1[i];
	    double t1[i];
	    double P2[max-i];
	    double t2[max-i];

	    for(int j=0;j<i;j++){
	    	P1[j]=(Rl[j]*Rl[j]+Ql[j]*Ql[j])*exp(t[j])/Z;
	        t1[j]=t[j];
	    }
	    for(int k=i;k<max;k++){
	        	P2[k-i]=(Rl[k]*Rl[k]+Ql[k]*Ql[k])*exp(t[k])/Z;
	        	t2[k-i]=t[k];
	    }

	    gsl_spline *p1 = gsl_spline_alloc (gsl_interp_cspline, i);
	    gsl_spline *p2 = gsl_spline_alloc (gsl_interp_cspline, max-i);
	    gsl_interp_accel *accp1 = gsl_interp_accel_alloc ();
	    gsl_interp_accel *accp2 = gsl_interp_accel_alloc ();


	    gsl_spline_init (p1, t1,P1 , i);
	    gsl_spline_init (p2, t2, P2, max-i);


	    double norm=gsl_spline_eval_integ (p1, t1[0], t[i-1], accp1)+gsl_spline_eval_integ (p2, t2[0], t2[max-1-i], accp2);

	    gsl_spline_free (p1);
	    gsl_spline_free (p2);
	    gsl_interp_accel_free (accp1);
	    gsl_interp_accel_free (accp2);

	    a=1;
	    return c*de/norm;
	};
};
void Orbital_rel::radial(){

      int i=Nt-1;

      while (t[i] >= tk){
    	  Rl[i]=a*Rl[i];
          Ql[i]=a*Ql[i];
    	  i--;

      }


};
void Orbital_rel::dens(){
	for (int i=0;i<Nt;i++){
		rho[i]=(Rl[i]*Rl[i]+Ql[i]*Ql[i])*exp(t[i])/Z;
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, Nt);

	    gsl_spline_init (spline, t, rho, Nt);


	    double norm=1/sqrt(gsl_spline_eval_integ (spline, t[0], t[Nt-1], acc));

        for(int i=0;i<Nt;i++){
        	Rl[i]=norm*Rl[i];
        	Ql[i]=norm*Ql[i];
        	rho[i]=(Rl[i]*Rl[i]+Ql[i]*Ql[i])/(4*pi*exp(t[i])*exp(t[i]));

        };

        gsl_spline_free (spline);
	    gsl_interp_accel_free (acc);

}
void Orbital_rel::print(){
	ofstream archivo("radial.txt");


	for (int i=0;i<Nt;i++){
		archivo<<exp(t[i])/Z<<"   "<<Rl[i]<<"   "<<Ql[i]<<"   "<<rho[i]<<endl;
	}
	archivo.close();

}
