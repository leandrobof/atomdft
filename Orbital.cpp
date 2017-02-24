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


Orbital::Orbital(int principal,int angular,double num_ocupacion,double *veff,double *t) {
n=principal;
l=angular;
e=-1;
noc=num_ocupacion;
a=1;
tk=0;
Nt=(tinf-t0)/h;
Rl=new double [Nt];
rho=new double [Nt];
acc = gsl_interp_accel_alloc ();
spline = gsl_spline_alloc (gsl_interp_cspline, Nt);
gsl_spline_init (spline, t, veff, Nt);

}

Orbital::~Orbital() {
	// TODO Auto-generated destructor stub
	delete [] Rl;
    delete [] rho;

}

void Orbital::inward(double *t,gsl_spline *spline,gsl_interp_accel *acc){

	//ofstream archivo("radial2.txt");
      Func Sch(e,l,spline,acc);
	  state_type  y(2);
	  int i=Nt-1;
	  final(exp (t[i]),y);

	  runge_kutta4<state_type> rk;
	  while (t[i] > tk)
	    {

	      rk.do_step( Sch, y , t[i] , -h );


	    //archivo<<exp (t)<<"   "<<y[0]<<"   "<<y[1]<<endl;
	    i--;
	    }

	  ci=y[0];dRi=y[1];

      //archivo.close();
};



void Orbital::outward(double *t,gsl_spline *spline,gsl_interp_accel *acc){
	//ofstream archivo("radial.txt");
	nodos=0;

	  Func Sch(e,l,spline,acc);
	  state_type y(2);
	  double ya;
      inicial(exp(t[0]),y);
	  int i=0;
      double x0,x1;
      runge_kutta4< state_type > rk;

      double r=exp(t[i]);
      x0=l*(l+1)/(r*r)-Z/r+gsl_spline_eval (spline, t[i], acc)-e;
      x1=x0;
      while (x0*x1>0 )
	    {
	      x0=x1;
    	  ya=y[0];

	        rk.do_step( Sch , y , t[i] , h );

          if (ya*y[0]<0){
        	  nodos++;
          }
          i++;
          r=exp(t[i]);
          x1=l*(l+1)/(r*r)-Z/r+gsl_spline_eval (spline, t[i], acc)-e;
          //cout<<exp(t[i])<<"   "<<y[0]<<"   "<<y[1]<<endl;

	    }
	  co=y[0];dRo=y[1];
	  tk=t[i];


	  //archivo.close();
};
double Orbital::f1(double t,double R,double dR,double veff,double dveff){
		double r=exp(t);
	    double M=1-al*al/4*(veff-e);
		return dR+(l*(l+1)+r*r*M*(veff-e))*R-al*al/(4*M)*dveff*r*(dR-R);
};




void Orbital::final(double r,state_type &R){
		double vext=W-1/r;
	    //double k=sqrt(l*(l+1)/(r*r)+(vext-e));
	    //double dk=(1/2)*1/k*(-l/2*(l+1)/pow(r,1.5)+1/(r*r));
         double k=sqrt(vext-e);
         double dk=1/2*(1/k)*1/(r*r);
	    R[0]=exp(-k*r);
        R[1]=	-r*k*R[0]*dk;
};
void Orbital::inicial(double r,state_type &R){
     //double gamma=(l*sqrt(l*l-al*al*Z*Z)+(l+1)*sqrt((l+1)*(l+1)-al*al*Z*Z))/(2*l+1);
	 //double gamma=sqrt(1-al*al*Z*Z);
     //R=pow(r,gamma);
     //dR=r*gamma*pow(r,gamma-1);
     R[0]=a*r;
     R[1]=a*r;
}

/*void Orbital::density(double *rho){
	for (int i=0;i<N;i++){
		rho[i]=Rl[i]*Rl[i];
	}
}*/


double Orbital::correct_e(){
    double de=(dRi-a*dRo);
	//e=e-de;
    a=1;
    return de;
};

void Orbital::estimate_a(){
	a=ci/co;
};


void Orbital::resolver(double *t,gsl_spline *spline,gsl_interp_accel *acc){


	double error=0.00000001;
	double sup=0.00;
	double inf=-Z*Z/(1.3*(n*n));
	double dif1=1;
	double dif2=1;
	e=inf;
	while(nodos !=n-1  ){
	      outward(t,spline,acc);
	      if (nodos> n-1){
        	  sup=e;
        	  e=(sup+inf)/2;

          }
          else if(nodos<n-1){
        	  inf=e;
        	  e=(sup+inf)/2;
          }


	}
	inward(t,spline,acc);
	estimate_a();
	double a=correct_e();
	double b;
	if(a > 0){
		sup=e;

		do{
		 e=(sup+inf)/2;
	     outward(t,spline,acc);
	     inf=e;

		}while(nodos !=n-1  );

		inward(t,spline,acc);
		estimate_a();
		b=correct_e();
		if (a*b>0){
            inf=inf-(sup-inf)/2;
		}

	}
	else{
		inf=e;

		do{
			e=(sup+inf)/2;
			outward(t,spline,acc);
			sup=e;

	    }while(nodos !=n-1 );
	  	inward(t,spline,acc);
		estimate_a();
		b=correct_e();
		if (a*b>0){
	        sup=sup+(sup-inf)/2;
			}
	}

    while (fabs(inf-sup)>error){
    	e=sup;
    	outward(t,spline,acc);
	    inward(t,spline,acc);
	    estimate_a();
	    dif1=correct_e();
	    e=(sup+inf)/2;
	    outward(t,spline,acc);
	    inward(t,spline,acc);
	    estimate_a();
	    dif2=correct_e();
	    if (dif1*dif2<0){
	    	inf=e;

	    }
	    else{
	    	sup=e;
	    }
	    //cout<<"energia:"<<H1s.get_energy()<<"      "<<"error:"<<dif2<<endl;

	}
	cout<<"energia: "<<e<<"      "<<"error:"<<dif2<<endl;

	outward(t,spline,acc);
	inward(t,spline,acc);
    estimate_a();
    radial(t,spline,acc);
    dens(t);


}

/*void Orbital::print_inward(double h){

	//ofstream archivo("radial2.txt");
	parametros p;
	p.e=e;
    p.l=l;

    gsl_odeiv2_system sys = { func, NULL,2,&p};

	  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,-h, 1e-8, 1e-8);

      double t=tinf;
	  double y[2];
	  final(exp (t),y[0],y[1]);
	  int i=Nt-2, s;

	  while (t > tk)
	    {
	      s = gsl_odeiv2_driver_apply_fixed_step (d, &t, -h, 1, y);


         Rl[i]=y[0];
         i--;
	    //archivo<<exp (t)<<"   "<<y[0]<<"   "<<y[1]<<endl;
	    }

	  gsl_odeiv2_driver_free (d);
      //archivo.close();
};

*/

void Orbital::radial(double *t,gsl_spline *spline,gsl_interp_accel *acc){
	//ofstream archivo("radial.txt");
	nodos=0;

	  Func Sch(e,l,spline,acc);
	  state_type y(2);

      inicial(exp(t[0]),y);
	  int i=0;
      double x0,x1;
      runge_kutta4< state_type > rk;

      double r=exp(t[i]);
      x0=l*(l+1)/(r*r)-Z/r+gsl_spline_eval (spline, t[i], acc)-e;
      x1=x0;
      while (x0*x1>0 ){
	      x0=x1;
    	  rk.do_step( Sch , y , t[i] , h );
          i++;
          r=exp(t[i]);
          x1=l*(l+1)/(r*r)-Z/r+gsl_spline_eval (spline, t[i], acc)-e;
          Rl[i]=y[0];

	  };

      i=Nt-1;
      final(exp (t[i]),y);

   	  while (t[i] > tk){
          rk.do_step( Sch, y , t[i] , -h );
     	  Rl[i]=y[0];
      	  i--;
      }


};





	  //archivo.close();

void Orbital::print(){
	ofstream archivo("radial.txt");
	double t=t0;
	double h=0.001;
	for (int i=0;i<Nt;i++){
		archivo<<exp(t)<<"   "<<Rl[i]<<"   "<<rho[i]<<endl;
	    t=t+h;
	}
	archivo.close();
}
void Orbital::dens(double *t){
	for (int i=0;i<Nt;i++){
		rho[i]=Rl[i]*Rl[i]/(exp(t[i]));
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, Nt);

	    gsl_spline_init (spline, t, rho, Nt);

	    //yi = gsl_spline_eval (spline, xi, acc);
	    double norm= gsl_spline_eval_integ (spline, t[0], t[Nt-1], acc);

        for(int i=0;i<Nt;i++){
        	rho[i]=rho[i]/(norm*exp(t[i]));
        };

        gsl_spline_free (spline);
	    gsl_interp_accel_free (acc);

}


double Orbital::operator[](int i){
	return rho[i]/(4*pi);
};
