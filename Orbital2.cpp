/*
 * Orbital.cpp
 *
 *  Created on: Nov 14, 2016
 *      Author: leandro
 */
#include "Orbital.h"

Orbital::Orbital(int principal,int angular,double guess_e,double A) {
n=principal;
l=angular;
e=guess_e;
a=A;
/*for(int i=0;i<N;i++){
    Rl[i]=0;
}*/
}

Orbital::~Orbital() {
	// TODO Auto-generated destructor stub
}

void Orbital::inward(double h){

	ofstream archivo("radial2.txt");
	parametros p;
	p.e=e;
    p.l=l;
	int N=(int) (log(9)/h);
    gsl_odeiv2_system sys = { func, NULL,2,&p};

	  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,-h, 1e-8, 1e-8);
      double t=log(10)+estimate_clastp();
	  double y[2];
	  final(exp (t),y[0],y[1]);
	  int i, s;

	  for (i = 0; i < N; i++)
	    {
	      s = gsl_odeiv2_driver_apply_fixed_step (d, &t, -h, 1, y);

	      if (s != GSL_SUCCESS)
	        {
	          printf ("error: driver returned %d\n", s);
	          break;
	        }

	      archivo<<exp (t)<<"   "<<y[0]<<"   "<<y[1]<<endl;
	    }

	  gsl_odeiv2_driver_free (d);
      archivo.close();
};



void Orbital::outward(double h){
	ofstream archivo("radial.txt");
	parametros p;
	p.e=e;
    p.l=l;
    int N=(int) ((estimate_clastp()-t0)/h);
    gsl_odeiv2_system sys = { func, NULL,2,&p};

	  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,h, 1e-8, 1e-8);
      double t=t0;
	  double y[2];
	  inicial(exp(t),y[0],y[1]);
	  int i, s;

	  for (i = 0; i < N+1; i++)
	    {
	      s = gsl_odeiv2_driver_apply_fixed_step (d, &t, h, 1, y);

	      if (s != GSL_SUCCESS)
	        {
	          printf ("error: driver returned %d\n", s);
	          break;
	        }

	      archivo<<exp(t)<<"   "<<y[0]<<"   "<<y[1]<<endl;
	    }

	  gsl_odeiv2_driver_free (d);
      archivo.close();
};
double Orbital::f1(double t,double R,double dR,double veff,double dveff){
		double r=exp(t);
	    double M=1-al*al/4*(veff-e);
		return dR+(l*(l+1)+r*r*M*(veff-e))*R-al*al/(4*M)*dveff*r*(dR-R);
};




void Orbital::final(double r,double &R,double &dR){
		double vext=W-1/r;
	    //double k=sqrt(l*(l+1)/(r*r)+(vext-e));
	    //double dk=(1/2)*1/k*(-l/2*(l+1)/pow(r,1.5)+1/(r*r));
         double k=sqrt(vext-e);
         double dk=1/2*(1/k)*1/(r*r);
	    R=exp(-k*r);
        dR=	-r*k*R*dk;
};
void Orbital::inicial(double r,double &R,double &dR){
     //double gamma=(l*sqrt(l*l-al*al*Z*Z)+(l+1)*sqrt((l+1)*(l+1)-al*al*Z*Z))/(2*l+1);
	 //double gamma=sqrt(1-al*al*Z*Z);
     //R=pow(r,gamma);
     //dR=r*gamma*pow(r,gamma-1);
     R=a*r;
     dR=a*r;
}

/*void Orbital::density(double *rho){
	for (int i=0;i<N;i++){
		rho[i]=Rl[i]*Rl[i];
	}
}*/
double Orbital::estimate_clastp(){
	return log(-1/e);
};
void correct_e(){

};
int func (double t, const double y[], double f[], void *params)
{
  //(void)(t); /* avoid unused parameter warning */
  parametros *param = ( parametros *)params;
  int l=param->l;
  double e=param->e;
  double r=exp(t);
  double veff=-Z/r;
  //double dveff=Z/(r*r);
  //double M=1-al*al/4*(veff-e);
  f[0] = y[1];
  f[1] = y[1]+(l*(l+1)+r*r*(veff-e))*y[0];
  return GSL_SUCCESS;
}
