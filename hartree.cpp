/*
 * hartree.cpp
 *
 *  Created on: Mar 8, 2017
 *      Author: leandro
 */
#include "hartree.h"
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

 void Hartree::update(double *r,double *rho,double *t,double h,int Z,int N){


/*Se crea un spline de rho para poder evaluar con rungekutta4 */

gsl_spline *n;
gsl_interp_accel *acc;
n=gsl_spline_alloc (gsl_interp_cspline, N);
acc = gsl_interp_accel_alloc ();
gsl_spline_init (n, t, rho, N);
Poisson P(n,acc,Z);
state_type y(2);

/*Condiciones iniciales,hay que calcular la integral de  r*rho.
Se crea un array r[]*rho[] se interpola con spline cubico y se integra para obtener y[0].*/

double rrho[N];
for (int i=0;i<N;i++){
	rrho[i]=r[i]*rho[i];
}
gsl_spline *rn;
gsl_interp_accel *acc2;
rn=gsl_spline_alloc (gsl_interp_cspline, N);
acc2=gsl_interp_accel_alloc ();
gsl_spline_init (rn, r, rrho, N);

y[0]=4*pi*gsl_spline_eval_integ (rn, r[0], r[N-1], acc2);
y[1]=0;
v[0]=y[0];
gsl_spline_free (rn);
gsl_interp_accel_free (acc2);
//*******************************************************


runge_kutta4 < state_type >   rk;
//Se resuelve usando integracion adaptativa,ver si es correcto.Error=10^-12.
for(int i=0;i<N-2;i++){                                     //N-1 array arreglar.

	rk.do_step( P , y , t[i] , h);
    v[i+1]=y[0];

}

gsl_spline_free (n);
gsl_interp_accel_free (acc);
//**********************

double I[N];
for(int i=0;i<N;i++){
	I[i]=v[i]*rho[i]*r[i]*r[i];
}
gsl_interp_accel *acc3 = gsl_interp_accel_alloc ();
gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
gsl_spline_init (spline, r, I,N );

Eh=2*pi*gsl_spline_eval_integ (spline, r[0],r[N-1], acc3);

gsl_spline_free (spline);
gsl_interp_accel_free (acc3);

}
Hartree::Hartree(double *r,double *rho,double *t,double h, int Z,int N){
	v=new double [N];
	update(r,rho,t,h,Z,N);
}
Hartree::~Hartree(){
	delete [] v;
}
double Hartree::operator[](int i){return v[i];}
