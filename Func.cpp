/*
 * Func.cpp
 *
 *  Created on: Dec 7, 2016
 *      Author: leandro
 */

#include "Func.h"

Func::Func(double a,int b,gsl_spline *v,gsl_interp_accel *ac) {
	// TODO Auto-generated constructor stub
	e=a;
	l=b;
	spline=v;
	acc=ac;
}

Func::~Func() {
	// TODO Auto-generated destructor stub
}


void Func::operator ()( const state_type &y, state_type &f,const double t){
  double r=exp(t);
  double veff=-Z/r+gsl_spline_eval (spline, t, acc);

  //double dveff=Z/(r*r);
  //double M=1-al*al/4*(veff-e);
  f[0] = y[1];
  f[1] = y[1]+(l*(l+1)+r*r*(veff-e))*y[0];

};
