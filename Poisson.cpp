/*
 * Poisson.cpp
 *
 *  Created on: Dec 27, 2016
 *      Author: leandro
 */

#include "Poisson.h"

Poisson::Poisson() {
	// TODO Auto-generated constructor stub

}

Poisson::~Poisson() {
	// TODO Auto-generated destructor stub
}

Poisson::Poisson(gsl_spline *spline,gsl_interp_accel *ac,int z){
 Z=z;
 n=spline;
 acc=ac;
}
void Poisson::operator ()( const state_type &y, state_type &f,const double t){
	double r=exp(t)/Z;
	f[0]=y[1];
	f[1]=-y[1]-4*pi*gsl_spline_eval (n, t, acc)*r*r;

};
