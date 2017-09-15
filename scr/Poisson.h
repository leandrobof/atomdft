/*
 * Poisson.h
 *
 *  Created on: Dec 27, 2016
 *      Author: leandro
 */

#ifndef POISSON_H_
#define POISSON_H_
#include "globals.h"
class Poisson {
public:
	Poisson();
	Poisson(gsl_spline *n,gsl_interp_accel *acc,int z);
	virtual ~Poisson();
	void operator ()( const state_type &y, state_type &f,const double r);
private:
	int Z;
	gsl_spline *n;
	gsl_interp_accel *acc;
};



#endif /* POISSON_H_ */
