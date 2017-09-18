/*
 * Func.h
 *
 *  Created on: Dec 7, 2016
 *      Author: leandro
 */
#include "globals.h"


#ifndef FUNC_H_
#define FUNC_H_

class Func {
protected:
	double e;
	int l;
	int Z;
	double *v;
	double *t;
    double *r;
public:
	Func(double a,int b,int z,double *veff,double *t,double *r);
	virtual void operator ()( double  *y, double *f,  int i);
	virtual ~Func();

};
class Dirac : public Func{
private:
	int k;
public:
	Dirac(double a,int b,int c,int z,double *veff,double *t,double *r);
	virtual void operator()( double *y, double *f, int i);
	virtual ~Dirac();
};

#endif /* FUNC_H_ */
