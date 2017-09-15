/*
 * xc_potential.h
 *
 *  Created on: Mar 27, 2017
 *      Author: leandro
 */

#ifndef XC_POTENTIAL_H_
#define XC_POTENTIAL_H_
#include <xc.h>
#include "globals.h"

class Xc{
protected:
	double *ec;
	double *ex;
	double *vx;
	double *vc;
	double *vxc;

    double Exc;
    xc_func_type correlacion;
    xc_func_type exchange;
public:
    Xc(int N);
    virtual ~Xc();
    virtual void update(double *r,double *rho,double *grad,int N)=0;
    double operator[](int);
    double energy(){return Exc;};
};
class Xc_lda :public Xc{
public:
	Xc_lda(int N,bool);
	virtual ~Xc_lda();
	virtual void update(double *r,double *rho,double *grad,int N);

};
class Xc_gga :public Xc{
private:
	double *vcsigma;
    double *vxsigma;
	double *sigma;
public:
	Xc_gga(int N);
	virtual ~Xc_gga();
	virtual void update(double *r,double *rho,double *grad,int N);

};

#endif /* XC_POTENTIAL_H_ */
