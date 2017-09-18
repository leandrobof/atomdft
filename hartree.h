/*
 * hartree.h
 *
 *  Created on: Mar 9, 2017
 *      Author: leandro
 */

#ifndef HARTREE_H_
#define HARTREE_H_



#include "Poisson.h"
#include "globals.h"



class Hartree{
private:
	double *v;
	double Eh;
public:
    Hartree(double *r,double *rho,double *t,double h,int Z,int N);
    ~Hartree();
    double energy(){return Eh;};
    void update(double *,double *,double *t,double h,int Z,int N);
    double operator[](int);
};

#endif /* HARTREE_H_ */
