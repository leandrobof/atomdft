/*
 * Func.cpp
 *
 *  Created on: Dec 7, 2016
 *      Author: leandro
 */

#include "Func.h"

Func::Func(double a,int b,int z,double *veff,double *j,double *R) {
	// TODO Auto-generated constructor stub
	e=a;
	l=b;

	v=veff;

    r=R;

}

Func::~Func() {

}


 void Func::operator ()( double *y, double *f,  int i){



  //double dveff=Z/(r*r);
  //double M=1-al*al/4*(veff-e);
  f[0] = y[1];
  f[1] = (l*(l+1)+2*r[i]*r[i]*(v[i]-e))*y[0]+y[1];

};
Dirac::Dirac(double a,int b,int d,int z,double *veff,double *c,double  *R) : Func(a,b,z,veff,c,R){
	k=d;
};
Dirac::~Dirac(){};
void Dirac::operator ()(  double *y, double *f, int i){



  //double dveff=Z/(r*r);
  //double M=1-al*al/4*(veff-e);
  f[0] = -k*y[0]+r[i]*((e-v[i])/c+2*c)*y[1];
  f[1] = -r[i]*((e-v[i])/c)*y[0]+k*y[1] ;

};
