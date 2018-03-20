/*
 * gauss.h
 *
 *  Created on: Apr 5, 2017
 *      Author: leandro
 */

#ifndef GAUSS_H_
#define GAUSS_H_
#include "Integrand.h"
#include <gsl/gsl_integration.h>

using namespace std;
class gauss {
private:
    int n;
    gsl_integration_glfixed_table *t;
    double a;
    double c1;
    double c2;
    double r1;
    double r2;
    double r;
    double z;
    double sin1;
    double sin2;
    double cos1;
    double cos2;
    double dV;
    double veff1;
    double veff2;
    double veff02;
    double DVb;
    double vconf2;
    double R1[3];
    double R2[3];

public:
	gauss();
	gauss(int N,double o);
	virtual ~gauss();
    void integrate2d(Potential_spline &veffa,Potential_spline &veffb,Potential_spline &vconfb,Potential_spline &veff0b,Orbital_spline **A,Orbital_spline **B, vector<Integrand*> S,vector <Integrand*> V,double*s,double*y);
    void update_a(double o){a=o;
    c1=(2./log(3.))*acosh(1.+0.4/a);};
};

#endif /* GAUSS_H_ */
