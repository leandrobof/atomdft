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


public:
	gauss();
	gauss(int N,double o);
	virtual ~gauss();
    void integrate2d(double b1,double b2,double c,double d, vector<Integrand*> S,vector <Integrand*> V,double*s,double*v);
    void update_a(double o){a=o;};
};

#endif /* GAUSS_H_ */
