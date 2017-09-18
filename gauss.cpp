/*
 * gauss.cpp
 *
 *  Created on: Apr 5, 2017
 *      Author: leandro
 */

#include "gauss.h"



gauss::gauss() {
	// TODO Auto-generated constructor stub
	n=20;

}
gauss::gauss(int N,double o) {
	n=N;
	t=gsl_integration_glfixed_table_alloc (n);
	a=o;
	c1=(2./log(3.))*acosh(1.+0.4/a);
	c2=-0.25;

}
gauss::~gauss() {
	// TODO Auto-generated destructor stub
}

void gauss::integrate2d(double b1,double b2,double c,double d, vector<Integrand*> S,vector <Integrand*> V,double*s,double*y){
     double x1[n];
     double x2[n];
     double w1[n];
     double w2[n];

     for(int i=0;i<n;i++){
    	 gsl_integration_glfixed_point (b1, b2, i, &x1[i], &w1[i], t);
    	 gsl_integration_glfixed_point (c, d, i, &x2[i], &w2[i], t);
     }

     for(int j=0;j<n;j++){
    	 for (int k=0;k<n;k++){
    		 double u=c1*atanh(x1[j]);
    		 double v=x2[k]+c2*sin(2*x2[k]);
             r1=a*(cosh(u)+cos(v));
             r2=a*(cosh(u)-cos(v));
             r=a*sinh(u)*sin(v);
             z=a*cosh(u)*cos(v);
             sin1=(r/r1);
             sin2=(r/r2);
             cos1=(z+a)/r1;
             cos2=(z-a)/r2;
    		 dV=sinh(u)*sin(v)*(c1/(1-x1[j]*x1[j]))*(1+2*c2*cos(2*x2[k]));
             for(int l=0;l<S.size();l++){
    			 if(S[l]!=NULL and V[l]!=NULL){
    				 s[l]=s[l]+a*w1[j]*w2[k]* S[l]->calculate(r1,r2,sin1,sin2,cos1,cos2)*dV;
    				 y[l]=y[l]+a*w1[j]*w2[k]* V[l]->calculate(r1,r2,sin1,sin2,cos1,cos2)*dV;
    			 }
    		 }

    	 }
     }


}

