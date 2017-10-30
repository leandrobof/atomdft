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

void gauss::integrate2d(Potential_spline &veff,Potential_spline &vconf,Orbital_spline **A,Orbital_spline **B, vector<Integrand*> S,vector <Integrand*> V,double*s,double*y){
     double x1[n];
     double x2[n];
     double w1[n];
     double w2[n];
     double cosh_u;
     double cos_v;
     double rmin=veff.min();
     double r1max=A[0]->max();
     double r2max=B[0]->max();

     for(int i=0;i<n;i++){
    	 gsl_integration_glfixed_point (0,1, i, &x1[i], &w1[i], t);
    	 gsl_integration_glfixed_point (0,pi, i, &x2[i], &w2[i], t);
     }

     for(int j=0;j<n;j++){
    	 for (int k=0;k<n;k++){
    		 double u=c1*atanh(x1[j]);
    		 double v=x2[k]+c2*sin(2*x2[k]);
    		 cosh_u=cosh(u);
    		 cos_v=cos(v);
    		 r1=a*(cosh_u+cos_v);
             r2=a*(cosh_u-cos_v);
             r=a*sinh(u)*sin(v);
             z=a*cosh_u*cos_v;
             sin1=(r/r1);
             sin2=(r/r2);
             cos1=(z+a)/r1;
             cos2=(z-a)/r2;
    		 dV=r*(c1/(1-x1[j]*x1[j]))*(1+2*c2*cos(2*x2[k]));                //r=a*sinh(u)*sin(v);

    		 if(r1>r1max or r2>r2max or r1<rmin or r2<rmin){
    			 veff1=0;
    			 vconf2=0;

    		 }
    		 else{
    			 veff1=veff(r1);
    			 vconf2=vconf(r2);
    		 }
    		 for(int i=0;i<3;i++){
    			 if(r1>r1max or r2>r2max){
    				 R1[i]=0;
    				 R2[i]=0;
    			 }
    			 else{
    				 R1[i]=(*A[i])(r1);
    				 R2[i]=(*B[i])(r2);
    			 }
    		 }
    		 for(int l=0;l<S.size();l++){
    			 if(S[l]!=NULL and V[l]!=NULL){
    				 s[l]=s[l]+w1[j]*w2[k]* S[l]->calculate(R1,R2,veff1,vconf2,sin1,sin2,cos1,cos2)*dV;
    				 y[l]=y[l]+w1[j]*w2[k]* V[l]->calculate(R1,R2,veff1,vconf2,sin1,sin2,cos1,cos2)*dV;
    			 }
    		 }

    	 }
     }


}

