/*
 * splines.cpp
 *
 *  Created on: Nov 30, 2016
 *      Author: leandro
 */
#include <xc.h>

void xc(double *rho,double *exc, int N){
    xc_func_type func;
    xc_func_init(&func, 1, XC_UNPOLARIZED);
    xc_lda_exc(&func, N, rho, exc);
    xc_func_end(&func);
}


void hartree(double *rho,double *vh,int N,double *t){
	for (int i=0;i<N;i++){
		double inside=0;
		double outside=0;
		double h=0.001;
		for(int j=0;j<i;j++){
			inside=inside+rho[j]*h*exp(t[i]);

		}
        for (int k=i;k<N;k++){
        	outside=outside+rho[k]*h;
        }
	vh[i]=inside/exp(t[i])+outside;
	}
}
void vext(double *vh,double *exc,double *v,int N){
	for (int i=0;i<N;i++){
		v=vh[i]+exc[i];
	}
}
