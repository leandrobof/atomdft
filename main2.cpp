/*
 * main.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: leandro
 */

#include "Orbital.h"
#include "splines.h"
using namespace std;


int main(int argc,char *argv[]){
int num=atoi(argv[1]);

int N=(tinf-t0)/h;
double t[N];
double veff[N];


double density[N];
double vxc[N],vh[N];
double energy0,energy1;

for (int i=0;i<N;i++){
	veff[i]=0.;
    density[i]=0.;
	t[i]=t0+i*h;

}
gsl_interp_accel *acc = gsl_interp_accel_alloc ();
gsl_spline *v = gsl_spline_alloc (gsl_interp_cspline, N);
Orbital H1s(num,0,1,veff,t);
gsl_spline_init (v, t, veff, N);
while(fabs(energy1-energy0)>0.00000001){
    H1s.resolver(t,v,acc);

    for (int i=0;i<N;i++){
	    density[i]=H1s[i];
    }
    xc(density,vxc,N);
    hartree(density,vh,N,t);
    vext(vh,vxc,veff,N);
    gsl_interp_accel_reset (acc);
    gsl_spline_init (v, t, veff, N);
}
H1s.print();
gsl_spline_free (v);
gsl_interp_accel_free (acc);
return 0;
}




