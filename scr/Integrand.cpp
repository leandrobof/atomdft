/*
 * Integrand.cpp
 *
 *  Created on: Apr 5, 2017
 *      Author: leandro
 */

#include "Integrand.h"

Integrand::Integrand() {
	// TODO Auto-generated constructor stub

}

Integrand::~Integrand() {
	// TODO Auto-generated destructor stub
}
/*void Integrand::update_a(double o){a=o;c1=(2./log(3.))*acosh(1.+0.4/a);}

Gaussianas::Gaussianas(double y,double z,double o) : Integrand(o){
    z1=y;
    z2=z;
    N1=sqrt(pow((2*z1/pi),3./2.));
    N2=sqrt(pow((2*z2/pi),3./2.));
}
inline double Integrand::r1(double u,double v){
	return a*(cosh(u)+cos(v));
}
inline double Integrand::r2(double u,double v){
	return a*(cosh(u)-cos(v));
}

double Gaussianas::gaussian1(double x){
	return N1*exp(-z1*x*x);
}

double Gaussianas::gaussian2(double x){
	return N2*exp(-z2*x*x);
}

double Gaussianas::real(double Rab){
	double p=z1+z2;
    double u=z1*z2/(z1+z2);

    return N1*N2*pow((pi/p),3/2.)*exp(-u*Rab*Rab);
};

double Gaussianas::calculate(double &r1,double &r2,double &sin1,double &sin2,double &cos1,double &cos2){
	double u=c1*atanh(x1);
	double v=x2+c2*sin(2*x2);
    double I;
    if(r1(u,v)>45 or r2(u,v)>45){
    	I=0;
    }
    else {
  	    I=r1(u,v)*r2(u,v)*(gaussian1(r1(u,v))*gaussian2(r2(u,v)))*sinh(u)*sin(v)*(c1/(1-x1*x1))*(1+2*c2*cos(2*x2));
    }
    return 2*pi*a*I;
}

Density::Density(gsl_spline *a1,gsl_interp_accel *b1,gsl_spline *a2,gsl_interp_accel *b2,double o,double max1,double max2): Integrand(o){
	rho1=a1;
	rho2=a2;
	accrho1=b1;
	accrho2=b2;
    r1max=max1;
    r2max=max2;

}
double Density::calculate(double &r1,double &r2,double &sin1,double &sin2,double &cos1,double &cos2){
	    double u=c1*atanh(x1);
		double v=x2+c2*sin(2*x2);
	    double I1,I2;
	    if(r1(u,v)>r1max ){
	    	I1=0;
	    }
	    else {
	  	    I1=r2(u,v)*r1(u,v)*(gsl_spline_eval(rho1, r1(u,v), accrho1))*sinh(u)*sin(v)*(c1/(1-x1*x1))*(1+2*c2*cos(2*x2));
	    }

	    if(r2(u,v)>r2max){
	    	I2=0;
	    }
	    else{
	    	I2=r1(u,v)*r2(u,v)*(gsl_spline_eval(rho2,r2(u,v),accrho2))*sinh(u)*sin(v)*(c1/(1-x1*x1))*(1+2*c2*cos(2*x2));
	    }
	    return 2*pi*a*(I1+I2);
};
*/
S::S() : Integrand(){
      l1=0;
      l2=0;

}

double S::calculate(double *R1,double *R2,double v1,double vconf,double &sin1,double &sin2,double &cos1,double &cos2){



	return	R1[l1]*R2[l2]*f(sin1,sin2,cos1,cos2);



};
S_ss::S_ss() : S(){l1=0;l2=0;};
S_sp::S_sp() : S(){l1=0;l2=1;};
S_pp_sig::S_pp_sig() : S(){l1=1;l2=1;};
S_pp_pi::S_pp_pi() : S(){l1=1;l2=1;};
S_sd::S_sd() : S(){l1=0;l2=2;};
S_pd_pi::S_pd_pi() : S(){l1=1;l2=2;};
S_pd_sig::S_pd_sig() : S(){l1=1;l2=2;};
S_dd_del::S_dd_del() : S(){l1=2;l2=2;};
S_dd_pi::S_dd_pi() : S(){l1=2;l2=2;};
S_dd_sig::S_dd_sig() : S(){l1=2;l2=2;};

V::V(double ener) : S(){
	e=ener;


};
//Mejorar la evaluacion de integrales.Las casos menores a rmin no se evaluan actualmente.
double V::calculate(double *R1,double *R2,double v1,double vconf,double &sin1,double &sin2,double &cos1,double &cos2){


	return  R1[l1]*R2[l2]*(v1-vconf)*f(sin1,sin2,cos1,cos2);

};
V_ss::V_ss(double en) : V(en){l1=0;l2=0;};
V_sp::V_sp(double en) : V(en){l1=0;l2=1;};
V_pp_sig::V_pp_sig(double en) : V(en){l1=1;l2=1;};
V_pp_pi::V_pp_pi(double en) : V(en){l1=1;l2=1;};
V_pd_pi::V_pd_pi(double en) : V(en){l1=1;l2=2;};
V_pd_sig::V_pd_sig(double en) : V(en){l1=1;l2=2;};
V_sd::V_sd(double en) : V(en){l1=0;l2=2;};
V_dd_del::V_dd_del(double en) : V(en){l1=2;l2=2;};
V_dd_pi::V_dd_pi(double en) : V(en){l1=2;l2=2;};
V_dd_sig::V_dd_sig(double en) : V(en){l1=2;l2=2;};
