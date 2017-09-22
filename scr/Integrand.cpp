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
S::S( Orbital_spline *a, Orbital_spline *b) : Integrand(){
	R1=a->spline();
	R2=b->spline();
	accR1=a->acc();
	accR2=b->acc();
	r1max=a->max();
	r2max=b->max();

}

double S::calculate(double &r1,double &r2,double &sin1,double &sin2,double &cos1,double &cos2){

	double I;
	if(r1>r1max or r2>r2max){
		I=0;
	}
	else{
		I=(gsl_spline_eval(R1,r1,accR1)*gsl_spline_eval(R2,r2,accR2))*f(sin1,sin2,cos1,cos2);
	}

    return I;
};
S_ss::S_ss( Orbital_spline *a, Orbital_spline *b) : S(a,b){};
S_sp::S_sp( Orbital_spline *a, Orbital_spline *b) : S(a,b){};
S_pp_sig::S_pp_sig( Orbital_spline *a, Orbital_spline *b) : S(a,b){};
S_pp_pi::S_pp_pi( Orbital_spline *a, Orbital_spline *b) : S(a,b){};
S_sd::S_sd( Orbital_spline *a, Orbital_spline *b) : S(a,b){};
S_pd_pi::S_pd_pi( Orbital_spline *a, Orbital_spline *b) : S(a,b){};
S_pd_sig::S_pd_sig( Orbital_spline *a, Orbital_spline *b) : S(a,b){};
S_dd_del::S_dd_del( Orbital_spline *a, Orbital_spline *b) : S(a,b){};
S_dd_pi::S_dd_pi( Orbital_spline *a, Orbital_spline *b) : S(a,b){};
S_dd_sig::S_dd_sig( Orbital_spline *a, Orbital_spline *b) : S(a,b){};

V::V(Orbital_spline *a, Orbital_spline *b,Potential_spline &va,Potential_spline &vb,Potential_spline &vc ): S(a,b){
	v1=va.spline();
	accv1=va.acc();
	v2=vb.spline();
	accv2=vb.acc();
	vconf=vc.spline();
	accvconf=vc.acc();
	rmin=va.min();
    e=b->energy();
};
//Mejorar la evaluacion de integrales.Las casos menores a rmin no se evaluan actualmente.
double V::calculate(double &r1,double &r2,double &sin1,double &sin2,double &cos1,double &cos2){

	double I;
	if(r1>r1max or r2>r2max or r1<rmin or r2<rmin){
		I=0;

	}
	else{
		I=(gsl_spline_eval(R1,r1,accR1)*gsl_spline_eval(R2,r2,accR2))*(gsl_spline_eval(v1,r1,accv1)-gsl_spline_eval(vconf,r2,accvconf))*f(sin1,sin2,cos1,cos2);
	}

    return I;
};
V_ss::V_ss(Orbital_spline *a, Orbital_spline *b,Potential_spline &va,Potential_spline &vb,Potential_spline &vc) : V(a,b,va,vb,vc){};
V_sp::V_sp(Orbital_spline *a, Orbital_spline *b,Potential_spline &va,Potential_spline &vb,Potential_spline &vc) : V(a,b,va,vb,vc){};
V_pp_sig::V_pp_sig(Orbital_spline *a, Orbital_spline *b,Potential_spline &va,Potential_spline &vb,Potential_spline &vc) : V(a,b,va,vb,vc){};
V_pp_pi::V_pp_pi(Orbital_spline *a, Orbital_spline *b,Potential_spline &va,Potential_spline &vb,Potential_spline &vc) : V(a,b,va,vb,vc){};
V_pd_pi::V_pd_pi(Orbital_spline *a, Orbital_spline *b,Potential_spline &va,Potential_spline &vb,Potential_spline &vc) : V(a,b,va,vb,vc){};
V_pd_sig::V_pd_sig(Orbital_spline *a, Orbital_spline *b,Potential_spline &va,Potential_spline &vb,Potential_spline &vc) : V(a,b,va,vb,vc){};
V_sd::V_sd(Orbital_spline *a, Orbital_spline *b,Potential_spline &va,Potential_spline &vb,Potential_spline &vc) : V(a,b,va,vb,vc){};
V_dd_del::V_dd_del(Orbital_spline *a, Orbital_spline *b,Potential_spline &va,Potential_spline &vb,Potential_spline &vc) : V(a,b,va,vb,vc){};
V_dd_pi::V_dd_pi(Orbital_spline *a, Orbital_spline *b,Potential_spline &va,Potential_spline &vb,Potential_spline &vc) : V(a,b,va,vb,vc){};
V_dd_sig::V_dd_sig(Orbital_spline *a, Orbital_spline *b,Potential_spline &va,Potential_spline &vb,Potential_spline &vc) : V(a,b,va,vb,vc){};
