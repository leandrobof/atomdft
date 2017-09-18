//============================================================================
// Name        : xc_potential.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "xc_potential.h"


Xc::Xc(int N){

	ex=new double [N];
    ec=new double [N];
    vx=new double [N];
    vc=new double [N];
    vxc=new double [N];


};
Xc_lda::Xc_lda(int N,bool relativistic): Xc(N){
	xc_func_init(&correlacion, XC_LDA_C_VWN, XC_UNPOLARIZED);
	xc_func_init(&exchange, XC_LDA_X, XC_UNPOLARIZED);


    if (relativistic==true){
    	xc_lda_x_set_params(&exchange, 4.0/3.0,XC_RELATIVISTIC,0.);
    }


};

Xc_gga::Xc_gga(int N) : Xc(N){
	xc_func_init(&correlacion, XC_GGA_C_PBE, XC_UNPOLARIZED);
	xc_func_init(&exchange, XC_GGA_X_PBE, XC_UNPOLARIZED);

    vxsigma=new double [N];
    vcsigma=new double [N];
    sigma=new double [N];




};
Xc::~Xc(){
	delete [] ex;
	delete [] ec;
	delete [] vx;
	delete [] vc;
	delete [] vxc;

	xc_func_end(&correlacion);
	xc_func_end(&exchange);
};
Xc_lda::~Xc_lda(){};
Xc_gga::~Xc_gga(){

	delete [] vxsigma;
	delete [] vcsigma;
	delete [] sigma;

};
void Xc_lda::update(double *r,double *rho,double *grad,int N){
	//xc_lda_x_set_params(&exchange, 4.0/3.0,XC_RELATIVISTIC,0.);

	//********correlacion*********
	xc_lda_exc_vxc(&correlacion, N, rho, ec,vc);


    //********intercambio*************


    xc_lda_exc_vxc(&exchange, N, rho, ex,vx);

    //*****************

    //********** integracion************

    double I[N];

    for (int i=0;i<N;i++){
    	//integrando
    	I[i]=rho[i]*r[i]*r[i]*(ex[i]+ec[i]);
        //Potencial de intercambio y correlacion.Vxc=exc+n*(dexc/dn).


    	vxc[i]=vx[i]+vc[i];

    }


    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
    gsl_spline_init (spline, r, I,N );
    Exc=4*pi*gsl_spline_eval_integ (spline, r[0],r[N-1], acc);


    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    }
void Xc_gga::update(double *r,double *rho,double *grad,int N){
	//xc_lda_x_set_params(&exchange, 4.0/3.0,XC_RELATIVISTIC,0.);
    for(int i=0;i<N;i++){
		sigma[i]=grad[i]*grad[i];

    }
	//********correlacion*********

	xc_gga_exc_vxc(&correlacion, N, rho, sigma, ec, vc, vcsigma);

    //********intercambio*************



    xc_gga_exc_vxc(&exchange, N, rho, sigma, ex, vx, vxsigma);
    //*****************

    //********** integracion************

    double I[N];
    double x[N];

    //Si f=rho*exc
    //vxc=df/dn- div*((df/dsigma)*(dsigma/dgrad_rho)
    //En coordenadas esfericas div.F=(1/r^2)*d(r^2*F_r)/dr+...

    for (int i=0;i<N;i++){
    	//integrando
    	I[i]=rho[i]*r[i]*r[i]*(ex[i]+ec[i]);
        //Potencial de intercambio y correlacion.

    	x[i]=2*grad[i]*(vxsigma[i]+vcsigma[i])*r[i]*r[i];

    	vxc[i]=vx[i]+vc[i];

    }


    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
    gsl_spline_init (spline, r, I,N );
    Exc=4*pi*gsl_spline_eval_integ (spline, r[0],r[N-1], acc);
   	gsl_spline_init (spline, r, x,N );

    for (int i=0;i<N;i++){
    		vxc[i]=vxc[i]-gsl_spline_eval_deriv (spline, r[i], acc)/(r[i]*r[i]);
    	}

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    }

double Xc::operator[](int i){return vxc[i];}
