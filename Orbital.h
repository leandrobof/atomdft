/*
 * Orbital.h
 *
 *  Created on: Nov 14, 2016
 *      Author: leandro
 */
#include "globals.h"



#ifndef ORBITAL_H_
#define ORBITAL_H_

class Orbital {
public:
	Orbital(int principal,int angular,double num_ocupacion,int z,double energy,double *t,int N);
	virtual ~Orbital();
	virtual void outward(gsl_spline *,gsl_interp_accel *);
	virtual void inward(gsl_spline *,gsl_interp_accel *);
    virtual double correct_e();
    virtual void final(double r,state_type &R);
    virtual void inicial(double r,state_type &R);
    virtual void dens();
    void estimate_a();
    double get_energy(){return e;};
    void resolver(gsl_spline *,gsl_interp_accel *);
    virtual void print();
    virtual void radial();
    double operator[](int);
    double operator()(int);
    void complete_outward(double *t,gsl_spline *,gsl_interp_accel *);

protected:
	int n;
	int l;
	double e;
	int Z;
	int Nt;
	int nodos;
	double h;
	double tinf;
	double t0;
	double tk;
	int max;
	double noc;
	double co;
	double ci;
	double a;
    double *t;
	double dRi;
    double dRo;
    double *Rl;
    double *rho;

};
class Orbital_rel : public Orbital{
private:
	int k;
	double *Ql;
public:
	Orbital_rel(int principal,int angular,double num_ocupacion,int s,int z,double energy,double *t,int N);
	~Orbital_rel();
	virtual void dens();
    virtual void radial();
	virtual void final(double r,state_type &R);
    virtual void inicial(double r,state_type &R);
    virtual double correct_e();
    virtual void outward(gsl_spline *,gsl_interp_accel *);
    virtual void inward(gsl_spline *,gsl_interp_accel *);
    virtual void print();
};

#endif /* ORBITAL_H_ */
