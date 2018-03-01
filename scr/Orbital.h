/*
 * Orbital.h
 *
 *  Created on: Nov 14, 2016
 *      Author: leandro
 */
#include "globals.h"
#include "Ode.h"


#ifndef ORBITAL_H_
#define ORBITAL_H_

double simpson(double *f,double h,int N);

class Orbital {
public:
	Orbital(int principal,int angular,double num_ocupacion,int z,double energy,double *t,int N);
	Orbital(const Orbital &orb);
	virtual ~Orbital();
	virtual void outward(double *,double *)=0;
	virtual void inward(double*,double *,double W)=0;
    virtual double correct_e(double *r)=0;
    virtual void final(double r,double *R,double W)=0;
    virtual void inicial(double r,double *R)=0;
    virtual void dens(double *r)= 0;
    void estimate_a();
    double energy(){return e;};
    int prin(){return n;};
    int ang(){return l;};
    void set_ocup(double o){noc=o;}
    void set_e(double en){e=en;}
    double ocup(){return noc;};
    void resolver(double *,double *,double W);
    virtual void print()=0;
    virtual void radial()=0;
    double operator[](int);
    double operator()(int);
    double grd(int i){return grad[i];};


protected:
	int n;
	int l;
	double e;
	int Z;
	int Nt;
	int nodos;
	double h;
	double *t;
	double tinf;
	double t0;
	double tk;
	double noc;
	double co;
	double ci;
	double a;
    double dRi;
    double dRo;
    double *Rl;
    double *rho;
    double *grad;
    int max;
};
class Orbital_norel :public Orbital{

public:
	Orbital_norel(int principal,int angular,double num_ocupacion,int z,double energy,double *t,int N);
	virtual ~Orbital_norel();
	virtual void dens(double *r);
    virtual void radial();
	virtual void final(double r,double *R,double W);
    virtual void inicial(double r,double *R);
    virtual double correct_e(double *r);
    virtual void outward(double *,double *);
    virtual void inward(double*,double *,double W);
    virtual void print();
};
class Orbital_rel : public Orbital{
private:
	int k;
	double *Ql;
public:
	Orbital_rel(int principal,int angular,double num_ocupacion,int s,int z,double energy,double *t,int N);
	Orbital_rel(const Orbital_rel &orb);
	Orbital_rel operator=( const Orbital_rel& orb );
	~Orbital_rel();
	virtual void dens(double *r);
    virtual void radial();
	virtual void final(double r,double *R,double W);
    virtual void inicial(double r,double *R);
    virtual double correct_e(double *r);
    virtual void outward(double *,double *);
    virtual void inward(double*,double *,double W);
    virtual void print();
};


class Orbital_spline{
private:
	gsl_spline *R;
	gsl_interp_accel *accR;
    char l;
    double rmax;
    double e;
public:
    Orbital_spline();
    Orbital_spline( Orbital &,double *,int N);
    Orbital_spline(double *,double *,double energy,int angl ,int N);
    Orbital_spline(const Orbital_spline &);
    Orbital_spline & operator=(const Orbital_spline &a);
    gsl_spline * spline(){return R;}
    gsl_interp_accel * acc(){return accR;}
    double max(){return rmax;}
    double operator()(double);
	virtual ~Orbital_spline();
    inline char ang(){return l;};
    double energy(){return e;}
};
class Potential_spline{
private:
	gsl_spline *v;
	gsl_interp_accel *accv;
	double rmin;
public:
	Potential_spline();
	Potential_spline(double *,double *,int N);
	void init(double *,double *,int N);
	virtual ~Potential_spline();
	gsl_spline * spline(){return v;}
	double operator()(double);
	gsl_interp_accel * acc(){return accv;};
    double min(){return rmin;};
};
#endif /* ORBITAL_H_ */
