/*
 * scf.h
 *
 *  Created on: 24/05/2017
 *      Author: leandro
 */

#ifndef SCF_H_
#define SCF_H_
#include "Orbital.h"
#include "hartree.h"
#include "xc_potential.h"

void dens_inicial(double *r,double *rho,int N,int Z);
void new_rho(double *rho, vector<Orbital> &Atom,double alfa,int N);
class Scf{
private:
	int Z;
	int N;
	double h;
	double *t;
	double *r;
	double *veff;
	double *vn;
	double *vconf;
	double *rho;
	double *grad;
	double alfa;
	double W;
	double a;
	double r0;
	bool Relativistic;

    Xc *vxc;
public:
	Scf(int z,double tmin,double tmax,double step,bool,bool);
    ~Scf();
    void initialize();
    void run(vector<Orbital*> &Atom,double W,double a,double r0,double alf);
    int z(){return Z;};
    int Nt(){return N;}
    double* tgrid(){return t;}
    double* rgrid(){return r;}
    double* Vconf(){return vconf;};
    double* Veff_noconf();
};


#endif /* SCF_H_ */
