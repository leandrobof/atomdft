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
#include <map>

void dens_inicial(double *r,double *rho,int N,int Z);
void new_rho(double *rho, vector<Orbital> &Atom,double alfa,int N);
class Scf{
private:
	int Z;
	int N;
	double h;
	double max;
	double min;
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
    bool gga;
    bool restart;
    bool save;
    string SK[3];
    Xc *vxc;
    std::map<char,int> x;
    std::map<string,int> index;
public:
	Scf(double tmin,double tmax,double step);
    ~Scf();
    int initialize(vector<Orbital*> &Atom,string archivo);
    void run(vector<Orbital*> &Atom,double W,double a,double r0,double alf);
    int z(){return Z;};
    int Nt(){return N;}
    void energy(vector<Orbital*> &Atom,double *e,double *nocup);
    void orbital(vector<Orbital*> &Atom,Orbital_spline **);
    double* tgrid(){return t;}
    double* rgrid(){return r;}
    double* Vconf(){return vconf;};
    double* Veff_noconf();
    void readpot(vector<Orbital*> &Atom);
    void savepot(vector<Orbital*> &Atom);
};


#endif /* SCF_H_ */
