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
#include <boost/tokenizer.hpp>

enum ocupation{Average=0,Energy=1};



void dens_inicial(double *r,double *rho,int N,int Z);


/**
 * @brief Actualiza la densidad y gradiente a partir de la funciones radiales de los nuevos orbitales y la antigua densidad.
 * $$\rho_{new}=(1-\alpha)\rho+\alpha\rho_{old}$$
 * @param rho = array  densidad.
 * @param grad = array  gradiente.
 * @param Atom = vector con Orbitales, de donde se construira la nueva densidad.
 * @param alfa = parametro de damping
 * @param N = Numero de elementos de rho y grad;
 */

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
    ocupation ocup_type;
    string SK[3];
    string prefix;
    Xc *vxc;
    std::map<char,int> x;
    std::map<string,int> index;
    vector <Orbital*> Atom;
public:
	/**
	 * Constructor
	 * Construye scf a valores por defecto.
	 *
	 * @param tmin valor minimo de grilla radial t
	 * @param tmax valor maximo de grilla radial t
	 * @param step longitud del paso .
	 */
    Scf();
    ~Scf();
    /**
     * Lee archivo y obtiene Z, la configuracion electronica, tipo de calculo:Relativista no relativista,
     * tipo de Funcional.Con esto crea los orbitales.
     * @param Atom  vector de orbitales vacio.Donde se crearan los orbitales.
     * @param archivo Archivo con informacion sobre el atomo
     * @return
     */

    int initialize(string archivo);

    /**
     * Realiza calculo de orbitales de forma autoconsistente para un determinado potencial de confinamiento de Woods-Saxon.
      \f$ V_{conf}(r)=W/(1+e^{-a r/r_0}) \f$
     *
     *
     * @param Atom vector con Orbitales atomicos a resolver.
     * @param W    parametro de confinamiento.
     * @param a    parametro de confinamiento.
     * @param r0   parametro de confinamiento.
     * @param alf  damping.
     */
    void run(double W,double a,double r0,double alf);
    int z(){return Z;};
    int Nt(){return N;}
    /**
     * Obtiene array con autovalores y numero de ocupacion.En el caso Relativista devuelve un valor promedio
     * de acuerdo con ....
     *
     * @param Atom
     * @param e    array donde se guardan  los autovalores.
     * @param nocup
     */
    void energy(double *e,double *nocup);

    /**
         * Obtiene array de Orbitales_splines.En el caso Relativista devuelve un orbital promediado
         * de acuerdo con ....
         *
         * @param Atom
         * @param orb array donde se guardaran los splines de los orbitales
         */

    void orbital(Orbital_spline **orb);
    /**
     *
     * @return puntero a grilla t.
     */
    double* tgrid(){return t;}

    /**Devuelve puntero a grilla en r.
     *
     */
    double* rgrid(){return r;}

    /**
     *
     * @return puntero a Vconfinamiento.
     */


    void Vconf(Potential_spline &);
    void Veff_noconf(Potential_spline &);
    /**
     * Lee archivo con autovalores y veff.
     * @param Atom
     */
    void readpot();

    /**
     * Guarda archivo con autovalores  e y veff.
     * @param Atom
     */
    void savepot();
};


#endif /* SCF_H_ */
