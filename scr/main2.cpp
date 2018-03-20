/*
 * main.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: leandro
 */

#include "Orbital.h"
#include "gauss.h"
#include "SKtable.h"
#include "scf.h"
#include <sstream>

using namespace std;


//double Hubbard(Scf &a,vector<Orbital*> &b,const int &c );


int main(int argc,char *argv[]){


vector<string> simbolo;
vector<double*> e;
vector<double*> ocupation;
vector<double*> U;

vector<Potential_spline*> Veff0;
vector<Potential_spline*> Veff;
vector<Potential_spline*> vconf;

vector<Orbital_spline **> A;

vector<Scf*> Atoms;


ifstream file(argv[1]);
string line;
int i=0;

while(getline(file,line)){
	double param_dens[3]={0,0,0};
	double param_basis[3]={0,0,0};

//***Lee  Atomo el archivo, y los parametros de confinamiento

	string s=line.substr(line.find_first_not_of(" "));
	string archivo=s.substr(0,s.find(" "));
	simbolo.push_back(archivo);

	if(line.find("basis")!=-1){
		string basis=line.substr(line.find("basis"));
		basis=basis.substr(basis.find("{")+1,basis.find("}")-basis.find("{")-1);
		istringstream b(basis);
		b>>param_basis[0]>>param_basis[1]>>param_basis[2];


	}


	for (int j=0;j<3;j++){
		param_dens[j]=param_basis[j];
	}

	if(line.find("dens")!=-1){
		string dens=line.substr(line.find("dens"));
		dens=dens.substr(dens.find("{")+1,dens.find("}")-dens.find("{")-1);
		istringstream d(dens);
		d>>param_dens[0]>>param_dens[1]>>param_dens[2];
	}

//********************************************

	A.push_back(new Orbital_spline* [3] {});
	e.push_back(new double[3]);
	ocupation.push_back(new double[3]);
	U.push_back(new double[3]);
	Veff0.push_back(new Potential_spline);
	Veff.push_back(new Potential_spline);
	vconf.push_back(new Potential_spline);
	Atoms.push_back(new Scf);
// Inicializa y corre atomo sin confinamiento.
	Atoms[i]->initialize(simbolo[i]);
    Atoms[i]->run(0,1,1,0.2);
    Atoms[i]->energy(e[i],ocupation[i]);
// Corre  confinamiento densidad y obtienen pot.
    Atoms[i]->run(param_dens[0],param_dens[1],param_dens[2],0.2);
    Atoms[i]->Veff_noconf(*(Veff0[i]));
// Corre confinamientos bases y obtiene bases y potenciales.
    Atoms[i]->run(param_basis[0],param_basis[1],param_basis[2],0.2);
	Atoms[i]->orbital(A[i]);
	Atoms[i]->Vconf(*(vconf[i]));
	Atoms[i]->Veff_noconf(*(Veff[i]));

	i++;
}


for (int i=0;i<Atoms.size();i++){
	for (int j=0;j<Atoms.size();j++){
		SKtable sk;
		if (i !=j){
			sk.create(A[i],A[j],false);
		}
		else {
			sk.create(A[i],A[j],true);
		};

		sk.run(A[i],A[j],*(Veff0[i]),*(Veff0[j]),*(Veff[j]),*(vconf[j]),e[i],U[i],ocupation[i],simbolo[i]+"-"+simbolo[j]+".skf");

		cout<<simbolo[i]+"-"+simbolo[j]+".skf"<<endl;
	}
}


for(int k=0;k<Atoms.size();k++){
	for(int i=0;i<3;i++){
		delete A[k][i];
	}
}

for (int i=0;i<Atoms.size();i++){
	delete Atoms[i];
    delete Veff[i];
    delete vconf[i];
    delete Veff0[i];
    delete [] e[i];
    delete [] ocupation[i];
    delete [] U[i];

}
file.close();

//*******************
return 0;
}


/*
double Hubbard(Scf &a,vector<Orbital*> &b,const int &c ){
	double h=0.001;
	double n=b[c]->ocup();
	double nf=n+0.5*h;
	double nb=n-0.5*h;
	double e=b[c]->energy();
	double ef;
	double eb;
	b[c]->set_ocup(nf);
	a.run(b,0,1,1,0.3);
	ef=b[c]->energy();
	b[c]->set_ocup(nb);
	a.run(b,0,1,1,0.3);
	eb=b[c]->energy();
	b[c]->set_ocup(n);
	a.run(b,0,1,1,0.3);
	return (ef-eb)/h;
}
*/

