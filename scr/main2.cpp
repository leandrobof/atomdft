/*
 * main.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: leandro
 */

#include "Orbital.h"
#include "gauss.h"
#include <map>
#include "scf.h"

using namespace std;

void create( Orbital_spline **a,Orbital_spline **b,vector<Integrand*> &c,vector<Integrand*> &d,Potential_spline &va,Potential_spline &vb,Potential_spline &vc);
void borrar(vector<Integrand*> &,vector<Integrand*> &);
double Hubbard(Scf &a,vector<Orbital*> &b,const int &c );
int main(int argc,char *argv[]){
int Z;
int n;
int l;
double n_ocup;
bool Relativistic=false;
bool gga=false;
vector <Orbital*> Atom;
std::map<char,int> x;
x['s']=0;
x['p']=1;
x['d']=2;
x['f']=3;

ifstream atomo(argv[1]);
string str;
string line2;
getline(atomo,str);
getline(atomo,line2);
atomo>>Z;

if(str.find("Relativistic")==-1){
	cout<<"Keyword \"Relativistic\" not found\n"<<endl;
return 1;
}
else{
	string str2=str.substr(str.find("Relativistic")+12);
	if(str2.find("True") != -1){
		Relativistic=true;
	}
	else if(str2.find("False")!=-1){
		Relativistic=false;
	}else {cout<<"\"Relativistic\" has only values \"True\" or \"False\"\n"<<endl;
	return 1;}
}

if(line2.find("GGA")==-1){
	cout<<"Keyword \"GGA\" not found\n"<<endl;
return 1;
}
else{
	string str2=line2.substr(line2.find("GGA")+3);
	if(str2.find("True") != -1){
		gga=true;
	}
	else if(str2.find("False")!=-1){
		gga=false;
	}else {cout<<"\"GGA\" has only values \"True\" or \"False\"\n"<<endl;
	return 1;}
}

double tinf=log(Z*50);
double t0=-8;
double h=0.001;




Scf scf(Z,t0,tinf,h,Relativistic,gga);

double *t=scf.tgrid();

scf.initialize();

//Lee que orbitales usar para tablas SK
// n1_s n2_p n3_d
string SK[3];
atomo>>SK[0];atomo>>SK[1];atomo>>SK[2];

//Lee entrada y crea vector de Orbitales.
if(Relativistic==false){
map<string,int> index;
cout<<Z<<endl;
string orbital;
int i=0;
if(atomo.is_open()){
	while(atomo.good()){
		atomo>>orbital;
		atomo>>n_ocup;
		if(atomo.eof())break;

		n=int(orbital[0]-'0');
		l=x[orbital[1]];
		Atom.push_back(new Orbital_norel(n,l,n_ocup,Z,-Z*Z/(2.*n*n),t,scf.Nt()));
        index[orbital]=i;
        i++;

	}
}
atomo.close();
}
else{
	double n_ocup_alfa,n_ocup_beta;
	string orbital;
	if(atomo.is_open()){
		while(atomo.good()){
			atomo>>orbital;
			atomo>>n_ocup_alfa;
			atomo>>n_ocup_beta;
			if(atomo.eof())break;

			n=int(orbital[0]-'0');
			l=x[orbital[1]];
			if(l>0){
				Atom.push_back(new Orbital_rel(n,l,n_ocup_alfa,1,Z,-Z*Z/(2.*n*n),t,scf.Nt()));
				Atom.push_back(new Orbital_rel(n,l,n_ocup_beta,-1,Z,-Z*Z/(2.*n*n),t,scf.Nt()));
			}
			else{Atom.push_back(new Orbital_rel(n,l,n_ocup_alfa+n_ocup_beta,1,Z,-Z*Z/(2.*n*n),t,scf.Nt()));}

		}
	}
	atomo.close();
}

scf.run(Atom,0,1,1,0.2);

/*
double e[3]={0.,0.,0.};
double ocupation[3]={0.,0.,0.};
double U[3]={0.,0.,0.};

for(int i=0;i<3;i++){
	if (index.count(SK[i])>0){                           //Comprueba si el orbitalseleccionado para las tabla es un orbital calculado.
		e[i]=Atom[index[SK[i]]]->energy();
		ocupation[i]=Atom[index[SK[i]]]->ocup();
        U[i]=Hubbard(scf,Atom,index[SK[i]]);          // Desmarcar para habilitar Hubbard
	}

}

scf.run(Atom,atof(argv[2]),atof(argv[3]),atof(argv[4]),0.2);

int N=scf.Nt();
double *r=scf.rgrid();
double *veff=scf.Veff_noconf();
double *vconf=scf.Vconf();
double veff0[N];
double veff2[N];
//potencial de referencia V0
for(int i=0;i<N;i++){
	veff0[i]=veff[i]-vconf[i];
}

Potential_spline Veff0(veff0,r,N);
//*****

scf.run(Atom,atof(argv[5]),atof(argv[6]),atof(argv[7]),0.3);

//*****
//Array  splines de los orbitales.
Orbital_spline *C[3]={NULL,NULL,NULL};

for(int i=0;i<3;i++){
	if(index.count(SK[i])>0){
		C[i]=new Orbital_spline(*(Atom[index[SK[i]]]),r,N);
	}

}

for(int i=0;i<N;i++){
	veff2[i]=veff[i]-vconf[i];
}

Potential_spline Veff(veff2,r,N);
Potential_spline Vconf(vconf,r,N);


vector<Integrand*> S(10,NULL);
vector<Integrand*> V(10,NULL);

double d=0.2;
create(C,C,S,V,Veff,Veff0,Vconf);


gauss g(50,d);

ofstream salida2("Ni-Ni.skf");

salida2<<"0.02"<<" "<<"599"<<endl;
salida2<<e[2]<<" "<<e[1]<<"  "<<e[0]<<"  0.0  "<<U[2]<<"  "<<U[1]<<"  "<<U[0]<<"  "<<ocupation[2]<<"  "<<ocupation[1]<<"  "<<ocupation[0]<<endl;
salida2<<"12.01, 19*0.0"<<endl;

for(int i=0;i<19;i++){
	salida2<<"20*1.0,"<<endl;
}

double overlap[10];
double H[10];
double s[10];
double v[10];

for(int i=0;i<10;i++){
		H[i]=0.00;
		overlap[i]=0.00;
		s[i]=0.0;
		v[i]=0.0;
}

while(2*d<12){
	g.integrate2d(0,1,0,pi,S,V,s,v);

    for(int i=0;i<10;i++){
    	if(V[i]!=NULL){
    		overlap[i]=s[i];
    		H[i]=V[i]->energy()*overlap[i]+v[i];
    	}
	}


	for(int i=9;i>=0;i--){
		salida2<<H[i]<<"  ";
	}

	for(int i=9;i>=0;i--){
			salida2<<overlap[i]<<"  ";
	}
    for(int i=0;i<10;i++){
    	s[i]=0.;
    	v[i]=0.;

    }


	salida2<<endl;
    d=d+0.01;
    g.update_a(d);


}

salida2<<"Spline"<<endl;
salida2<<"2 1.28"<<endl;
salida2<<"2.151029456234113 3.917667206325493 -0.4605879014976964"<<endl;
salida2<<"1.2 1.24    3.344853 -8.185615473079642 8.803750000000022 1.68154567477936"<<endl;
salida2<<"1.24 1.28    0.016 -0.006590813456982203 -0.02356970905317782 -0.09209220073124012 0.2061755069509315 -0.1001089592255145"<<endl;
borrar(S,V);
salida2.close();


delete C[0];
delete C[1];
delete C[2];
*/
for (int i=0;i<Atom.size();i++){
	delete Atom[i];
}

//*******************
return 0;
}


void create( Orbital_spline **a,Orbital_spline **b,vector<Integrand*> &c,vector<Integrand*> &d,Potential_spline &va,Potential_spline &vb,Potential_spline &vc){
    if(a[0]!=NULL and b[0]!=NULL){
    	c[0]=new S_ss(a[0],b[0]);
    	d[0]=new V_ss(a[0],b[0],va,vb,vc);
    }

    if(a[0]!=NULL and b[1]!=NULL){
    	c[1]=new S_sp(a[0],b[1]);
    	d[1]=new V_sp(a[0],b[1],va,vb,vc);
    }

    if(a[0]!=NULL and b[2]!=NULL){
    	c[2]=new S_sd(a[0],b[2]);
    	d[2]=new V_sd(a[0],b[2],va,vb,vc);
    }

    if(a[1]!=NULL and b[1]!=NULL){
    	c[3]=new S_pp_pi(a[1],b[1]);
    	d[3]=new V_pp_pi(a[1],b[1],va,vb,vc);
    	c[4]=new S_pp_sig(a[1],b[1]);
    	d[4]=new V_pp_sig(a[1],b[1],va,vb,vc);
    }

    if(a[1]!=NULL and b[2]!=NULL){
    	c[5]=new S_pd_pi(a[1],b[2]);
    	d[5]=new V_pd_pi(a[1],b[2],va,vb,vc);
    	c[6]=new S_pd_sig(a[1],b[2]);
    	d[6]=new V_pd_sig(a[1],b[2],va,vb,vc);
    }

    if(a[2]!=NULL and b[2]!=NULL){
    	c[7]=new S_dd_del(a[2],b[2]);
    	d[7]=new V_dd_del(a[2],b[2],va,vb,vc);
    	c[8]=new S_dd_pi(a[2],b[2]);
    	d[8]=new V_dd_pi(a[2],b[2],va,vb,vc);
    	c[9]=new S_dd_sig(a[2],b[2]);
    	d[9]=new V_dd_sig(a[2],b[2],va,vb,vc);
    }


};



void borrar(vector<Integrand*> &c,vector<Integrand*> &d){
	for(int i=0;i<c.size();i++){
		delete c[i];
		delete d[i];
	}
}
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
