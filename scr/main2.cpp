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
void SK_table(Orbital_spline **A,Orbital_spline **B,Potential_spline &Veff,Potential_spline &Veff0,Potential_spline &Vconf,double *e,double *U,double *ocupation,string archivo);

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
	string str3=str2.substr(str2.find_first_not_of(" "),str2.find_last_not_of(" "));
	if(str3=="True"){
		Relativistic=true;
	}
	else if(str3=="False"){
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
	string str3=str2.substr(str2.find_first_not_of(" "),str2.find_last_not_of(" "));
	if(str3=="True" ){
		gga=true;
	}
	else if(str3=="False"){
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
map<string,int> index;

//Lee entrada y crea vector de Orbitales.
if(Relativistic==false){

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
	int i=0;
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
			index[orbital]=i;
			if(l>0){
				Atom.push_back(new Orbital_rel(n,l,n_ocup_alfa,1,Z,-Z*Z/(2.*n*n),t,scf.Nt()));
				Atom.push_back(new Orbital_rel(n,l,n_ocup_beta,-1,Z,-Z*Z/(2.*n*n),t,scf.Nt()));
			    i++;
			}
			else{Atom.push_back(new Orbital_rel(n,l,n_ocup_alfa+n_ocup_beta,1,Z,-Z*Z/(2.*n*n),t,scf.Nt()));}
            i++;
		}
	}
	atomo.close();
}

scf.run(Atom,0,1,1,0.4);


double e[3]={0.,0.,0.};
double ocupation[3]={0.,0.,0.};
double U[3]={0.,0.,0.};

if(Relativistic==false){
	for(int i=0;i<3;i++){
		if (index.count(SK[i])>0){                          //Comprueba si el orbitalseleccionado para las tabla es un orbital calculado.
			e[i]=Atom[index[SK[i]]]->energy();
			ocupation[i]=Atom[index[SK[i]]]->ocup();
			//U[i]=Hubbard(scf,Atom,index[SK[i]]);          // Desmarcar para habilitar Hubbard

		}
	}
}
else{
	for(int i=0;i<3;i++){
		if (index.count(SK[i])>0){                           //Comprueba si el orbitalseleccionado para las tabla es un orbital calculado.
			l=x[SK[i][1]];
			if(l>0){
			    e[i]=(l/(2.*l+1.))*Atom[index[SK[i]]+1]->energy()+((l+1.)/(2.*l+1.))*Atom[index[SK[i]]]->energy();
			    ocupation[i]=Atom[index[SK[i]]]->ocup()+Atom[index[SK[i]]+1]->ocup();
			    //U[i]=Hubbard(scf,Atom,index[SK[i]]);          // Desmarcar para habilitar Hubbard


			}
		    else{
		    	e[i]=Atom[index[SK[i]]]->energy();
		    	ocupation[i]=Atom[index[SK[i]]]->ocup();


		    }

		}
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

if (Relativistic==true){
	for(int i=0;i<3;i++){
		if(index.count(SK[i])>0){
			double o[N];
			l=x[SK[i][1]];
			if(l>0){
			for(int j=0;j<N;j++){
				o[j]=(*(Atom[index[SK[i]]+1]))(j)*(l/2.*l+1.)+(*(Atom[index[SK[i]]]))(j)*(l+1.)/(2.*l+1.);
			}
			}
			else{
				for(int j=0;j<N;j++){
					o[j]=(*(Atom[index[SK[i]]]))(j);
				}
			}

			C[i]=new Orbital_spline(o,r,e[i],l,N);
		}
	}
}
else{
for(int i=0;i<3;i++){
	if(index.count(SK[i])>0){
		C[i]=new Orbital_spline(*(Atom[index[SK[i]]]),r,N);
	}

}}

for(int i=0;i<N;i++){
	veff2[i]=veff[i]-vconf[i];
}

Potential_spline Veff(veff2,r,N);
Potential_spline Vconf(vconf,r,N);

string simbolo(argv[1]);

SK_table(C,C,Veff,Veff0,Vconf,e,U,ocupation,simbolo+"-"+simbolo+".skf");


delete C[0];
delete C[1];
delete C[2];

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

void SK_table(Orbital_spline **A,Orbital_spline **B,Potential_spline &Veff,Potential_spline &Veff0,Potential_spline &Vconf,double *e,double *U,double *ocupation,string archivo){
vector<Integrand*> S(10,NULL);
vector<Integrand*> V(10,NULL);

double d=0.2;
create(A,B,S,V,Veff,Veff0,Vconf);


gauss g(50,d);

ofstream salida2(archivo);

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


borrar(S,V);
salida2.close();
}

