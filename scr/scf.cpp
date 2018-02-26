/*
 * scf.cpp
 *
 *  Created on: 24/05/2017
 *      Author: leandro
 */
#include "scf.h"

void dens_inicial(double *r,double *rho,int N,int Z){
	double x;
	double a=0.7280642371;
	double b=-0.5430794693;
	double gamma=0.3612163121;
	double Zeff;
	double v;
	for(int i=0;i<N;i++){
		x=r[i]*pow((128*Z/9*pi*pi),1/3.);
		Zeff=Z*(1+a*sqrt(x)+b*x*exp(-gamma*sqrt(x)))*(1+a*sqrt(x)+b*x*exp(-gamma*sqrt(x)))*exp(-2*a*sqrt(x));
        v=Zeff/r[i];
        rho[i]=-1/(3*pi*pi)*pow(-2*v,3/2.);
	}
}


void new_rho(double *rho, double *grad,vector<Orbital*> &Atom,double alfa,int N){
	double nrho;
	double ngrad;
	for(int i=0;i<N;i++){
		nrho=0;
		ngrad=0;
		for (int j=0;j<Atom.size();j++){
	         nrho=nrho+(Atom[j])->operator[](i);
             ngrad=ngrad+(Atom[j]->grd(i));
	    }
		rho[i]=(1-alfa)*rho[i]+alfa*nrho;
		grad[i]=(1-alfa)*grad[i]+alfa*ngrad;

	}
};
Scf::Scf(double tmin,double tmax,double step){
Z=0;
h=step;
max=tmax;
min=tmin;
alfa=0.3;
W=0;
a=1;
r0=1;
Relativistic=false;
gga=false;
restart=false;
save=false;
ocup_type=Average;
x['s']=0;
x['p']=1;
x['d']=2;
x['f']=3;

}

int Scf::initialize(string archivo){
	ifstream atomo(archivo);

	string line1;
	string line2;
    string line3;
	string line4;
	getline(atomo,line1);
	getline(atomo,line2);
	getline(atomo,line3);
	getline(atomo,line4);
	atomo>>Z;
	atomo>>SK[0];atomo>>SK[1];atomo>>SK[2];
	atomo.ignore(256, '\n');
	N=(log(Z*max)-min)/h;
	t=new double [N];
	r=new double [N];
	veff=new double [N];
	vn=new double [N];
	vconf=new double [N];
	rho=new double [N];
	grad=new double [N];

/*****Genera malla uniforme t, y potencial*****************/
	for (int i=0;i<N;i++){
		t[i]=min+i*h;
		r[i]=exp(t[i])/Z;
		vn[i]=-Z/r[i];
		rho[i]=0;
		vconf[i]=0.;
		veff[i]=vn[i]+vconf[i];
	}




	if(line1.find("Relativistic")==-1){
		cout<<"Keyword \"Relativistic\" not found\n"<<endl;
	return 1;
	}
	else{
		string str2=line1.substr(line1.find("Relativistic")+12);
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

    if(line3.find("restart")==-1){
			cout<<"Keyword \"restart\" not found\n"<<endl;
		return 1;
		}
		else{
			string str2=line3.substr(line3.find("restart")+7);
			string str3=str2.substr(str2.find_first_not_of(" "),str2.find_last_not_of(" "));
			if(str3=="True" ){
				restart=true;
			}
			else if(str3=="False"){
				restart=false;
			}else {cout<<"\"restart\" has only values \"True\" or \"False\"\n"<<endl;
			return 1;}
		}

    if(line4.find("save")==-1){
 			cout<<"Keyword \"save\" not found\n"<<endl;
 		return 1;
 		}
 		else{
 			string str2=line4.substr(line4.find("save")+4);
 			string str3=str2.substr(str2.find_first_not_of(" "),str2.find_last_not_of(" "));
 			if(str3=="True" ){
 				save=true;
 			}
 			else if(str3=="False"){
 				save=false;
 			}else {cout<<"\"save\" has only values \"True\" or \"False\"\n"<<endl;
 			return 1;}
 		}

//****************************************************************
int n;
int l;
double n_ocup;

if(Relativistic==false){


	string orbital;
	string configuration;
	int i=0;
	if(atomo.is_open()){
		getline(atomo,configuration);
			int pos_ini=0;
			int pos_fin=0;

			while((pos_fin=configuration.find(" ",pos_ini)) != string::npos){

				n_ocup=atof((configuration.substr(pos_ini+2,(pos_fin-pos_ini-2))).c_str());
				orbital=configuration.substr(pos_ini,2);
				n=int(orbital[0]-'0');
				l=x[orbital[1]];
				Atom.push_back(new Orbital_norel(n,l,n_ocup,Z,-Z*Z/(2.*n*n),t,N));
				index[orbital]=i;
				pos_ini=pos_fin+1;
                i++;

			}

			;


	}
	atomo.close();
}
else{

	string orbital;
	string configuration;
	int i=0;
	if(atomo.is_open()){
		getline(atomo,configuration);
			int pos_ini=0;
			int pos_fin=0;

			while((pos_fin=configuration.find(" ",pos_ini)) != string::npos){

				n_ocup=atof((configuration.substr(pos_ini+2,(pos_fin-pos_ini-2))).c_str());
				orbital=configuration.substr(pos_ini,2);
				n=int(orbital[0]-'0');
				l=x[orbital[1]];
				index[orbital]=i;
				if(l>0){
					double n_ocup_1=0;
					double n_ocup_2=0;
					switch(ocup_type){

					case Average:
						n_ocup_1=(2.*l)/(2.*l+1.)*n_ocup/2.;
					    n_ocup_2=(2.*l+2.)/(2.*l+1.)*n_ocup/2.;
					break;

					case Energy:
						for (int j=n_ocup;j>0;j--){
							if(n_ocup_1<2*l){
								n_ocup_1++;
							}
							else{
								n_ocup_2++;
							}
						}
					break;
					}

					Atom.push_back(new Orbital_rel(n,l,n_ocup_2,1,Z,-Z*Z/(2.*n*n),t,N));
					Atom.push_back(new Orbital_rel(n,l,n_ocup_1,-1,Z,-Z*Z/(2.*n*n),t,N));
				    i++;
					}

				else{Atom.push_back(new Orbital_rel(n,l,n_ocup,1,Z,-Z*Z/(2.*n*n),t,N));}
	            i++;
	            pos_ini=pos_fin+1;

			}
		}
		atomo.close();
	}



/***Genera malla uniforme t, y potencial*****************/


if(gga==false){
	vxc=new Xc_lda(N,Relativistic);
}
else{
	vxc=new Xc_gga(N);
}

if(restart==true){
	readpot();
}
};

Scf::~Scf(){
	delete [] t;
	delete [] r;
	delete [] veff;
	delete [] vn;
	delete [] vconf;
	delete [] rho;

    delete [] grad;
    delete vxc;
    for (int i=0;i<Atom.size();i++){
    	delete Atom[i];
    }

}



void Scf::run(double w,double al,double ro,double alf){
	/**************************************************************/
	alfa=alf;
	W=w;
	a=al;
	r0=ro;

	for(int i=0;i<N;i++){
		veff[i]=veff[i]-vconf[i];
		vconf[i]=W/(1+exp(-a*(r[i]-r0)));
		veff[i]=veff[i]+vconf[i];

	}

	/*Se crea spline y Orbital.Se resuelve el atomo hidrogenoide veff=-Z/r */



	for(int i=0;i<Atom.size();i++){
		Atom[i]->resolver(veff,r,W);

	}


	/*************************************************************/
	//Se obtiene rho inicial a partir del atomo hidrogenoide

	new_rho(rho,grad,Atom, alfa, N);
	//Se crea e inicializa los potenciales de hartree e intercambio y correlacion.y se crea veff
	Hartree vh(r,rho,t,h,Z,N);

    vxc->update(r,rho,grad,N);
	for(int i=0;i<N;i++){
	        veff[i]=vn[i]+vh[i]+vxc->operator[](i)+vconf[i];
	     }



	for(int i=0;i<Atom.size();i++){
		Atom[i]->resolver(veff,r,W);

	}

	double *e1;
	double *e2;
	e1=new double [Atom.size()];
	e2=new double [Atom.size()];



	for(int i=0;i<Atom.size();i++){
		e1[i]=0.;
		e2[i]=Atom[i]->energy();

	}

	double error=0.00000001;
	int iteraciones=0;
	while(fabs(e2[Atom.size()-1]-e1[Atom.size()-1])> error ){
	    new_rho(rho,grad,Atom, alfa, N);
		vh.update(r,rho,t,h,Z,N);
		vxc->update(r,rho,grad,N);
	    for(int i=0;i<N;i++){
	        veff[i]=vn[i]+vh[i]+vxc->operator[](i)+vconf[i];
	     }
	     e1[Atom.size()-1]=e2[Atom.size()-1];


	     for(int i=0;i<Atom.size();i++){
	    	 Atom[i]->resolver(veff,r,W);
	    	 e2[i]=Atom[i]->energy();

	     }
	     iteraciones++;
	     if(fabs(e2[Atom.size()-1]-e1[Atom.size()-1])<0.00001){
	    	 alfa=0.5;
	     }
	     if (iteraciones>100)break;
	}
	cout<<"iteraciones: "<<iteraciones<<endl;
	cout<<"Eh"<<"   "<<vh.energy()<<endl;
	cout<<"Exc"<<"  "<<vxc->energy()<<endl;

	for (int i=0;i<Atom.size();i++){
		cout<<Atom[i]->prin()<<angular[Atom[i]->ang()]<<"  "<<Atom[i]->energy()<<endl;

	}
	//***Se imprime  Rln de todos los orbitales en archivo radial.txt

	ofstream salida("radial2.txt");
	for(int i=0;i<N;i++){
		salida<<exp(t[i])/Z<<"   ";
		for(int j=0;j<Atom.size();j++){
			salida<<Atom[j]->operator()(i)<<"   ";
		}

		salida<<endl;
	}
	salida.close();

	/******Se libera spline y acc*********************/
    if(save==true){
    	savepot();
    }

	delete [] e1;
	delete [] e2;

}
void Scf::energy(double *e,double *ocupation){
    int l;
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
}
void Scf::orbital(Orbital_spline **C){
	int l;
	if (Relativistic==true){
		double e[3];
		double nocup[3];
		energy(e,nocup);
		for(int i=0;i<3;i++){
			if(index.count(SK[i])>0){
				double o[N];
				l=x[SK[i][1]];
				if(l>0){
					for(int j=0;j<N;j++){
						o[j]=(*(Atom[index[SK[i]]+1]))(j)*l/(2.*l+1.)+(*(Atom[index[SK[i]]]))(j)*(l+1.)/(2.*l+1.);
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
		}
	}
}

double* Scf::Veff_noconf(){
	/*for(int i=0;i<N;i++){
		veff[i]=veff[i]-vconf[i];
	}*/
return veff;};

void Scf::readpot(){
	double en;
	ifstream potfile("pot");

	for (int i=0;i<Atom.size();i++){
		    potfile>>en;
			Atom[i]->set_e(en);
		}

	for(int i=0;i<N;i++){
		potfile>>veff[i]>>rho[i];

	}

    potfile.close();
}

void Scf::savepot(){
	ofstream pot_file("pot");
	for (int i=0;i<Atom.size();i++){
			pot_file<<Atom[i]->energy()<<endl;
		}

	for(int i=0;i<N;i++){
			pot_file<<veff[i]<<" "<<rho[i]<<endl;
	}

	pot_file.close();

}
