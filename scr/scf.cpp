/*
 * scf.cpp
 *
 *  Created on: 24/05/2017
 *      Author: leandro
 */
#include "scf.h"

using namespace boost;

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
Scf::Scf(){
Z=0;
h=0.008;
max=50;
min=-8;
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
    vector<string> param;
	vector<string> config;
    char* str;
    prefix=archivo;
	string line;


    	map<string,int> options;
    	options["Relativistic"]=1;
    	options["GGA"]=2;
    	options["restart"]=3;
    	options["restart"]=4;
    	options["z"]=5;
    	options["sk_orb"]=6;





    	char_separator<char> sep("=, '");
    	while(getline(atomo,line)){
    		if(line.find("config")!=-1){
    			string subline=line.substr(line.find("="));
    			tokenizer<char_separator<char> > tokens(subline, sep);
    		    for(tokenizer<char_separator<char> >::iterator beg=tokens.begin(); beg!=tokens.end();++beg){
    		    	config.push_back(*beg) ;
    		    }
    		}
    		else{
    			tokenizer<char_separator<char> > tokens(line, sep);
    			for(tokenizer<char_separator<char> >::iterator beg=tokens.begin(); beg!=tokens.end();++beg){
    		    		  param.push_back(*beg) ;
    			}
    		}
    	}

    	atomo.close();

    for(int i=0;i<param.size();i++){
    	switch(options[param[i]]){

    	case 1:if(param[i+1]=="True"){
    			Relativistic=true;
    		}
    		break;

    	case 2:	if(param[i+1]=="True"){
    			gga=true;
    		}
    		break;

    	case 3:
    		if(param[i+1]=="True"){
    			restart=true;
    		}
    		break;

    	case 4:if(param[i+1]=="True"){
    		save=true;
    	};break;

    	case 5:Z=atoi(param[i+1].c_str());break;

    	case 6:
    		for (int j=0;j<3;j++){
    			SK[j]=param[i+1+j];

    		};
    		break;


    	default:break;
    	}
    }




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



//****************************************************************
int n;
int l;
double n_ocup;

if(Relativistic==false){


	string orbital;

	for(int i=0;i<config.size();i++){

				n_ocup=atof(config[i].substr(2).c_str());
				orbital=config[i].substr(0,2);
				n=int(orbital[0]-'0');
				l=x[orbital[1]];
				Atom.push_back(new Orbital_norel(n,l,n_ocup,Z,-Z*Z/(2.*n*n),t,N));
				index[orbital]=i;








	}

}
else{

	string orbital;
	int j=0;
	for(int i=0;i<config.size();i++){

		        n_ocup=atof(config[i].substr(2).c_str());
				orbital=config[i].substr(0,2);
				n=int(orbital[0]-'0');
				l=x[orbital[1]];
				index[orbital]=j;
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
				    j++;
					}

				else{Atom.push_back(new Orbital_rel(n,l,n_ocup,1,Z,-Z*Z/(2.*n*n),t,N));}
	            j++;


		}

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
    cout<<"Wood-Saxon potencial:  W: "<<W<<"  a:"<<a<<"  r0:"<<r0<<endl;

    cout<<"Functional :";
	if (gga==true){
		cout<<"PBE"<<endl;
	}else{
		cout<<"LDA"<<endl;
	}

	cout<<"Relativistic: ";
	if (Relativistic==true){
		cout<<"TRUE"<<endl;
	}
	else{
		cout<<"FALSE"<<endl;
	}


	cout<<"restart: ";
		if (restart==true){
			cout<<"TRUE"<<endl;
		}else{
			cout<<"FALSE"<<endl;
		}

	cout<<"Save: ";
	if (restart==true){
		cout<<"TRUE"<<endl;
	}
	else{
		cout<<"FALSE"<<endl;
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

void Scf::Vconf(Potential_spline &vc){
	vc.init(vconf,r,N);
}
void Scf::Veff_noconf(Potential_spline &vefnc){
	double v[N];
	for(int i=0;i<N;i++){
		v[i]=veff[i]-vconf[i];
	}
	vefnc.init(v,r,N);
};

void Scf::readpot(){
	double en;
	ifstream potfile(prefix+".pot");

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
	ofstream pot_file(prefix+".pot");
	for (int i=0;i<Atom.size();i++){
			pot_file<<Atom[i]->energy()<<endl;
		}

	for(int i=0;i<N;i++){
			pot_file<<veff[i]<<" "<<rho[i]<<endl;
	}

	pot_file.close();

}
