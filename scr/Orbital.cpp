/*
 * Orbital.cpp
 *
 *  Created on: Nov 14, 2016
 *      Author: leandro
 */
#include "Orbital.h"
#include "Func.h"




Orbital::Orbital(int principal,int angular,double num_ocupacion,int z,double e_orbital,double *r,int N) {


n=principal;
l=angular;
e=e_orbital;
Z=z;
Nt=N;
t=r;
h=t[1]-t[0];
tinf=t[N];
t0=t[0];
tk=0;
a=1;
noc=num_ocupacion;
Rl=new double [Nt];
rho=new double [Nt];
grad=new double [Nt];
max=Nt;

}
Orbital::Orbital(const Orbital &orb) {


n=orb.n;
l=orb.l;
e=orb.e;
Z=orb.Z;
Nt=orb.Nt;
t=orb.t;
h=orb.h;
tinf=orb.tinf;
t0=orb.t0;
tk=0;
a=1;
noc=orb.noc;
max=Nt;
Rl=new double [Nt];
rho=new double [Nt];

}


Orbital::~Orbital() {
	// TODO Auto-generated destructor stub
	delete [] Rl;
    delete [] rho;

}
Orbital_norel::Orbital_norel(int principal,int angular,double num_ocupacion,int z,double e_orbital,double *r,int N): Orbital(principal, angular, num_ocupacion, z, e_orbital,r, N){};
Orbital_norel::~Orbital_norel(){

};
void Orbital_norel::inward(double *veff,double *r,double W){


	  Func Sch(e,l,Z,veff,t,r);
	  double  y[2];
	  y[0]=0;
	  y[1]=0;
	  int i=Nt-1;
	  while(y[0]<1.0e-15){
	  	  final(r[i],y,W);
	      i--;
	  }
	  max=i;
	  Ode<Func> abm;
	  for(int l=0;l<4;l++){
	      	  final(r[i],y,W);
	      	  abm.inicializar(Sch,y,i);
	      	  Rl[i]=y[0];
	      	  i--;


	  }

      final(r[i],y,W);
	  Rl[i]=y[0];
	  while (t[i] > tk)
	    {

	      abm.step( Sch, y ,-h ,i  );



	    i--;
	    Rl[i]=y[0];
	    }

	  ci=y[0];dRi=y[1];
	  for(int k=max;k<Nt;k++){
	      	  Rl[k]=0.;
	  }

};



void Orbital_norel::outward(double *veff,double *r){


	nodos=0;

	Func Sch(e,l,Z,veff,t,r);
	double y[2];
    double ya;
	Ode<Func> abm;
	for (int i=0;i<4;i++){
	  inicial(r[i],y);
	  abm.inicializar(Sch,y,i);
	  Rl[i]=y[0];

	 }
	double x0,x1;
    int i=4;
	x0=veff[i]-e;
	x1=x0;
	inicial(r[i],y);
	Rl[i]=y[0];
	while (x0*x1>0 ){
    	x0=x1;
    	ya=y[0];
        abm.step( Sch , y ,h, i );
        if (ya*y[0]<0){
        	nodos++;
    	}
    	i++;
    	x1=veff[i]-e;
        Rl[i]=y[0];
	    }

	  tk=t[i];
	  co=y[0];dRo=y[1];


};


void Orbital_norel::final(double r,double *R,double W){
		/* Asumiendo V->0 cuando r->0 ,segun
		  O Čertík ,J E. Pask ,J Vackář Computer Physics Communications 184 (2013) 1777–1791*/

        double k=sqrt(-2*(e-W));

	    R[0]=a*exp(-k*r);
        R[1]=	-k*R[0]*r;
};

void Orbital_norel::inicial(double r,double *R){
	 /*O Čertík ,J E. Pask ,J Vackář Computer Physics Communications 184 (2013) 1777–1791*/

	 R[0]=pow(r,l+1);
     R[1]=(l+1)*pow(r,l+1);

}




double Orbital_norel::correct_e(double *r){
	//correccion utilizando teoria de perturbaciones.
	/*O Čertík ,J E. Pask ,J Vackář Computer Physics Communications 184 (2013) 1777–1791*/

	double de=co*(dRo-a*dRi)*Z/exp(tk);


    int i=0;
    while(t[i]<tk){
    		i++;
    }

    double P1[i];

    double P2[Nt-i];


    for(int j=0;j<i;j++){
    	P1[j]=Rl[j]*Rl[j]*r[j];

    }
    for(int k=i;k<Nt;k++){
        	P2[k-i]=Rl[k]*Rl[k]*r[k];

    }






    double norm=simpson(P1,h,i-1)+simpson(P2, h, Nt-1-i);



    a=1;
    return de/(2*norm);
};

void Orbital::estimate_a(){
	a=co/ci;
};


void Orbital::resolver(double *veff,double *r,double W){
    double e_old=0.;
    double error=0.000000001;
    double de=0.;
    double inf=-Z*Z/2.;
    double sup=W;
    outward(veff,r);
    int nod=0;
    int perturbaciones=0;
    while(nodos !=n-l-1  ){
    	      if (nodos> n-l-1){
            	  sup=e;
            	  e=(sup+inf)/2;

              }
              else {
            	  inf=e;
            	  e=(sup+inf)/2;

              }
    outward(veff,r);
    nod++;
    }

while(fabs(e-e_old)>error ){

    e_old=e;
	inward(veff,r,W);
    estimate_a();
    radial();
    de=correct_e(r);
    perturbaciones++;
    e=e+de;
    outward(veff,r);
    while(nodos !=n-l-1){
        	de=de/2;
        	e=e_old+de;
            outward(veff,r);

        };
    a=1;
}

//cout<<e<<endl;
//cout<<"evaluaciones: "<<nod<<"   "<<"perturbaciones"<<perturbaciones<<endl;
inward(veff,r,W);
estimate_a();
radial();
a=1;
dens(r);

}


void Orbital_norel::radial(){

      int i=Nt-1;

      while (t[i] >= tk){
    	  Rl[i]=a*Rl[i];
          i--;

      }
 //Por convencion para que las integrales den con el signo correcto la cola de la funcion debe ser positiva.
      if(Rl[i] < 0){
    	  for(int j=0;j<Nt;j++){
    		  Rl[j]=-Rl[j];
    	  }

            }

};







void Orbital_norel::print(){
	ofstream archivo("radial.txt");


	for (int i=0;i<Nt;i++){
		archivo<<exp(t[i])/Z<<"   "<<Rl[i]<<"   "<<rho[i]<<endl;

	}
	archivo.close();
}
void Orbital_norel::dens(double *r){
	for (int i=0;i<Nt;i++){
		rho[i]=Rl[i]*Rl[i]*r[i];
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, Nt);

	    gsl_spline_init (spline, t, rho, Nt);


	    double norm=1/sqrt(gsl_spline_eval_integ (spline, t[0], t[Nt-1], acc));

        for(int i=0;i<Nt;i++){
        	Rl[i]=norm*Rl[i];
        	rho[i]=noc*Rl[i]*Rl[i]/(4*pi*r[i]*r[i]);

        };
        gsl_spline_init (spline, t, rho, Nt);

        for(int i=0;i<Nt;i++){
        	grad[i]=gsl_spline_eval_deriv(spline,t[i],acc)/r[i];
        }

        gsl_spline_free (spline);
	    gsl_interp_accel_free (acc);

};



double Orbital::operator[](int i){
	return rho[i];
};
double Orbital::operator()(int i){
	return Rl[i];
};
Orbital_rel::Orbital_rel(int principal,int angular,double num_ocupacion,int s,int z,double energy,double *r,int N) : Orbital(principal,angular,num_ocupacion,z,energy,r,N){
	if (l == 0){k=-1;}
	else if(s == 1 ){k=-l-1;}
	else{k=l;}
    Ql=new double [Nt];
} ;

Orbital_rel::~Orbital_rel(){
	delete [] Ql;
};

void Orbital_rel::inicial(double r,double *R){
	double b=sqrt(k*k-(Z/c)*(Z/c));
	R[0]=pow(r,b);
	R[1]=R[0]*c*(b+k)/Z;
}
void Orbital_rel::final(double r,double *R,double W){
	double lamb=sqrt(-2*e-e*e/(c*c));
    R[0]=exp(-lamb*r);
    R[1]=-sqrt(-(e/(e+2*c*c)))*R[0];
}

void Orbital_rel::outward(double *veff,double *r){


	  	nodos=0;

	  	Dirac Ec(e,l,k,Z,veff,t,r);
	  	double y[2];
	    double ya;
	  	Ode<Dirac> abm;
	  	int j=0;
	  	inicial(r[j],y);
	  	while(y[0]<1.0e-16){
	  	    Rl[j]=0; Ql[j]=0;
	  	    inicial(r[j],y);
	  	    j++;
 	    }
	  	for (int i=0;i<4;i++){
	  	  inicial(r[i+j],y);
	  	  abm.inicializar(Ec,y,i+j);
	  	  Rl[i+j]=y[0]; Ql[i+j]=y[1];

	  	}
	  	double x0,x1;
	    int i=4+j;
	  	x0=veff[i]-e;
	  	x1=x0;
	  	inicial(r[i],y);
	  	Rl[i]=y[0]; Ql[i]=y[1];
	  	while (x0*x1>0 ){
	      	x0=x1;
	      	ya=y[0];
	        abm.step( Ec , y ,h, i );
	          if (ya*y[0]<0){
	          	nodos++;
	      	}
	      	i++;
	      	x1=veff[i]-e;
	      	Rl[i]=y[0]; Ql[i]=y[1];
	  	    }

	  	  tk=t[i];
	  	  co=y[0];dRo=y[1];


	  };

void Orbital_rel::inward(double *veff,double *r,double W){


	  Dirac Ec(e,l,k,Z,veff,t,r);
	  double  y[2];
	  y[0]=0;
	  y[1]=0;
	  int i=Nt-1;
	  while(y[0]<1.0e-15){
	  	  final(r[i],y,W);
	      i--;
	  }
	  max=i;
	  Ode<Dirac> abm;
	  for(int l=0;l<4;l++){
	      	  final(r[i],y,W);
	      	  abm.inicializar(Ec,y,i);
	      	  Rl[i]=y[0]; Ql[i]=y[1];
	      	  i--;


	  }

      final(r[i],y,W);
      Rl[i]=y[0]; Ql[i]=y[1];

      while (t[i] > tk)
	    {

	      abm.step( Ec, y ,-h ,i  );



	    i--;
	    Rl[i]=y[0]; Ql[i]=y[1];
	    }

	  ci=y[0];dRi=y[1];
	  for(int k=max;k<Nt;k++){
	      	  Rl[k]=0.; Ql[k]=0;
	  }

};


double Orbital_rel::correct_e(double *r){

		//correccion utilizando teoria de perturbaciones.
		/*O Čertík ,J E. Pask ,J Vackář Computer Physics Communications 184 (2013) 1777–1791*/

		double de=co*(dRo-a*dRi);   // En la ecuacion Dirac no va con exp(tk)/Z.


		int i=0;
		while(t[i]<tk){
		    		i++;
		    }

		    double P1[i];

		    double P2[max-i];




	    for(int j=0;j<i;j++){
	    	P1[j]=(Rl[j]*Rl[j]+Ql[j]*Ql[j])*r[j];

	    }
	    for(int k=i;k<max;k++){
	        	P2[k-i]=(Rl[k]*Rl[k]+Ql[k]*Ql[k])*r[k];

	    }




	    double norm=simpson(P1,h,i-1)+simpson(P2, h,max-1-i);



	    a=1;
	    return c*de/norm;
	};


void Orbital_rel::radial(){

      int i=Nt-1;

      while (t[i] >= tk){
    	  Rl[i]=a*Rl[i];
          Ql[i]=a*Ql[i];
    	  i--;

      }
      if(Rl[i] < 0){
    	  for(int j=0;j<Nt;j++){
    		  Rl[j]=-Rl[j];
    	  }
      }



};
void Orbital_rel::dens(double *r){
	for (int i=0;i<Nt;i++){
		rho[i]=(Rl[i]*Rl[i]+Ql[i]*Ql[i])*r[i];
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, Nt);

	    gsl_spline_init (spline, t, rho, Nt);


	    double norm=1/sqrt(gsl_spline_eval_integ (spline, t[0], t[Nt-1], acc));

        for(int i=0;i<Nt;i++){
        	Rl[i]=norm*Rl[i];
        	Ql[i]=norm*Ql[i];
        	rho[i]=noc*(Rl[i]*Rl[i]+Ql[i]*Ql[i])/(4*pi*r[i]*r[i]);

        };

        gsl_spline_free (spline);
	    gsl_interp_accel_free (acc);

}
void Orbital_rel::print(){
	ofstream archivo("radial.txt");


	for (int i=0;i<Nt;i++){
		archivo<<exp(t[i])/Z<<"   "<<Rl[i]<<"   "<<Ql[i]<<"   "<<rho[i]<<endl;
	}
	archivo.close();

}

Orbital_spline::Orbital_spline( Orbital &a,double *t,int N){
	l=angular[a.ang()];
	rmax=t[N-1];
	e=a.energy();
	double r[N+1];
	double orb[N+1];
    r[0]=0;
    orb[0]=0;
	for (int i=1;i<N+1;i++){
    	r[i]=t[i-1];
		orb[i]=a(i-1);
    }
    R=gsl_spline_alloc (gsl_interp_cspline, N+1);
    accR=gsl_interp_accel_alloc ();
    gsl_spline_init (R, r, orb,N+1);
};
Orbital_spline::Orbital_spline( double *a,double *t,double energy,int angl,int N){
	l=angl;
	rmax=t[N-1];
	e=energy;
	double r[N+1];
	double orb[N+1];
    r[0]=0;
    orb[0]=0;
	for (int i=1;i<N+1;i++){
    	r[i]=t[i-1];
		orb[i]=a[i-1];
    }
    R=gsl_spline_alloc (gsl_interp_cspline, N+1);
    accR=gsl_interp_accel_alloc ();
    gsl_spline_init (R, r, orb,N+1);
};
Orbital_spline::Orbital_spline(const Orbital_spline &a){
	l=a.l;
	e=a.e;
	R=a.R;
	accR=a.accR;
    rmax=a.rmax;
};
Orbital_spline & Orbital_spline::operator=(const Orbital_spline &a){
	l=a.l;
	R=a.R;
	accR=a.accR;
    e=a.e;
    rmax=a.rmax;
    return *this;
};
Orbital_spline::~Orbital_spline(){
	gsl_spline_free (R);
	gsl_interp_accel_free (accR);
}
 double Orbital_spline::operator()(double x){
	return gsl_spline_eval(R,x,accR);
}
Potential_spline::Potential_spline(double *y,double *x,int N){
	v=gsl_spline_alloc (gsl_interp_cspline, N);
	accv=gsl_interp_accel_alloc ();
	gsl_spline_init (v, x, y,N);
    rmin=x[0];
}
Potential_spline::~Potential_spline(){
	gsl_spline_free (v);
	gsl_interp_accel_free (accv);
}
double simpson(double *f,double h,int N){
	double impar=0;
	double par=0;
	for (int i=0;2*i+1<N-2;i++){
		impar+=f[2*i+1];

	}
    for(int i =1;2*i<N-3;i++){
    	par+=f[2*i];
    }
    return (h/3)*(f[0]+4*impar+2*par+f[N-1]);

}
