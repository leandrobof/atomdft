/*
 * Ode.h
 *
 *  Created on: Jul 31, 2017
 *      Author: leandro
 */

#ifndef ODE_H_
#define ODE_H_
#include <iostream>
template <typename T> class Ode {

private:

	double y[2][5];
	double predictor[2];
	double corrector[2];
    int contador;
    double dxdt[2];

public:
	Ode(void);
	virtual ~Ode();
    void inicializar(T &func,double *x,int i);
	void step(T &func,double *x,double h, int &i);
};
template <class T>  Ode<T>::Ode(void) {

    contador=0;

    }


template <class T>  void Ode<T>::inicializar(T &func,double *x, int i){
	func(x,dxdt, i);
	y[0][contador]=dxdt[0];
	y[1][contador]=dxdt[1];
	contador+=1;

}
template <class T>  Ode<T>::~Ode() {
	// TODO Auto-generated destructor stub
}
template <class T> void Ode<T>::step(T &func,double *x,double h, int &i)
{
    func(x,dxdt,i);

for (int j=0;j<2;j++){
    y[j][4]=dxdt[j];
	predictor[j]=x[j]+h/720*(y[j][4]*1901-2774*y[j][3]+y[j][2]*2616-1274*y[j][1]+251*y[j][0]);
}
if(h>0){func(predictor,dxdt,i+1);}
else{func(predictor,dxdt,i-1);};
for (int j=0;j<2;j++){
    corrector[j]=x[j]+h/720*(251*dxdt[j]+646*y[j][4]-264*y[j][3]+106*y[j][2]-19*y[j][1]);

	y[j][0]=y[j][1];
	y[j][1]=y[j][2];
	y[j][2]=y[j][3];
	y[j][3]=y[j][4];
	x[j]=(1./502)*(475*corrector[j]+27*predictor[j]);
}
}

#endif /* ODE_H_ */
