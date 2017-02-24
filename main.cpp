/*
 * main.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: leandro
 */

#include "Orbital.h"

using namespace std;
void gridlog(double *t,double *r,int Ngrid);
double vn(double r);
double dvn(double r);

int main(int argc,char *argv[]){
double a=atof(argv[1]);
//double b=atof(argv[2]);
Orbital H1s(1,0,a);
double h=0.001;
double dif=1.;
double error=0.00000001;
//while(dif > error){
    H1s.outward(h);
    H1s.inward(h);
    H1s.estimate_a();
    H1s.outward(h);
    dif=H1s.correct_e();
    cout<<"energia:"<<H1s.get_energy()<<"      "<<"error:"<<dif<<endl;
//}
//cout<<"energia:"<<H1s.get_energy()<<"      "<<"error:"<<dif<<endl;
return 0;
}



void gridlog(double *t,double *r,int Ngrid){
	for (int i=0;i<Ngrid;i++){
		r[i]=exp (t[i]);
	}
}
double vn(double r){
	return -Z/r;
}
double dvn(double r){
	return Z/(r*r);
}
