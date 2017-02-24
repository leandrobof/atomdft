/*
 * globals.cpp
 *
 *  Created on: Nov 16, 2016
 *      Author: leandro
 */
int func (double t, const double y[], double f[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  parametros *param = ( parametros *)params;
  int l=param->l;
  double e=param->e;
  double r=exp(t);
  double veff=-Z/r;
  double dveff=Z/(r*r);
  double M=1-al*al/4*(veff-e);
  f[0] = y[1];
  f[1] = y[1]+(l*(l+1)+r*r*M*(veff-e))*y[0]-al*al/(4*M)*dveff*r*(y[1]-y[0]);
  return GSL_SUCCESS;
}



