/*
 * globals.h
 *
 *  Created on: Nov 14, 2016
 *      Author: leandro
 */
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <iostream>
#include <fstream>
#include <cstddef>
#include <vector>


using namespace std;

#ifndef GLOBALS_H_
#define GLOBALS_H_
typedef vector< double > state_type;
const double c=137.0359895;
const double pi=3.14159265358979323846;
const double al=1/c;
const double W=0;
const char angular[5]={'s','p','d','f','g'};



#endif /* GLOBALS_H_ */
