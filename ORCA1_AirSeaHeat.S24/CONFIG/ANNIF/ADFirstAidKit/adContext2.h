#ifndef ADCONTEXT_INCLUDED
#define ADCONTEXT_INCLUDED

#include "complex.h"

void adContextTgt_init(double epsilon, double seed) ;
void adContextTgt_initComplex16(char* varname, double complex *indep, double complex *indepd) ;
void adContextTgt_startConclude() ;
void adContextTgt_concludeComplex16(char* varname, double complex dep, double complex depd) ;
void adContextTgt_conclude() ;
void adContextAdj_init(double seed) ;
void adContextAdj_initComplex16(char* varname, double complex *dep, double complex *depb) ;
void adContextAdj_startConclude() ;
void adContextAdj_concludeComplex16(char* varname, double complex dep, double complex depb) ;
void adContextAdj_conclude() ;

#endif
