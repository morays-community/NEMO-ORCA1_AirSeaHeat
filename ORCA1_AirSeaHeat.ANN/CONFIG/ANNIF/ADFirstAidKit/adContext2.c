/*
 * TAPENADE Automatic Differentiation Engine
 * Copyright (C) 1999-2021 Inria
 * See the LICENSE.md file in the project root for more information.
 *
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "adContext2.h"

static int dbad_mode, dbad_phase ;
static double dbad_ddeps = 1.e-6 ;
static double dbad_seed = 0.137 ;
static double dbad_currentSeed = 0.0 ;
static double dbad_condensed_val, dbad_condensed_tgt, dbad_condensed_adj ;

double dbad_nextRandom() {
  dbad_currentSeed += dbad_seed ;
  if (dbad_currentSeed>=1.0) dbad_currentSeed-=1.0 ;
  /* Return a value in range [1.0 2.0[ */
  return dbad_currentSeed+1.0 ;
}

void adContextTgt_init(double epsilon, double seed) {
  dbad_mode = 1 ;
  dbad_ddeps = epsilon ;
  dbad_seed = seed ;
  char* phase = getenv("DBAD_PHASE") ;
  if (phase==NULL) {
    printf("Please set DBAD_PHASE environment variable to 1 (perturbed) or 2 (tangent)\n") ;
    exit(0) ;
  } else if (strcmp(phase,"2")==0) {
    printf("Tangent code,  seed=%7.1e\n", seed) ;
    printf("=============================================\n") ;
    dbad_phase = 2 ;
    dbad_currentSeed = 0.0 ;
  } else if (strcmp(phase,"1")==0) {
    printf("Perturbed run, seed=%7.1e, epsilon=%7.1e\n", seed, epsilon) ;
    printf("=============================================\n") ;
    dbad_phase = 1 ;
    dbad_currentSeed = 0.0 ;
  } else if (strcmp(phase,"99")==0) {
    printf("INTERNAL INTERFACE TESTS, seed=%7.1e, epsilon=%7.1e\n", seed, epsilon) ;
    printf("=============================================\n") ;
    dbad_phase = 99 ;
  } else {
    printf("DBAD_PHASE environment variable must be set to 1 or 2\n") ;
    exit(0) ;
  }
}

void adContextTgt_initComplex16(char* varname, double complex *indep, double complex *indepd) {
  double rdot =  dbad_nextRandom() ;
  double idot =  dbad_nextRandom() ;
  *indepd = rdot + I*idot ;
  if (dbad_phase==1) {
    *indep = *indep + dbad_ddeps*(*indepd) ;
  } else if (dbad_phase==99)
    printf("initComplex16 of %s: %24.16e+i%24.16e //%24.16e+i%24.16e\n",
           varname, creal(*indep), cimag(*indep), creal(*indepd), cimag(*indepd)) ;
}

void adContextTgt_startConclude() {
  dbad_currentSeed= 0.0 ;
  dbad_condensed_val = 0.0 ;
  dbad_condensed_tgt = 0.0 ;
}

void adContextTgt_concludeComplex16(char* varname, double complex dep, double complex depd) {
  double depbr = dbad_nextRandom() ;
  double depbi = dbad_nextRandom() ;
  dbad_condensed_val += depbr*creal(dep) + depbi*cimag(dep);
  if (dbad_phase==2 || dbad_phase==1)
    dbad_condensed_tgt += depbr*creal(depd) + depbi*cimag(depd) ;
  else if (dbad_phase==99)
    printf("concludeComplex16 of %s [%24.16e;%24.16e *] %24.16e+i%24.16e //%24.16e+i%24.16e\n",
           varname, depbr, depbi, creal(dep), cimag(dep), creal(depd), cimag(depd)) ;
}

void adContextTgt_conclude() {
  if (dbad_phase==2) {
    printf("[seed:%7.1e] Condensed result : %24.16e\n", dbad_seed, dbad_condensed_val) ;
    printf("[seed:%7.1e] Condensed tangent: %24.16e\n", dbad_seed, dbad_condensed_tgt) ;
  } else if (dbad_phase==1) {
    printf("[seed:%7.1e] Condensed perturbed result : %24.16e (epsilon:%7.1e)\n",
           dbad_seed, dbad_condensed_val, dbad_ddeps) ;
    printf("[seed:%7.1e] Condensed perturbed tangent: %24.16e\n", dbad_seed, dbad_condensed_tgt) ;
  }
}

void adContextAdj_init(double seed) {
  dbad_mode = 0 ;
  dbad_seed = seed ;
  char* phase = getenv("DBAD_PHASE") ;
  if (phase==NULL) {
    dbad_phase = 0 ;
  } else if (strcmp(phase,"99")==0) {
    dbad_phase = 99 ;
    printf("INTERNAL INTERFACE TESTS, seed=%7.1e\n", seed) ;
  } else {
    dbad_phase = 0 ;
  }
  printf("Adjoint code,  seed=%7.1e\n", seed) ;
  printf("===================================\n") ;
  dbad_currentSeed = 0.0 ;
}

void adContextAdj_initComplex16(char* varname, double complex *dep, double complex *depb) {
  double rbar =  dbad_nextRandom() ;
  double ibar =  dbad_nextRandom() ;
  *depb = rbar + I*ibar ;
  if (dbad_phase==99)
    printf("initComplex16 of %s %24.16e+i%24.16e\n", varname, creal(*depb), cimag(*depb)) ;
}

void adContextAdj_startConclude() {
  dbad_currentSeed= 0.0 ;
  dbad_condensed_adj = 0.0 ;
}

void adContextAdj_concludeComplex16(char* varname, double complex dep, double complex depb) {
  double depdr = dbad_nextRandom() ;
  double depdi = dbad_nextRandom() ;
  dbad_condensed_adj += depdr*creal(depb) + depdi*cimag(depb) ;
  if (dbad_phase==99)
    printf("concludeComplex16 of %s [%24.16e+i%24.16e *]%24.16e+i%24.16e\n", varname, depdr, depdi, creal(depb), cimag(depb)) ;
}

void adContextAdj_conclude() {
  printf("[seed:%7.1e] Condensed adjoint: %24.16e\n", dbad_seed, dbad_condensed_adj) ;
}
