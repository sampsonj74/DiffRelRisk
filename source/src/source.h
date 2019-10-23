#ifndef __SOURCE_H__
#define __SOURCE_H__

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define DEBUG 0
#define MAX( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define MIN( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define CHECK_MEM(obj) if (obj == NULL) {error("Memory");}
#define BADRETURN 1.0e100
#define DBL_MISS 1.0e100
#define ISMISS_GT 1.0e90
#define ISNOTMISS_LT 1.0e90
#define DEFAULTEPS 1.0e-6
#define MINVALUE0 1.0e100
#define MAXVALUE0 -1.0e100
#define NUMERICZERO 1e-16
#define MINEXPARG -100.0
#define MAXEXPARG 100.0
#define EXPMINARG 3.720076e-44
#define EXPMAXARG 2.688117e43
#define MINVALUE -100.0
#define MAXVALUE 100.0
#define MACHINE_EPS 2.220446049250313e-16

/* Integer input args */
#define IARG_NCOVARS 0
#define IARG_NSUBS 1
#define IARG_MAXITER 2
#define IARG_METHOD_RR 3
#define IARG_FIX_R 4
#define IARG_NSTRATA 5
#define IARG_NALLPARAM 6
#define IARG_NPOSSVALS 7
#define IARG_DEBUG 8

/* Double input args */
#define DARG_PROBEPS 0
#define DARG_INITSTEPSIZE 1
#define DARG_STOPTOL 2
#define DARG_F0 3

void C_optim_force(int*, double*, double*, double*, double*, double*, double*, double*,
                   double*, double*, int*, int*, double*, double*, double*, double*, 
                   int *, char **, 
                   double*, int*, double *); 

static const R_CMethodDef callMethods[] = {
  {"C_optim_force", (DL_FUNC)&C_optim_force, 21},
  {NULL, NULL, 0}
};

struct obj_struct {
  
  double *x1, *x2, *x3, *x4; 
  double *n1mx1, *n2mx2, *n3mx3, *n4mx4; /* n1-x1, n2-x2, n3-x3 */
  int  *Y;            /* outcome */
  double **Z;         /* matrix of covars. NOTE this wil be a vector of (double *) */
  int *tr;            /* treatment */
  int ncovars, ncolZ; /* Number of covars, number of columns of Z */
  int nsubs;          /* length of Y */

  /* For probabilities */
  double minProb, maxProb, minErrorProb, maxErrorProb;
  double logMinProb, logMaxProb; 

  int errorProbFlag;     /* Flag to check that computed probabilities are within [minErrorProb, maxErrorProb] */
  int maxiter;           /* maximum number of iterations in optim_force */
  double probEps;        /* Force probs to be in [probEps, 1-probEps] to prevent log errors */
  double initStepSize;   /* initial step size (similar to mult in R code) */
  double *allParam0;     /* vector of initial parameters */
  double *allParam;      /* vector of final parameters */
  int method_rr;         /* 0=OP0, 1=OP1, 2=CO, 3=COR, 4=TP */
  int fix_r;             /* same as fix.r in R code */
  int nstrata;           /* length of x1, x2, x3, x4, n1mx1, etc */
  int nallParam;         /* length of allParam */

  /* vector indices of beginning parameters for objects (adjusted for C) */
  int sp1To2s, twosp1, sp2, twosp1To3s, npm1;  
  int npossVals;         /* number of possible values to find a minimum */
  double *minParmValue;  /* minimum values for allParam */
  double *maxParmValue;  /* maximum values for allParam */
  double stopTol;        /* Stopping tolerance */
  int *parmOrder;        /* The order of the parms when finding the MLEs (adjusted for C when passed in) */

  int debug;
  double f0;             /* For testing/debugging */
  double fmin;           /* Final function (minimum) value */

  char *paramType;       /* parameter types: 
                            p3: "3"
                            p1: "a"
                            r:  "r"
                            b:  "b"
                            R1: "1"
                            R2: "2"
                            k:  "k"
                            t:  "t" */
  double minValue, maxValue; /* Default bounds for parms, not important */
};
typedef struct obj_struct OBJ;


void R_init_DiffRelRisk(DllInfo *dll)
{
    R_registerRoutines(dll, callMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}

#endif
