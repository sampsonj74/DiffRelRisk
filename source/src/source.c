/* History: Jul 15 2019 Initial coding
            Sep 18 2019 Allow initial estimates to give a non-feasible point and
                        iterate until good estimates are found
*/

#include "./source.h"


static void print_dVec(vec, n, name)
double *vec;
int n;
char name[10];
{
  int i;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %g ", vec[i]);
  }
  Rprintf("\n \n");
}
/*
static void print_iVec(vec, n, name)
int *vec;
int n;
char name[10];
{
  int i;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %d ", vec[i]);
  }
  Rprintf("\n \n");
}
*/
void print_dMat(mat, nr, nc, name)
double **mat;
int nr, nc;
char name[10];
{
  int i, j;
  Rprintf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) Rprintf(" %g ", mat[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}


static double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} 

/* Function to allocate a double matrix */
double ** dMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag;
double initVal;
{
  double **mat, **ptr;
  int i;

  mat = (double **) malloc(nrow*sizeof(double *));
  CHECK_MEM(mat);
  if (ncol > 0) {
    for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = dVec_alloc(ncol, initFlag, initVal);
  } else {
    for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = NULL;
  }

  return(mat);

} /* END: dMat_alloc */



/*
static int * iVec_alloc(n, initFlag, initVal)
int n, initFlag, initVal;
{
  int i, *ret, *p;

  ret = (int *) malloc(n*sizeof(int));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} 
*/

/*
void colsumsVecNoMiss(vec, nr, nc, ret)
double *vec, *ret;
int nr, nc;
{
  int i, j;
  double sum, *pv, *pr;

  pv = vec;
  for (i=0, pr=ret; i<nc; i++, pr++) {
    sum = 0.0;
    for (j=0; j<nr; j++) sum += *pv++;
    *pr = sum;
  }
   
}

void C_colsums(vec, pnr, pnc, ret)
double *vec, *ret;
int *pnr, *pnc;
{

  colsumsVecNoMiss(vec, *pnr, *pnc, ret);
  return;

} 
*/

static void obj_init(obj)
OBJ *obj;
{

  obj->x1            = NULL;
  obj->x2            = NULL;
  obj->x3            = NULL; 
  obj->x4            = NULL;
  obj->n1mx1         = NULL; 
  obj->n2mx2         = NULL;
  obj->n3mx3         = NULL;
  obj->n4mx4         = NULL;
  obj->Y             = NULL;
  obj->Z             = NULL;
  obj->tr            = NULL;
  obj->ncovars       = 0; 
  obj->ncolZ         = 0;
  obj->nsubs         = 0;
  obj->minProb       = DEFAULTEPS;
  obj->maxProb       = 1.0 - DEFAULTEPS;
  obj->minErrorProb  = -DEFAULTEPS;
  obj->maxErrorProb  = 1.0 + DEFAULTEPS;
  obj->logMinProb    = log(obj->minProb);
  obj->logMaxProb    = log(obj->maxProb);
  obj->errorProbFlag = 1;
  obj->maxiter       = 10000;
  obj->probEps       = DEFAULTEPS;
  obj->initStepSize  = 1.0;
  obj->allParam0     = NULL;
  obj->allParam      = NULL;
  obj->method_rr     = -1;
  obj->fix_r         = -1;
  obj->nstrata       = 0;
  obj->nallParam     = 0;
  obj->sp1To2s       = -1;
  obj->twosp1        = -1; 
  obj->sp2           = -1;
  obj->twosp1To3s    = -1; 
  obj->npm1          = -1;
  obj->npossVals     = 100;
  obj->minParmValue  = NULL;
  obj->maxParmValue  = NULL;
  obj->stopTol       = 1.0e-3;
  obj->parmOrder     = NULL;

  obj->debug         = 0;
  obj->f0            = 0.0;
  obj->fmin          = BADRETURN;

  obj->minValue      = MINVALUE;
  obj->maxValue      = MAXVALUE;
  obj->paramType     = NULL;
 
} /* obj_init */

static void obj_free(obj)
OBJ *obj;
{

  if (obj->Z) free(obj->Z);

} /* END: obj_free */ 

static double ** getCovarMat(Zvec, nr, nc)
double *Zvec;
int nr, nc;
{
  double **ret;
  int i;

  /* Get a vector of double * */
  ret = dMat_alloc(nr, 0, 0, 0.0);
  for (i=0; i<nr; i++) ret[i] = &Zvec[nc*i];

  return(ret);

} /* END: getCovarMat */

static void obj_setup(obj, iargs, dargs, x1, x2, x3, x4, n1x1, n2x2, n3x3, n4x4, Y, tr, Z, 
                      allParam0, allParam, minParmValue, maxParmValue, parmOrder, parmType)
OBJ *obj;
int *iargs, *Y, *tr, *parmOrder;
double *dargs, *x1, *x2, *x3, *x4, *n1x1, *n2x2, *n3x3, *n4x4, *Z, *allParam0, *allParam,
       *minParmValue, *maxParmValue;
char *parmType;
{
  int method, s;
  double eps;

  method = iargs[IARG_METHOD_RR];
  s      = iargs[IARG_NSTRATA];
  
  if ((method < 0) || (method > 4)) error("ERROR with method");
  if (s < 0) error("ERROR with nstrata");
  obj->errorProbFlag = 1;

  if ((method < 2) || (method == 4)) {
    obj->x1            = x1;
    obj->x2            = x2;
    obj->x3            = x3; 
    obj->n1mx1         = n1x1; 
    obj->n2mx2         = n2x2;
    obj->n3mx3         = n3x3;
    if (method == 4) {
      obj->x4          = x4;
      obj->n4mx4       = n4x4;
    }
  }
  if ((method == 2) || (method == 3)) {
    obj->Y             = Y;
    obj->tr            = tr;
    obj->ncovars       = iargs[IARG_NCOVARS]; 
    obj->ncolZ         = 1 + obj->ncovars;
    obj->nsubs         = iargs[IARG_NSUBS];
    if (obj->ncovars) {
      obj->Z = getCovarMat(Z, obj->nsubs, obj->ncolZ);
    }
    obj->errorProbFlag = 0;
  }

  eps                = dargs[DARG_PROBEPS];
  obj->minProb       = eps;
  obj->maxProb       = 1.0 - eps;
  obj->minErrorProb  = 0.0;
  obj->maxErrorProb  = 1.0;
  obj->logMinProb    = log(obj->minProb);
  obj->logMaxProb    = log(obj->maxProb);
  obj->maxiter       = iargs[IARG_MAXITER];
  obj->probEps       = eps;
  obj->initStepSize  = dargs[DARG_INITSTEPSIZE];
  obj->allParam0     = allParam0;
  obj->allParam      = allParam;
  obj->nallParam     = iargs[IARG_NALLPARAM];
  obj->method_rr     = method;
  obj->fix_r         = iargs[IARG_FIX_R];
  obj->nstrata       = s;
  
  obj->sp1To2s       = s + 1   - 1;
  obj->twosp1        = 2*s + 1 - 1; 
  obj->sp2           = s + 2   - 1;
  obj->twosp1To3s    = 2*s + 1 - 1; 
  obj->npm1          = obj->nallParam - 1;
  obj->npossVals     = iargs[IARG_NPOSSVALS];
  obj->minParmValue  = minParmValue;
  obj->maxParmValue  = maxParmValue;
  obj->stopTol       = dargs[DARG_STOPTOL];
  obj->parmOrder     = parmOrder;

  obj->debug         = iargs[IARG_DEBUG];

  if (obj->nallParam < 3) error("ERROR with nallParam");
  if (obj->npossVals < 2) error("ERROR with npossVals");
  if (eps <= 0.0) error("ERROR with eps");
  obj->f0 = dargs[DARG_F0];
  obj->paramType = parmType;
  

} /* END: obj_setup */

static void copy_dVec(v, n, ret)
double *v, *ret;
int n;
{
  int i;
  
  for (i=0; i<n; i++) ret[i] = v[i];

} /* END: copy_dVec */

static double dotProd(v1, v2, n)
double *v1, *v2;
int n;
{
  int i;
  double ret=0.0;

  for (i=0; i<n; i++) ret += v1[i]*v2[i];

  return(ret);

} /* dotProd */

static double dotProdDropIndex(v1, v2, n, drop)
double *v1, *v2;
int n, drop;
{
  int i;
  double ret=0.0;

  for (i=0; i<n; i++) {
    if (i != drop) ret += v1[i]*v2[i];
  }

  return(ret);

} /* dotProdDropIndex */


static int isOutside_dRange(x, minValue, maxValue)
double x, minValue, maxValue;
{
  if ((x < minValue) || (x > maxValue)) return(1);
 
  return(0);

} /* check_dRange */

static double bndLogProb(p, minp, maxp, logMinp, logMaxp)
double p, minp, maxp, logMinp, logMaxp;
{
  double ret;

  if (p < minp) {
    ret = logMinp;
  } else if (p > maxp) {
    ret = logMaxp;
  } else {
    ret = log(p);
  }

  return(ret);

} /* END: bndLogProb */



static double fL_OP0(pa, p3, r, n, obj)
double *pa, *p3, r;
int n;
OBJ *obj;
{
  double ret=0.0, tmp, p1, p2, rOver2, pai, p3i, minProb, maxProb, minValue, maxValue;
  double *x1, *x2, *x3, *n1mx1, *n2mx2, *n3mx3, logMinProb, logMaxProb;
  int i, setToMissIfError;

  minProb          = obj->minProb;
  maxProb          = obj->maxProb;
  logMinProb       = obj->logMinProb;
  logMaxProb       = obj->logMaxProb;
  minValue         = obj->minErrorProb;
  maxValue         = obj->maxErrorProb;
  setToMissIfError = obj->errorProbFlag;
  x1               = obj->x1;
  x2               = obj->x2;
  x3               = obj->x3; 
  n1mx1            = obj->n1mx1;
  n2mx2            = obj->n2mx2; 
  n3mx3            = obj->n3mx3;  

  rOver2 = 0.5*r;
  for (i=0; i<n; i++) {
    pai = pa[i];
    p3i = p3[i]; 
    tmp = rOver2*p3i;
    p1  = pai - tmp;
    p2  = pai + tmp;
    if (setToMissIfError) { 
      if (isOutside_dRange(p1, minValue, maxValue)) return(DBL_MISS);
      if (isOutside_dRange(p2, minValue, maxValue)) return(DBL_MISS);
    }
    tmp = x3[i]*bndLogProb(p3i, minProb, maxProb, logMinProb, logMaxProb) + 
           n3mx3[i]*bndLogProb(1.0 - p3i, minProb, maxProb, logMinProb, logMaxProb) +
           x1[i]*bndLogProb(p1, minProb, maxProb, logMinProb, logMaxProb) + 
           n1mx1[i]*bndLogProb(1.0 - p1, minProb, maxProb, logMinProb, logMaxProb) +
           x2[i]*bndLogProb(p2, minProb, maxProb, logMinProb, logMaxProb) + 
           n2mx2[i]*bndLogProb(1.0 - p2, minProb, maxProb, logMinProb, logMaxProb);
    ret += tmp;

  }
  ret = -ret;

  return(ret);

} /* END: fL_OP0 */

static double fL_OP1(t, p3, r, n, obj)
double t, *p3, r;
int n;
OBJ *obj;
{
  double ret=0.0, tmp, tmp2, p1, p2, rOver2, p3i, minProb, maxProb, minValue, maxValue;
  double *x1, *x2, *x3, *n1mx1, *n2mx2, *n3mx3, logMinProb, logMaxProb;
  int i, setToMissIfError;

  minProb          = obj->minProb;
  maxProb          = obj->maxProb;
  logMinProb       = obj->logMinProb;
  logMaxProb       = obj->logMaxProb;
  minValue         = obj->minErrorProb;
  maxValue         = obj->maxErrorProb;
  setToMissIfError = obj->errorProbFlag;
  x1               = obj->x1;
  x2               = obj->x2;
  x3               = obj->x3; 
  n1mx1            = obj->n1mx1;
  n2mx2            = obj->n2mx2; 
  n3mx3            = obj->n3mx3;  

  rOver2 = 0.5*r;
  for (i=0; i<n; i++) {
      p3i  = p3[i]; 
      tmp  = rOver2*p3i;
      tmp2 = t*p3i; 
      p1   = tmp2 - tmp;
      p2   = tmp2 + tmp;
      if (setToMissIfError) { 
        if (isOutside_dRange(p1, minValue, maxValue)) return(DBL_MISS);
        if (isOutside_dRange(p2, minValue, maxValue)) return(DBL_MISS);
      }
      ret += x3[i]*bndLogProb(p3i, minProb, maxProb, logMinProb, logMaxProb) + 
             n3mx3[i]*bndLogProb(1.0 - p3i, minProb, maxProb, logMinProb, logMaxProb) +
             x1[i]*bndLogProb(p1, minProb, maxProb, logMinProb, logMaxProb) + 
             n1mx1[i]*bndLogProb(1.0 - p1, minProb, maxProb, logMinProb, logMaxProb) +
             x2[i]*bndLogProb(p2, minProb, maxProb, logMinProb, logMaxProb) + 
             n2mx2[i]*bndLogProb(1.0 - p2, minProb, maxProb, logMinProb, logMaxProb);
  }
  ret = -ret;

  return(ret);

} /* END: fL_OP1 */

static double fL_CO(r, k, b, obj)
double r, k, *b;
OBJ *obj;
{
  int i, nsubs, ncovars, *Y, ncZ, *tr, setToMissIfError;
  double **Z, ret, expb, kmrOver2, kprOver2, vec, tri, p;
  double minProb, maxProb, logMinProb, logMaxProb, minValue, maxValue;

  ret              = 0.0;
  Y                = obj->Y;
  Z                = obj->Z;
  ncovars          = obj->ncovars;
  ncZ              = obj->ncolZ;
  nsubs            = obj->nsubs;
  kmrOver2         = k - 0.5*r;
  kprOver2         = k + 0.5*r;
  tr               = obj->tr;
  minProb          = obj->minProb;
  maxProb          = obj->maxProb;
  logMinProb       = obj->logMinProb;
  logMaxProb       = obj->logMaxProb;
  minValue         = obj->minErrorProb;
  maxValue         = obj->maxErrorProb;
  setToMissIfError = obj->errorProbFlag;

  if (!ncovars) {
    expb = exp(*b);
  } else {
    Z = obj->Z;
  }

  for (i=0; i<nsubs; i++) {
    tri = tr[i];
    if (tri == 1) {
      vec = kmrOver2;
    } else if (tri == 2) {
      vec = kprOver2;
    } else {
      vec = 1.0;
    }
    if (ncovars) {
      p = exp(dotProd(Z[i], b, ncZ))*vec;
    } else {
      p = expb*vec;
    }
    
    if (setToMissIfError) { 
      if (isOutside_dRange(p, minValue, maxValue)) return(DBL_MISS);
    }

    if (Y[i]) {
      ret += bndLogProb(p, minProb, maxProb, logMinProb, logMaxProb);
    } else {
      ret += bndLogProb(1.0 - p, minProb, maxProb, logMinProb, logMaxProb);
    }
  }
  ret = -ret;

  return(ret);

} /* END: fL_CO */

static double fL_COR(R1, R2, b, obj)
double R1, R2, *b;
OBJ *obj;
{
  double ret;

  ret = fL_CO(R2-R1, 0.5*(R1 + R2), b, obj);

  return(ret);

} /* END: fL_COR */

static double fL_TP(t, r, p01, p02, n, obj)
double *t, r, *p01, *p02;
int n;
OBJ *obj;
{
  double ret=0.0, tmp, p1, p2, p01i, p02i, ti, minp, maxp, minValue, maxValue;
  double *x1, *x2, *x3, *x4, *n1mx1, *n2mx2, *n3mx3, *n4mx4, logMinp, logMaxp;
  int i, setToMissIfError;

  minp             = obj->minProb;
  maxp             = obj->maxProb;
  logMinp          = obj->logMinProb;
  logMaxp          = obj->logMaxProb;
  minValue         = obj->minErrorProb;
  maxValue         = obj->maxErrorProb;
  setToMissIfError = obj->errorProbFlag;
  x1               = obj->x1;
  x2               = obj->x2;
  x3               = obj->x3; 
  x4               = obj->x4;
  n1mx1            = obj->n1mx1;
  n2mx2            = obj->n2mx2; 
  n3mx3            = obj->n3mx3;
  n4mx4            = obj->n4mx4;
  
  for (i=0; i<n; i++) {
    ti   = t[i];
    p01i = p01[i];
    p02i = p02[i];
    p1   = 0.5*p01i*(ti - r);
    p2   = 0.5*p02i*(ti + r);
   
    if (setToMissIfError) { 
      if (isOutside_dRange(p1, minValue, maxValue)) return(DBL_MISS);
      if (isOutside_dRange(p2, minValue, maxValue)) return(DBL_MISS);
    }

    tmp = x4[i]*bndLogProb(p02i, minp, maxp, logMinp, logMaxp) + 
           n4mx4[i]*bndLogProb(1.0 - p02i, minp, maxp, logMinp, logMaxp) +
           x3[i]*bndLogProb(p01i, minp, maxp, logMinp, logMaxp) + 
           n3mx3[i]*bndLogProb(1.0 - p01i, minp, maxp, logMinp, logMaxp) +
           x1[i]*bndLogProb(p1, minp, maxp, logMinp, logMaxp) + 
           n1mx1[i]*bndLogProb(1.0 - p1, minp, maxp, logMinp, logMaxp) +
           x2[i]*bndLogProb(p2, minp, maxp, logMinp, logMaxp) + 
           n2mx2[i]*bndLogProb(1.0 - p2, minp, maxp, logMinp, logMaxp);
    ret += tmp;
  }
  ret = -ret;

  return(ret);

} /* END: fL_TP */

static double fL(allParam, obj)
double *allParam;
OBJ *obj;
{
  int method_rr=obj->method_rr, nstrata=obj->nstrata;
  double ret;

  /* All the objects pa, p3, etc are contained in allParam, so just pass in pointers
     to the first element of the object in vector allParam. Note that nstrata will
     determine how much of the vector will be used for each object. 
     Corresponding R code:
       func.vals <- fL.OP.mat(NA, allParamMat[, oneTos, drop=FALSE], allParamMat[, sp1To2s, drop=FALSE], 
                                allParamMat[, twosp1], x1mat, x2mat, x3mat, nx1, nx2, nx3)
      } else if (OP1) {
        func.vals <- fL.OP.mat(allParamMat[, 1], NA, allParamMat[, twoTosp1, drop=FALSE], 
                     allParamMat[, sp2], x1mat, x2mat, x3mat, nx1, nx2, nx3)
      } else if (CO) {
        func.vals <- fL.CO.mat(allParamMat[, nallParam], allParamMat[, 1], allParamMat[, -vec1, drop=FALSE], 
                               Y, Z, tr1, tr2, oneMinusY)
      } else if (COR) {
        func.vals <- fL.COR.mat(allParamMat[, nallParam], allParamMat[, 1], allParamMat[, -vec1, drop=FALSE], 
                               Y, Z, tr1, tr2, oneMinusY)
      } else if (TP) {
        func.vals <- try(fL.TP.mat(allParamMat[, oneTos, drop=FALSE], allParamMat[, nallParam],
                               allParamMat[, sp1To2s, drop=FALSE], allParamMat[, twosp1To3s, drop=FALSE],
                               x1, x2, x3, x4, n1, n2, n3, n4))
      } 
  */
  if (!method_rr) {
    ret = fL_OP0(allParam, &allParam[obj->sp1To2s], allParam[obj->twosp1], nstrata, obj);
  } else if (method_rr == 1) {
    ret = fL_OP1(allParam[0], &allParam[1], allParam[obj->sp2], nstrata, obj);
  } else if (method_rr == 2) {
    ret = fL_CO(allParam[obj->npm1], allParam[0], &allParam[1], obj);
  } else if (method_rr == 3) {
    ret = fL_COR(allParam[obj->npm1], allParam[0], &allParam[1], obj);
  } else if (method_rr == 4) {
    ret = fL_TP(allParam, allParam[obj->npm1], &allParam[obj->sp1To2s], 
                &allParam[obj->twosp1To3s], nstrata, obj);
  }

  return(ret);

} /* END: fL */

static void optimRange_OP0(obj, curIndex, minv, maxv)
OBJ *obj;
int curIndex;
double *minv, *maxv;
{
  int i, s=obj->nstrata, np=obj->nallParam;
  double *parms=obj->allParam, r, pa, p3, tmp, tmp1, tmp2;
  char *ptype=obj->paramType, type;


  /* parms are c(pa, p3, r) */
  type = ptype[curIndex];


  if (type == 'a') {  /* pa */
    r     = fabs(parms[np-1]);
    p3    = 0.5*parms[curIndex+s]; 
    tmp1  = r*p3;
    *minv = MAX(tmp1, obj->minProb);
    *maxv = MIN(1.0-tmp1, obj->maxProb);
  } else if (type == 'r') {
    for (i=0; i<s; i++) {
      pa   = parms[i];
      p3   = parms[i+s];
      /* update min.v */
      tmp1 = -2.0*(1.0-pa)/p3;
      tmp2 = -2.0*pa/p3;
      tmp  = MAX(tmp1, tmp2);
      if (!i) {
        *minv = tmp;
      } else {
        *minv = MAX(*minv, tmp);
      }
      /* Update max.v */
      tmp  = MIN(-tmp1, -tmp2);
      if (!i) {
        *maxv = tmp;
      } else {
        *maxv = MIN(*maxv, tmp);
      }
    }
  } else if (type == '3') {  /* p3 */
    i     = curIndex - s;
    if (i < 0) error("INTERNAL CODING ERROR in optimRange_OP0: i < 0");
    pa    = parms[i];
    r     = fabs(parms[np-1]);
    *minv = obj->minProb;
    if (r > NUMERICZERO) {
      tmp1  = 2.0*(1.0-pa)/r;
      tmp2  = 2.0*pa/r;
      tmp   = MIN(tmp1, tmp2);
      *maxv = MIN(tmp, obj->maxProb);
    } else {
      *maxv = obj->maxProb;
    }
  }
/*
Rprintf("index=%d, type=%c, minv=%g, maxv=%g\n", curIndex, type, *minv, *maxv);
*/

} /* END: optimRange_OP0 */

static void optimRange_OP1(obj, curIndex, minv, maxv)
OBJ *obj;
int curIndex;
double *minv, *maxv;
{
  int i, s=obj->nstrata, np=obj->nallParam;
  double *parms=obj->allParam, t, r, p3, tmp, tmp1, tmp2;
  char *ptype=obj->paramType, type;

  /* parms are c(t, p3, r) */
  type = ptype[curIndex];
  if (type == 't') {  /* t */
    r     = parms[np-1];
    p3    = parms[curIndex+1]; 
    *minv = 0.5*fabs(r);
    tmp   = r*p3/2.0;
    *maxv = MIN((1.0+tmp)/p3, (1.0-tmp)/p3);
  } else if (type == 'r') {
    t = parms[0];
    for (i=0; i<s; i++) {
      p3   = parms[i+1];
      /* update min.v */
      tmp1 = -2.0*(1.0-t*p3)/p3;
      tmp2 = -2.0*t;
      tmp  = MAX(tmp1, tmp2);
      if (!i) {
        *minv = tmp;
      } else {
        *minv = MAX(*minv, tmp);
      }
      /* Update max.v */
      tmp  = MIN(-tmp1, -tmp2);
      if (!i) {
        *maxv = tmp;
      } else {
        *maxv = MIN(*maxv, tmp);
      }
    }
  } else if (type == '3') {  /* p3 */
    t     = parms[0];
    r     = fabs(parms[np-1]);
    *minv = obj->minProb;
    tmp   = 1.0/(t + 0.5*r);
    *maxv = MIN(tmp, obj->maxProb);
  }

} /* END: optimRange_OP1 */

static double myExp(x)
double x;
{
  double ret;

  if (x > MAXEXPARG) {
    ret = EXPMAXARG;
  } else if (x < MINEXPARG) {
    ret = EXPMINARG;
  } else {
    ret = exp(x);
  }

  return(ret);

} /* END: myExp */

static double getMaxExpZb(Z, b, nr, nc)
double **Z, *b;
int nr, nc;
{
  int i;
  double ret=-1.0, tmp;

  if (Z) {
    for (i=0; i<nr; i++) {
      tmp = myExp(dotProd(Z[i], b, nc));
      if (tmp > ret) ret = tmp;
    }
  } else {
    /* Intercept only */
    ret = exp(b[0]);
  }

  return(ret);

} /* END: getMaxExpZb */

static void optimRange_CO(obj, curIndex, minv, maxv)
OBJ *obj;
int curIndex;
double *minv, *maxv;
{
  int i, np=obj->nallParam, ncZ=obj->ncolZ, nsubs=obj->nsubs, *tr, tri, drop;
  double *parms=obj->allParam, r, tmp, tmp1, tmp2, kmr2, kpr2, **Z=obj->Z, *b, ppp, minp;
  char *ptype=obj->paramType, type;
  double k, minv0, maxv0, *Zi;

  /* parms are c(k, b, r) */
  type = ptype[curIndex];

  if (type == 'k') {
    r    = parms[np-1];
    b    = &parms[1];
    tmp1 = fabs(r);
    tmp  = 1.0001*0.5*tmp1;
    if (tmp1 < obj->minProb) tmp += obj->minProb;
    *minv = tmp;
    tmp   = getMaxExpZb(Z, b, nsubs, ncZ);
    if (tmp < NUMERICZERO) tmp = NUMERICZERO;
    tmp1  = 1.0/tmp;
    tmp2  = 0.5*r;
    tmp   = MIN(tmp1-tmp2, tmp1+tmp2);
    *maxv = MIN(tmp, 3.0);
  } else if (type == 'b') {
    if (Z) {
      minv0 = -1.0e10;
      maxv0 = 1.0e10;
      tr    = obj->tr;
      minp  = obj->minProb;
      r     = parms[np-1];
      k     = parms[0];
      kmr2  = k - 0.5*r;
      kpr2  = k + 0.5*r;
      b     = &parms[1];
      drop  = curIndex - 1;
      for (i=0; i<nsubs; i++) {
        tri = tr[i];
        Zi  = Z[i];
        tmp = dotProdDropIndex(Zi, b, ncZ, drop);
        ppp = myExp(tmp); 
        if (tri == 1) {
          ppp *= kmr2;
        } else if (tri == 2) {
          ppp *= kpr2;
        }
        if (ppp < minp) ppp = minp; 
        tmp  = Zi[drop];
        tmp1 = log(1.0/ppp)/tmp;
        if (tmp < 0.0) {
          minv0 = MAX(minv0, tmp1);
        } else if (tmp > 0.0) {
          maxv0 = MIN(maxv0, tmp1);   
        }
      }
      if (minv0 > maxv0) {
        minv0 = obj->minValue;
        maxv0 = obj->maxValue;
      }
      *minv = minv0;
      *maxv = maxv0;
    } else {
      *minv = obj->minValue;
      *maxv = obj->maxValue;
    }
  } else if (type == 'r') {
    minv0 = -5.0;
    maxv0 = 5.0;
    k     = parms[0];
    b     = &parms[1];
    tr    = obj->tr;
    kmr2  = 2.0*k;
    for (i=0; i<nsubs; i++) {
      tri = tr[i];
      if (!tri) {
        tmp1 = -5.0;
        tmp2 = 5.0;
      } else {
        Zi  = Z[i];
        tmp = myExp(dotProd(Zi, b, ncZ));
        tmp = 2.0*(1.0/tmp - k);
        if (tri == 2) {
          tmp1 = -kmr2;
          tmp2 = tmp;
        } else {
          tmp1 = -tmp;
          tmp2 = kmr2;
        }
      }
      minv0 = MIN(minv0, tmp1);
      maxv0 = MIN(maxv0, tmp2);
    }  
    if (minv0 > maxv0) {
      minv0 = obj->minValue;
      maxv0 = obj->maxValue;
    }
    *minv = minv0;
    *maxv = maxv0;
  }

} /* END: optimRange_CO */

static double getMaxOneOverAbsVec(b, n)
double *b;
int n;
{
  int i;
  double ret, tmp;

  tmp = fabs(b[0]);
  if (tmp < NUMERICZERO) tmp = NUMERICZERO;
  ret = 1.0/tmp;
  if (n > 1) {
    for (i=1; i<n; i++) {
      tmp = fabs(b[i]);
      if (tmp < NUMERICZERO) tmp = NUMERICZERO;  
      tmp = 1.0/tmp;
      if (tmp > ret) ret = tmp;
    }
  }

  return(ret);

} /* END: getMaxOneOverVec */

static void optimRange_COR(obj, curIndex, minv, maxv)
OBJ *obj;
int curIndex;
double *minv, *maxv;
{
  int i, np=obj->nallParam, ncZ=obj->ncolZ, nsubs=obj->nsubs, *tr, tri, drop;
  double *parms=obj->allParam, tmp, tmp1, tmp2, **Z=obj->Z, *b, ppp;
  char *ptype=obj->paramType, type;
  double minv0, maxv0, *Zi, R1, R2, minp=obj->minProb;

  /* parms are c(R2, b, R1) */
  type = ptype[curIndex];

  if ((type == '1') || (type == '2')) {
    *minv = obj->minProb;
    b     = &parms[1];
    *maxv = getMaxOneOverAbsVec(b, ncZ);
  } else if (type == 'b') {
    if (Z) {
      minv0 = -1.0e10;
      maxv0 = 1.0e10;
      R2    = parms[np-1];
      R1    = parms[0];
      b     = &parms[1];
      tr    = obj->tr;
      drop  = curIndex - 1;
      for (i=0; i<nsubs; i++) {
        tri = tr[i];
        Zi  = Z[i];
        tmp = dotProdDropIndex(Zi, b, ncZ, drop);
        ppp = myExp(tmp); 
        if (tri == 1) {
          ppp *= R1;
        } else if (tri == 2) {
          ppp *= R2;
        }
        if (ppp < minp) ppp = minp; 
        tmp  = Zi[drop];
        tmp1 = log(1.0/ppp)/tmp;
        if (tmp < 0.0) {
          minv0 = MAX(minv0, tmp1);
        } else if (tmp > 0.0) {
          maxv0 = MIN(maxv0, tmp1);   
        }
      }
      if (minv0 > maxv0) {
        minv0 = obj->minValue;
        maxv0 = obj->maxValue;
      }
      *minv = minv0;
      *maxv = maxv0;
    } else {
      *minv = obj->minValue;
      *maxv = obj->maxValue;
    }
  } 

} /* END: optimRange_COR */

static void optimRange_TP(obj, curIndex, minv, maxv)
OBJ *obj;
int curIndex;
double *minv, *maxv;
{
  if (curIndex < obj->nallParam - 1) {
    *minv = obj->minProb;
    *maxv = obj->maxProb;
  } else {
    *minv = obj->minValue;
    *maxv = obj->maxValue;
  }

} /* END: optimRange_TP */

static void optimRange(obj, curIndex, minv, maxv)
OBJ *obj;
int curIndex;
double *minv, *maxv;
{
  int method=obj->method_rr;

  if (!method) {
    optimRange_OP0(obj, curIndex, minv, maxv);
  } else if (method == 1) {
    optimRange_OP1(obj, curIndex, minv, maxv);
  } else if (method == 2) {
    optimRange_CO(obj, curIndex, minv, maxv);
  } else if (method == 3) {
    optimRange_COR(obj, curIndex, minv, maxv);
  } else if (method == 4) {
    optimRange_TP(obj, curIndex, minv, maxv);
  } else {
    error("INTERNAL CODING ERROR in optimRange");
  }

} /* END: optimRange */

static double getStepAndBounds(method, parmValue, parmType, mult, npossVals, minv, maxv, lower, upper)
double parmValue, mult, minv, maxv, *lower, *upper;
int method, npossVals;
char parmType;
{
  double h, tmp, mult2;

  if (method < 2) {
    tmp    = 0.1*mult;
    *lower = MAX(parmValue - tmp, minv);
    *upper = MIN(parmValue + tmp, maxv); 
  } else if (method == 2) {
    tmp    = 2.0*mult;
    *lower = MAX(parmValue - tmp, minv);
    *upper = MIN(parmValue + tmp, maxv); 
  } else if (method == 3) {
    if ((parmType == '1') || (parmType == '2')) {
      mult2 = 0.1;
    } else {
      mult2 = 2.0;
    }
    tmp    = mult2*mult;
    *lower = MAX(parmValue - tmp, minv);
    *upper = MIN(parmValue + tmp, maxv); 
  } else {
    tmp = 0.1*mult;
    if (parmType != 'r') {
      *lower = MAX(parmValue - tmp, minv);
      *upper = MIN(parmValue + tmp, maxv); 
    } else {
      *lower = parmValue - tmp;
      *upper = parmValue + tmp;  
    }
  }

  if (*lower > *upper) {
    tmp    = *lower;
    *lower = *upper;
    *upper = tmp;
  }
  h = (*upper - *lower)/(npossVals - 1); 
  if (h < 0.0) {
    error("INTERNAL CODING ERROR in getStepAndBounds: h < 0");
  }
  return(h);

} /* END: getStepAndBounds */


static double findMinOrig(obj, allParam, curIndex, h, N, lower, upper, fmin, midj, isChange)
OBJ *obj;
double *allParam, h, lower, upper, fmin;
int curIndex, N, *isChange, midj;
{
  /* Vector allParam will be updated */
  int i, stopFlag;
  double possVal, allParamCur, fval0, fval, midpt, minx;

  allParamCur = allParam[curIndex];
  *isChange   = 0;
  stopFlag    = 0;
  possVal     = allParamCur;
  fval0       = fmin;
  midpt       = lower + midj*h;
  possVal     = midpt;

  /* First try to the left */
  for (i=0; i<N; i++) {
    allParam[curIndex] = possVal;
    fval               = fL(allParam, obj);
    if ((fval < ISNOTMISS_LT) && (fval < fval0 - MACHINE_EPS)) {
      stopFlag = 1;
      minx     = possVal;
      fval0    = fval;
    } else if (stopFlag && (fval > fval0 + MACHINE_EPS)) {
      break;
    } 
    possVal = possVal - h; 
  }

  if (!stopFlag) {
    /* No better point found to the left, so try to the right */
    possVal = midpt;

    for (i=0; i<N; i++) {
      possVal            = possVal + h;
      allParam[curIndex] = possVal;
      fval               = fL(allParam, obj);
      if ((fval < ISNOTMISS_LT) && (fval < fval0 - MACHINE_EPS)) {
        stopFlag = 1;
        minx     = possVal;
        fval0    = fval;
      } else if (stopFlag && (fval > fval0 + MACHINE_EPS)) {
        break;
      }
    }
  }
  if (stopFlag) {
    *isChange          = 1;
    allParam[curIndex] = minx;
  } else {
    allParam[curIndex] = allParamCur;
  }

  return(fval0);

} /* END: findMinOrig */

static double findMin2(obj, allParam, curIndex, h, N, lower, upper, minv, maxv, fmin, isChange)
OBJ *obj;
double *allParam, h, lower, upper, minv, maxv, fmin;
int curIndex, N, *isChange;
{
  /* Vector allParam will be updated */
  int i, stopFlag;
  double possVal, allParamCur, fval0, fval, minx;

  allParamCur = allParam[curIndex];
  *isChange   = 0;
  fval0       = DBL_MISS; 
  stopFlag    = 0;
  possVal     = lower;

  for (i=0; i<N; i++) {
    allParam[curIndex] = possVal;
    fval               = fL(allParam, obj);
    if ((fval < ISNOTMISS_LT) && (fval < fval0)) {
      stopFlag = 1;
      minx     = possVal;
      fval0    = fval;
    } 
    possVal +=h;
  }

  if ((allParamCur >= minv - MACHINE_EPS) && (allParamCur <= maxv + MACHINE_EPS)) {
    possVal            = allParamCur;
    allParam[curIndex] = possVal;
    fval               = fL(allParam, obj);
    if ((fval < ISNOTMISS_LT) && (fval < fval0)) {
      stopFlag = 1;
      minx     = possVal;
      fval0    = fval;
    } 
  }
 
  if (stopFlag) {
    *isChange          = 1;
    allParam[curIndex] = minx;
  } else {
    /* This must be done to make sure that all parms are within a certain range
       in future iterations */
    allParam[curIndex] = lower;
  }

  return(fval0);

} /* END: findMin2 */


static int optim_force_main(obj)
OBJ *obj;
{
  int i, nattempt, curIndex, nparam, ok, npossVals, isChange, anyChange, debug, method;
  double fMin, newMin, mult, *allParam, stopTol, minv, maxv, lower, upper, h, prevMin, sameTol;
  char *ptype;
  

  allParam     = obj->allParam;
  nparam       = obj->nallParam;
  mult         = obj->initStepSize;
  ok           = 0;
  npossVals    = obj->npossVals;
  stopTol      = obj->stopTol;
  debug        = obj->debug;  
  ptype        = obj->paramType;
  method       = obj->method_rr;
  fMin         = obj->fmin;
  prevMin      = DBL_MISS;
  sameTol      = 1000.0*MACHINE_EPS;
  if (obj->fix_r) nparam = nparam - 1;

/*
int prtFlag;
if (obj->fix_r) {
prtFlag = 1;
} else {
prtFlag = 0;
}
if (prtFlag) {
  print_dVec(allParam, obj->nallParam, "parms");
  Rprintf("%s\n", ptype);
}
*/

  for (nattempt=1; nattempt<=obj->maxiter; nattempt++) {    
    anyChange = 0;
    for (i=0; i<nparam; i++) {
      curIndex = i;

      /* Get the bounds and step size h for this parm */
      optimRange(obj, curIndex, &minv, &maxv);


      h = getStepAndBounds(method, allParam[curIndex], ptype[curIndex], mult, npossVals, 
                           minv, maxv, &lower, &upper);
/*
if (prtFlag) {
Rprintf("natt=%d, i=%d, minv=%g, maxv=%g, lower=%g, upper=%g, mult=%g\n", 
         nattempt, i+1, minv, maxv, lower, upper, mult);
}
*/

      /* Minimize for this parm */
      fMin = findMin2(obj, allParam, curIndex, h, npossVals, lower, upper, minv, maxv, fMin, &isChange);
      anyChange += isChange;

/*
if (prtFlag) {
Rprintf("minx=%g, fmin=%10.6f\n", allParam[curIndex], fMin);
}
*/
    }  
    newMin = fL(allParam, obj);
    if (debug > 1) {
      Rprintf("iter=%d, newMin=%g, mult=%g, anyChange=%d\n", nattempt, newMin, mult, anyChange);
      print_dVec(allParam, obj->nallParam, "parms");
    }

    if (fabs(newMin - prevMin) < sameTol) {
      mult *= 0.9;
    }
 
/*
if (obj->fix_r ) {
Rprintf("%g, %g\n", prevMin, newMin);
}
*/


    prevMin = newMin;

    if (mult < stopTol) {
      ok = 1;  
      break;
    }
    fMin = newMin;
  }
  obj->fmin = newMin;  

  return(ok);

} /* END: optim_force_main */

static int optim_force(obj)
OBJ *obj;
{
  int debug=obj->debug, ok;
  double f0;

  /* Copy initial estimates */
  copy_dVec(obj->allParam0, obj->nallParam, obj->allParam);
  if (debug) print_dVec(obj->allParam, obj->nallParam, "parms");

  /* Initial function value */
  f0 = fL(obj->allParam, obj);
  if (debug) Rprintf("fMin = %g\n", f0);
  obj->fmin = f0;

  ok = optim_force_main(obj);

  return(ok);

} /* END: optim_force */

void C_optim_force(iargs, dargs, x1, x2, x3, x4, n1x1, n2x2, n3x3, n4x4, Y, tr, Z, 
                   allParam0, minParmValue, maxParmValue, parmOrder, parmType,
                   retParms, retOK, retMin) 
int *iargs, *Y, *tr, *retOK, *parmOrder;
double *dargs, *x1, *x2, *x3, *x4, *n1x1, *n2x2, *n3x3, *n4x4, *Z, *allParam0, 
       *minParmValue, *maxParmValue, *retParms, *retMin;
char **parmType;
{
  OBJ obj;

  obj_init(&obj);
  obj_setup(&obj, iargs, dargs, x1, x2, x3, x4, n1x1, n2x2, n3x3, n4x4, Y, tr, Z, 
            allParam0, retParms, minParmValue, maxParmValue, parmOrder, *parmType);
  *retOK  = optim_force(&obj);
  *retMin = obj.fmin;

  obj_free(&obj);

  return;

} /* END: C_optim_force */
