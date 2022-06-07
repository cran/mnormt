#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }
double F77_SUB(unifrnd)(void) { return unif_rand(); }

double F77_SUB(pnormr)(double *x, double *mu, double *sigma, int *lower_tail, int *log_p) { return pnorm(*x, *mu, *sigma, *lower_tail, *log_p); }
double F77_SUB(qnormr)(double *p, double *mu, double *sigma, int *lower_tail, int *log_p) { return qnorm(*p, *mu, *sigma, *lower_tail, *log_p); }

/* .Fortran calls */
extern void F77_NAME(sadmvn)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(sadmvt)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(smvbvt)(void *, void *, void *, void *, void *, void *); 
 
extern void F77_NAME(stvtl)(void *, void *, void *, void *,  void *); 

extern void F77_NAME(rtmng)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
 
static const R_FortranMethodDef FortEntries[] = {
    {"sadmvn",          (DL_FUNC) &F77_NAME(sadmvn), 11},
    {"sadmvt",          (DL_FUNC) &F77_NAME(sadmvt), 12},
    {"smvbvt",          (DL_FUNC) &F77_NAME(sadmvn), 6},
    {"stvtl",           (DL_FUNC) &F77_NAME(stvtl), 5},
    {"rtmng",           (DL_FUNC) &F77_NAME(rtmng), 9},
    {NULL, NULL, 0}
};


static const R_CMethodDef CEntries[] = {
    {NULL, NULL, 0}
};

void R_init_heavy(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
