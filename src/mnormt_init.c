#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(sadmvn)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(sadmvt)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(smvbvt)(void *, void *, void *, void *, void *, void *); 
 
extern void F77_NAME(stvtl)(void *, void *, void *, void *,  void *); 
 
static const R_FortranMethodDef FortEntries[] = {
    {"sadmvn",          (DL_FUNC) &F77_NAME(sadmvn), 11},
    {"sadmvt",          (DL_FUNC) &F77_NAME(sadmvt), 12},
    {"smvbvt",          (DL_FUNC) &F77_NAME(sadmvn), 6},
    {"stvtl",           (DL_FUNC) &F77_NAME(stvtl), 5},
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
