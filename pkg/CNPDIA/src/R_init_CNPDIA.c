#ifndef R_R_H
# include <R.h>
#endif

#ifndef R_EXT_DYNLOAD_H_
# include <R_ext/Rdynload.h>
#endif


#include <Rinternals.h>
#include <stdlib.h> // for NULL

/* register native routines ------------------------------------------------ */


/* NO .C calls */
/* NO .Call calls */
 
/* .Fortran calls */
void F77_NAME(initomexdiap)(void (* steadyparms)(int *, double *));
void F77_NAME(initforcp)   (void (* steadyforcs)(int *, double *));
void F77_NAME(initmpbdiap) (void (* steadyparms)(int *, double *));

void F77_NAME(omexdiamodp) (int *, double *, double *, double *, double *, int *);
void F77_NAME(omexdiamodbw)(int *, double *, double *, double *, double *, int *);
void F77_NAME(mpbdiamodp)  (int *, double *, double *, double *, double *, int *);
void F77_NAME(mpbdiamodbw) (int *, double *, double *, double *, double *, int *);
 
R_FortranMethodDef FEntries[] = {
    {"initomexdiap",    (DL_FUNC) &F77_SUB(initomexdiap),   1},
    {"initforcp",       (DL_FUNC) &F77_SUB(initforcp),      1},
    {"initmpbdiap",     (DL_FUNC) &F77_SUB(initmpbdiap),    1},
    {"omexdiamodp",     (DL_FUNC) &F77_SUB(omexdiamodp),    6},
    {"omexdiamodbw",    (DL_FUNC) &F77_SUB(omexdiamodbw),   6},
    {"mpbdiamodp",      (DL_FUNC) &F77_SUB(mpbdiamodp),     6},
    {"mpbdiamodbw",     (DL_FUNC) &F77_SUB(mpbdiamodbw),    6},
    {NULL, NULL, 0}
};

/* Initialization ---------------------------------------------------------- */
void R_init_CNPDIA(DllInfo *dll) {

  R_registerRoutines(dll, NULL, NULL, FEntries, NULL);

  // the following line protects against accidentially finding entry points

  R_useDynamicSymbols(dll, FALSE); // disable dynamic searching
}
