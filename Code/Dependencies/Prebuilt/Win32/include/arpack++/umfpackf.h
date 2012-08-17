/*
  ARPACK++ v1.2 2/20/2000
  c++ interface to ARPACK code.

  MODULE umfpackf.h
  UMFPACK FORTRAN routines.

  ARPACK Authors
     Richard Lehoucq
     Danny Sorensen
     Chao Yang
     Dept. of Computational & Applied Mathematics
     Rice University
     Houston, Texas
*/

#ifndef UMFPACKF_H
#define UMFPACKF_H

#include "arch.h"

extern "C"
{

  // Single precision real routines.

  void F77NAME(ums21i)(ARint *keep, float *cntl, ARint *icntl);

  void F77NAME(ums2fa)(const ARint *n, const ARint *ne, 
                       const ARint *job, const ARlogical *transa,
                       const ARint *lvalue, const ARint *lindex,
                       float *value, ARint *index, ARint *keep,
                       const float *cntl, const ARint *icntl,
                       ARint *info, float *rinfo); 

  void F77NAME(ums2so)(const ARint *n, const ARint *job,
                       const ARlogical *transc, const ARint *lvalue,
                       const ARint *lindex, float *value,
                       ARint *index, const ARint *keep, 
                       const float *b, float *x, float *w,
                       const float *cntl, const ARint *icntl,
                       ARint *info, float *rinfo);


  // Double precision real routines.

  void F77NAME(umd21i)(ARint *keep, double *cntl, ARint *icntl);

  void F77NAME(umd2fa)(const ARint *n, const ARint *ne, 
                       const ARint *job, const ARlogical *transa,
                       const ARint *lvalue, const ARint *lindex,
                       double *value, ARint *index, ARint *keep,
                       const double *cntl, const ARint *icntl,
                       ARint *info, double *rinfo); 

  void F77NAME(umd2so)(const ARint *n, const ARint *job,
                       const ARlogical *transc, const ARint *lvalue,
                       const ARint *lindex, double *value,
                       ARint *index, const ARint *keep, 
                       const double *b, double *x, double *w,
                       const double *cntl, const ARint *icntl,
                       ARint *info, double *rinfo);


  // Single precision complex routines.

#ifdef ARCOMP_H

  void F77NAME(umc21i)(ARint *keep, float *cntl, ARint *icntl);

  void F77NAME(umc2fa)(const ARint *n, const ARint *ne, 
                       const ARint *job, const ARlogical *transa,
                       const ARint *lvalue, const ARint *lindex,
                       arcomplex<float> *value, ARint *index, ARint *keep,
                       const float *cntl, const ARint *icntl,
                       ARint *info, float *rinfo); 

  void F77NAME(umc2so)(const ARint *n, const ARint *job,
                       const ARlogical *transc, const ARint *lvalue,
                       const ARint *lindex, arcomplex<float> *value,
                       ARint *index, const ARint *keep, 
                       const arcomplex<float> *b, arcomplex<float> *x, 
                       arcomplex<float> *w, const float *cntl, 
                       const ARint *icntl, ARint *info, float *rinfo);


  // Double precision complex routines.

  void F77NAME(umz21i)(ARint *keep, double *cntl, ARint *icntl);

  void F77NAME(umz2fa)(const ARint *n, const ARint *ne, 
                       const ARint *job, const ARlogical *transa,
                       const ARint *lvalue, const ARint *lindex,
                       arcomplex<double> *value, ARint *index, ARint *keep,
                       const double *cntl, const ARint *icntl,
                       ARint *info, double *rinfo); 

  void F77NAME(umz2so)(const ARint *n, const ARint *job,
                       const ARlogical *transc, const ARint *lvalue,
                       const ARint *lindex, arcomplex<double> *value,
                       ARint *index, const ARint *keep, 
                       const arcomplex<double> *b, arcomplex<double> *x, 
                       arcomplex<double> *w, const double *cntl, 
                       const ARint *icntl, ARint *info, double *rinfo);

#endif // ARCOMP_H

}
#endif // UMFPACKF_H

