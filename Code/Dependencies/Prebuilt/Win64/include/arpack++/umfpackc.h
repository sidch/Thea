/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE umfpackc.h.
   Interface to UMFPACK FORTRAN routines.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arch.h"
#include "umfpackf.h"

#ifndef UMFPACKC_H
#define UMFPACKC_H


void AfterUm21i(int keep[], int icntl[], bool simest, bool reducible) 

// A function that changes some components of icntl and keep vectors. 

{

  icntl[2]  = 0;              // icntl(3)=0 to avoid printing umfpack messages.
  icntl[3]  = (int)reducible; // icntl(4)=1 if matrix must be permuted to
                              // block triangular form.
  icntl[5]  = (int)simest;    // icntl(6)=1 if matrix is structurally symmetric

  // Parameters defined in arch.h file.

  if (UICNTL5 != 0) icntl[4] = UICNTL5; // # cols examined during pivot search.
  if (UICNTL7 != 0) icntl[6] = UICNTL7; // block size for the blas.
  if (UKEEP7  != 0)  keep[6] = UKEEP7;  // dense columns absolute threshold.
  if (UKEEP8  != 0)  keep[7] = UKEEP8;  // dense columns relative threshold.

} // AfterUm21i.


// UM21I

inline void um21i(ARint keep[], float cntl[], ARint icntl[],
                  double threshold = 0.1, bool simest = false, 
                  bool reducible = true) 

// An extended version of UMFPACK's ums21i.

{

  // Calling the FORTRAN function.

  F77NAME(ums21i)(keep, cntl, icntl);

  // Changing cntl vector.

  cntl[0] = (float)threshold;                // relative pivot tolerance.
  if (UCNTL2 != 0) cntl[1] = (float)UCNTL2;  // amalgamation parameter.

  // Changing keep and icntl.

  AfterUm21i(keep, icntl, simest, reducible);

} // um21i (float)

inline void um21i(ARint keep[], double cntl[], ARint icntl[],
                  double threshold = 0.1, bool simest = false, 
                  bool reducible = true) 

// An extended version of UMFPACK's umd21i.

{

  // Calling the FORTRAN function.

  F77NAME(umd21i)(keep, cntl, icntl);

  // Changing cntl vector.

  cntl[0] = threshold;                       // relative pivot tolerance.
  if (UCNTL2 != 0) cntl[1] = (double)UCNTL2; // amalgamation parameter.

  // Changing keep and icntl.

  AfterUm21i(keep, icntl, simest, reducible);

} // um21i (double)

#ifdef ARCOMP_H

inline void um21i(ARint keep[], arcomplex<float> cntl[], ARint icntl[],
                  double threshold = 0.1, bool simest = false, 
                  bool reducible = true) 

// Another extended version of UMFPACK's ums21i (for complex problems).

{

  // Calling the FORTRAN function.

  float* fcntl = (float*)cntl;
  F77NAME(ums21i)(keep, fcntl, icntl);

  // Changing cntl vector.

  fcntl[0] = (float)threshold;               // relative pivot tolerance.
  if (UCNTL2 != 0) fcntl[1] = (float)UCNTL2; // amalgamation parameter.

  // Changing keep and icntl.

  AfterUm21i(keep, icntl, simest, reducible);

} // um21i (arcomplex<float>)

inline void um21i(ARint keep[], arcomplex<double> cntl[], ARint icntl[],
                  double threshold = 0.1, bool simest = false, 
                  bool reducible = true) 

// Another extended version of UMFPACK's umd21i (for complex problems).

{

  // Calling the FORTRAN function.

  double* fcntl = (double*)cntl;
  F77NAME(umd21i)(keep, fcntl, icntl);

  // Changing cntl vector.

  fcntl[0] = threshold;                       // relative pivot tolerance.
  if (UCNTL2 != 0) fcntl[1] = (double)UCNTL2; // amalgamation parameter.

  // Changing keep and icntl.

  AfterUm21i(keep, icntl, simest, reducible);

} // um21i (arcomplex<double>)

#endif


// UM2FA

inline void um2fa(const ARint &n, const ARint &ne, const ARint &job,
                  const ARlogical &transa, const ARint &lvalue,
                  const ARint &lindex, float value[], ARint index[],
                  ARint keep[], const float cntl[], 
                  const ARint icntl[], ARint info[], float rinfo[]) {
  F77NAME(ums2fa)(&n, &ne, &job, &transa, &lvalue, &lindex, 
                  value, index, keep, cntl, icntl, info, rinfo);
} // um2fa (float)

inline void um2fa(const ARint &n, const ARint &ne, const ARint &job,
                  const ARlogical &transa, const ARint &lvalue,
                  const ARint &lindex, double value[], ARint index[],
                  ARint keep[], const double cntl[], 
                  const ARint icntl[], ARint info[], double rinfo[]) {
  F77NAME(umd2fa)(&n, &ne, &job, &transa, &lvalue, &lindex, 
                  value, index, keep, cntl, icntl, info, rinfo);
} // um2fa (double)

#ifdef ARCOMP_H

inline void um2fa(const ARint &n, const ARint &ne, const ARint &job,
                  const ARlogical &transa, const ARint &lvalue,
                  const ARint &lindex, arcomplex<float> value[], 
                  ARint index[], ARint keep[], 
                  const arcomplex<float> cntl[], const ARint icntl[], 
                  ARint info[], arcomplex<float> rinfo[]) {
  F77NAME(umc2fa)(&n, &ne, &job, &transa, &lvalue, &lindex, value, index, 
                  keep, (float*)cntl, icntl, info, (float*)rinfo);
} // um2fa (arcomplex<float>)

inline void um2fa(const ARint &n, const ARint &ne, const ARint &job,
                  const ARlogical &transa, const ARint &lvalue,
                  const ARint &lindex, arcomplex<double> value[], 
                  ARint index[], ARint keep[], 
                  const arcomplex<double> cntl[], const ARint icntl[], 
                  ARint info[], arcomplex<double> rinfo[]) {
  F77NAME(umz2fa)(&n, &ne, &job, &transa, &lvalue, &lindex, value, index, 
                  keep, (double*)cntl, icntl, info, (double*)rinfo);
} // um2fa (arcomplex<double>)

#endif


// UM2SO

inline void um2so(const ARint &n, const ARint &job, 
                  const ARlogical &transc, const ARint &lvalue, 
                  const ARint &lindex, float value[], ARint index[],
                  const ARint keep[], const float b[], float x[], 
                  float w[], const float cntl[], const ARint icntl[], 
                  ARint info[], float rinfo[]) {
  F77NAME(ums2so)(&n, &job, &transc, &lvalue, &lindex, value,
                  index, keep, b, x, w, cntl, icntl, info, rinfo);
} // um2so (float)

inline void um2so(const ARint &n, const ARint &job, 
                  const ARlogical &transc, const ARint &lvalue, 
                  const ARint &lindex, double value[], ARint index[],
                  const ARint keep[], const double b[], double x[], 
                  double w[], const double cntl[], const ARint icntl[], 
                  ARint info[], double rinfo[]) {
  F77NAME(umd2so)(&n, &job, &transc, &lvalue, &lindex, value,
                  index, keep, b, x, w, cntl, icntl, info, rinfo);
} // um2so (double)

#ifdef ARCOMP_H

inline void um2so(const ARint &n, const ARint &job, 
                  const ARlogical &transc, const ARint &lvalue, 
                  const ARint &lindex, arcomplex<float> value[], 
                  ARint index[], const ARint keep[], const 
                  arcomplex<float> b[], arcomplex<float> x[], 
                  arcomplex<float> w[], const arcomplex<float> cntl[], 
                  const ARint icntl[], ARint info[], 
                  arcomplex<float> rinfo[]) {
  F77NAME(umc2so)(&n, &job, &transc, &lvalue, &lindex, value, index, keep,
                  b, x, w, (float*)cntl, icntl, info, (float*)rinfo);
} // um2so (arcomplex<float>)

inline void um2so(const ARint &n, const ARint &job, 
                  const ARlogical &transc, const ARint &lvalue, 
                  const ARint &lindex, arcomplex<double> value[], 
                  ARint index[], const ARint keep[], 
                  const arcomplex<double> b[], arcomplex<double> x[], 
                  arcomplex<double> w[], const arcomplex<double> cntl[], 
                  const ARint icntl[], ARint info[], 
                  arcomplex<double> rinfo[]) {
  F77NAME(umz2so)(&n, &job, &transc, &lvalue, &lindex, value, index, keep,
                  b, x, w, (double*)cntl, icntl, info, (double*)rinfo);
} // um2so (arcomplex<double>)

#endif


#endif // UMFPACKC_H
