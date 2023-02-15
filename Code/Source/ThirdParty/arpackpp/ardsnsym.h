/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARDSNSym.h.
   Arpack++ class ARdsNonSymStdEig definition
   (dense matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARDSNSYM_H
#define ARDSNSYM_H

#include <cstddef>

#include "arch.h"
#include "arsnsym.h"
#include "ardnsmat.h"


template<class ARFLOAT>
class ARdsNonSymStdEig:
  public virtual ARNonSymStdEig<ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> > {

 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

 // a.2) Constructors and destructor.

  ARdsNonSymStdEig() { }
  // Short constructor.

  ARdsNonSymStdEig(int nevp, ARdsNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   const char* whichp = "LM", int ncvp = 0,
                   ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARdsNonSymStdEig(int nevp, ARdsNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   ARFLOAT sigma, const char* whichp = "LM", int ncvp = 0,
                   ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARdsNonSymStdEig(const ARdsNonSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARdsNonSymStdEig() { }
  // Destructor.

 // b) Operators.

  ARdsNonSymStdEig& operator=(const ARdsNonSymStdEig& other);
  // Assignment operator.

}; // class ARdsNonSymStdEig.


// ------------------------------------------------------------------------ //
// ARdsNonSymStdEig member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARdsNonSymStdEig<ARFLOAT>::
ChangeShift(ARFLOAT sigmaRp)
{

  this->sigmaR    = sigmaRp;
  this->sigmaI    = 0.0;
  this->mode      = 3;
  this->iparam[7] = this->mode;

  this->objOP->FactorAsI(this->sigmaR);
  this->Restart();

} // ChangeShift.


template<class ARFLOAT>
inline void ARdsNonSymStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> >::
    SetRegularMode(this->objOP, &ARdsNonSymMatrix<ARFLOAT, ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARdsNonSymStdEig<ARFLOAT>::SetShiftInvertMode(ARFLOAT sigmap)
{

  ARStdEig<ARFLOAT, ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> >::
    SetShiftInvertMode(sigmap, this->objOP,
                       &ARdsNonSymMatrix<ARFLOAT, ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARdsNonSymStdEig<ARFLOAT>::
ARdsNonSymStdEig(int nevp, ARdsNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 const char* whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &A,
                         &ARdsNonSymMatrix<ARFLOAT, ARFLOAT>::MultMv,
                         whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARdsNonSymStdEig<ARFLOAT>::
ARdsNonSymStdEig(int nevp, ARdsNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 ARFLOAT sigmap, const char* whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->DefineParameters(A.ncols(), nevp, &A,
                         &ARdsNonSymMatrix<ARFLOAT, ARFLOAT>::MultInvv,
                         whichp, ncvp, tolp, maxitp, residp, ishiftp);
  this->ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARdsNonSymStdEig<ARFLOAT>& ARdsNonSymStdEig<ARFLOAT>::
operator=(const ARdsNonSymStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDSNSYM_H
