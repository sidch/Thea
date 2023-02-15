/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARUSSym.h.
   Arpack++ class ARumSymStdEig definition
   (UMFPACK version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARUSSYM_H
#define ARUSSYM_H

#include <cstddef>

#include "arch.h"
#include "arssym.h"
#include "arusmat.h"


template<class ARFLOAT>
class ARumSymStdEig:
  public virtual ARSymStdEig<ARFLOAT, ARumSymMatrix<ARFLOAT> > {

 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

 // a.2) Constructors and destructor.

  ARumSymStdEig() { }
  // Short constructor.

  ARumSymStdEig(int nevp, ARumSymMatrix<ARFLOAT>& A,
                const char* whichp = "LM", int ncvp = 0,
                ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARumSymStdEig(int nevp, ARumSymMatrix<ARFLOAT>& A,
                ARFLOAT sigma, const char* whichp = "LM", int ncvp = 0,
                ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARumSymStdEig(const ARumSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumSymStdEig() { }
  // Destructor.

 // b) Operators.

  ARumSymStdEig& operator=(const ARumSymStdEig& other);
  // Assignment operator.

}; // class ARumSymStdEig.


// ------------------------------------------------------------------------ //
// ARumSymStdEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARumSymStdEig<ARFLOAT>::ChangeShift(ARFLOAT sigmaRp)
{

  this->sigmaR    = sigmaRp;
  this->sigmaI    = 0.0;
  this->mode      = 3;
  this->iparam[7] = this->mode;

  this->objOP->FactorAsI(this->sigmaR);
  this->Restart();

} // ChangeShift.


template<class ARFLOAT>
inline void ARumSymStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARumSymMatrix<ARFLOAT> >::
    SetRegularMode(this->objOP, &ARumSymMatrix<ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARumSymStdEig<ARFLOAT>::SetShiftInvertMode(ARFLOAT sigmap)
{

  ARStdEig<ARFLOAT, ARFLOAT, ARumSymMatrix<ARFLOAT> >::
    SetShiftInvertMode(sigmap, this->objOP, &ARumSymMatrix<ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARumSymStdEig<ARFLOAT>::
ARumSymStdEig(int nevp, ARumSymMatrix<ARFLOAT>& A,
              const char* whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)
{

  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &A, &ARumSymMatrix<ARFLOAT>::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARumSymStdEig<ARFLOAT>::
ARumSymStdEig(int nevp, ARumSymMatrix<ARFLOAT>& A,
              ARFLOAT sigmap, const char* whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->DefineParameters(A.ncols(), nevp, &A, &ARumSymMatrix<ARFLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  this->ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARumSymStdEig<ARFLOAT>& ARumSymStdEig<ARFLOAT>::
operator=(const ARumSymStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUSSYM_H
