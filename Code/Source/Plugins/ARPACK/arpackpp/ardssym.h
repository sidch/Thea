/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARDSSym.h.
   Arpack++ class ARdsSymStdEig definition
   (dense matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARDSSYM_H
#define ARDSSYM_H

#include <cstddef>

#include "arch.h"
#include "arssym.h"
#include "ardsmat.h"


template<class ARFLOAT>
class ARdsSymStdEig:
  public virtual ARSymStdEig<ARFLOAT, ARdsSymMatrix<ARFLOAT> > {

 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

 // a.2) Constructors and destructor.

  ARdsSymStdEig() { }
  // Short constructor.

  ARdsSymStdEig(int nevp, ARdsSymMatrix<ARFLOAT>& A,
                const char* whichp = "LM", int ncvp = 0,
                ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARdsSymStdEig(int nevp, ARdsSymMatrix<ARFLOAT>& A,
                ARFLOAT sigma, const char* whichp = "LM", int ncvp = 0,
                ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARdsSymStdEig(const ARdsSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARdsSymStdEig() { }
  // Destructor.

 // b) Operators.

  ARdsSymStdEig& operator=(const ARdsSymStdEig& other);
  // Assignment operator.

}; // class ARdsSymStdEig.


// ------------------------------------------------------------------------ //
// ARdsSymStdEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARdsSymStdEig<ARFLOAT>::
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
inline void ARdsSymStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARdsSymMatrix<ARFLOAT> >::
    SetRegularMode(this->objOP, &ARdsSymMatrix<ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARdsSymStdEig<ARFLOAT>::SetShiftInvertMode(ARFLOAT sigmap)
{

  ARStdEig<ARFLOAT, ARFLOAT, ARdsSymMatrix<ARFLOAT> >::
    SetShiftInvertMode(sigmap, this->objOP, &ARdsSymMatrix<ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARdsSymStdEig<ARFLOAT>::
ARdsSymStdEig(int nevp, ARdsSymMatrix<ARFLOAT>& A,
              const char* whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)
{

  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &A, &ARdsSymMatrix<ARFLOAT>::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARdsSymStdEig<ARFLOAT>::
ARdsSymStdEig(int nevp, ARdsSymMatrix<ARFLOAT>& A,
              ARFLOAT sigmap, const char* whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->DefineParameters(A.ncols(), nevp, &A, &ARdsSymMatrix<ARFLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  this->ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARdsSymStdEig<ARFLOAT>& ARdsSymStdEig<ARFLOAT>::
operator=(const ARdsSymStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDSSYM_H
