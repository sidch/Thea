/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARBSSym.h.
   Arpack++ class ARbdSymStdEig definition
   (band matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARBSSYM_H
#define ARBSSYM_H

#include <cstddef>
#include "arch.h"
#include "arssym.h"
#include "arbsmat.h"


template<class ARFLOAT>
class ARbdSymStdEig:
  public virtual ARSymStdEig<ARFLOAT, ARbdSymMatrix<ARFLOAT> > {

 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

 // a.2) Constructors and destructor.

  ARbdSymStdEig() { }
  // Short constructor.

  ARbdSymStdEig(int nevp, ARbdSymMatrix<ARFLOAT>& A,
                char* whichp = "LM", int ncvp = 0,
                ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARbdSymStdEig(int nevp, ARbdSymMatrix<ARFLOAT>& A,
                ARFLOAT sigma, char* whichp = "LM", int ncvp = 0,
                ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARbdSymStdEig(const ARbdSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARbdSymStdEig() { }
  // Destructor.

 // b) Operators.

  ARbdSymStdEig& operator=(const ARbdSymStdEig& other);
  // Assignment operator.

}; // class ARbdSymStdEig.


// ------------------------------------------------------------------------ //
// ARbdSymStdEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARbdSymStdEig<ARFLOAT>::
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
inline void ARbdSymStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARbdSymMatrix<ARFLOAT> >::
    SetRegularMode(this->objOP, &ARbdSymMatrix<ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARbdSymStdEig<ARFLOAT>::SetShiftInvertMode(ARFLOAT sigmap)
{

  ARStdEig<ARFLOAT, ARFLOAT, ARbdSymMatrix<ARFLOAT> >::
    SetShiftInvertMode(sigmap, this->objOP, &ARbdSymMatrix<ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARbdSymStdEig<ARFLOAT>::
ARbdSymStdEig(int nevp, ARbdSymMatrix<ARFLOAT>& A,
              char* whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)
{

  this->NoShift();
  DefineParameters(A.ncols(), nevp, &A, &ARbdSymMatrix<ARFLOAT>::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARbdSymStdEig<ARFLOAT>::
ARbdSymStdEig(int nevp, ARbdSymMatrix<ARFLOAT>& A,
              ARFLOAT sigmap, char* whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)

{

  DefineParameters(A.ncols(), nevp, &A, &ARbdSymMatrix<ARFLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARbdSymStdEig<ARFLOAT>& ARbdSymStdEig<ARFLOAT>::
operator=(const ARbdSymStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBSSYM_H
