/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARUSNSym.h.
   Arpack++ class ARumNonSymStdEig definition
   (umfpack version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARUSNSYM_H
#define ARUSNSYM_H

#include <cstddef>

#include "arch.h"
#include "arsnsym.h"
#include "arunsmat.h"


template<class ARFLOAT>
class ARumNonSymStdEig:
  public virtual ARNonSymStdEig<ARFLOAT, ARumNonSymMatrix<ARFLOAT, ARFLOAT> > {

 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

 // a.2) Constructors and destructor.

  ARumNonSymStdEig() { }
  // Short constructor.

  ARumNonSymStdEig(int nevp, ARumNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   char* whichp = "LM", int ncvp = 0,
                   ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARumNonSymStdEig(int nevp, ARumNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   ARFLOAT sigma, char* whichp = "LM", int ncvp = 0,
                   ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARumNonSymStdEig(const ARumNonSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumNonSymStdEig() { }
  // Destructor.

 // b) Operators.

  ARumNonSymStdEig& operator=(const ARumNonSymStdEig& other);
  // Assignment operator.

}; // class ARumNonSymStdEig.


// ------------------------------------------------------------------------ //
// ARumNonSymStdEig member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARumNonSymStdEig<ARFLOAT>::ChangeShift(ARFLOAT sigmaRp)
{

  this->sigmaR    = sigmaRp;
  this->sigmaI    = 0.0;
  this->mode      = 3;
  this->iparam[7] = this->mode;

  this->objOP->FactorAsI(this->sigmaR);
  this->Restart();

} // ChangeShift.


template<class ARFLOAT>
inline void ARumNonSymStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARumNonSymMatrix<ARFLOAT, ARFLOAT> >::
    SetRegularMode(this->objOP, &ARumNonSymMatrix<ARFLOAT, ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARumNonSymStdEig<ARFLOAT>::SetShiftInvertMode(ARFLOAT sigmap)
{

  ARStdEig<ARFLOAT, ARFLOAT, ARumNonSymMatrix<ARFLOAT, ARFLOAT> >::
    SetShiftInvertMode(sigmap, this->objOP, 
                       &ARumNonSymMatrix<ARFLOAT, ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARumNonSymStdEig<ARFLOAT>::
ARumNonSymStdEig(int nevp, ARumNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 char* whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->NoShift();
  DefineParameters(A.ncols(), nevp, &A, 
                   &ARumNonSymMatrix<ARFLOAT, ARFLOAT>::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARumNonSymStdEig<ARFLOAT>::
ARumNonSymStdEig(int nevp, ARumNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 ARFLOAT sigmap, char* whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)

{

  DefineParameters(A.ncols(), nevp, &A, 
                   &ARumNonSymMatrix<ARFLOAT, ARFLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARumNonSymStdEig<ARFLOAT>& ARumNonSymStdEig<ARFLOAT>::
operator=(const ARumNonSymStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUSNSYM_H
