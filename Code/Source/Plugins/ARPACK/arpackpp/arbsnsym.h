/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARBSNSym.h.
   Arpack++ class ARbdNonSymStdEig definition
   (band matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARBSNSYM_H
#define ARBSNSYM_H

#include <cstddef>
#include "arch.h"
#include "arsnsym.h"
#include "arbnsmat.h"


template<class ARFLOAT>
class ARbdNonSymStdEig:
  public virtual ARNonSymStdEig<ARFLOAT, ARbdNonSymMatrix<ARFLOAT, ARFLOAT> > {

 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

 // a.2) Constructors and destructor.

  ARbdNonSymStdEig() { }
  // Short constructor.

  ARbdNonSymStdEig(int nevp, ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   const char* whichp = "LM", int ncvp = 0,
                   ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARbdNonSymStdEig(int nevp, ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   ARFLOAT sigma, const char* whichp = "LM", int ncvp = 0,
                   ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARbdNonSymStdEig(const ARbdNonSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARbdNonSymStdEig() { }
  // Destructor.

 // b) Operators.

  ARbdNonSymStdEig& operator=(const ARbdNonSymStdEig& other);
  // Assignment operator.

}; // class ARbdNonSymStdEig.


// ------------------------------------------------------------------------ //
// ARbdNonSymStdEig member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARbdNonSymStdEig<ARFLOAT>::
ChangeShift(ARFLOAT sigmaRp)
{

   this->sigmaR    = sigmaRp;
   this->sigmaI    = 0.0;
   this->mode      = 3;
   this->iparam[7] =  this->mode;

   this->objOP->FactorAsI( this->sigmaR);
   this->Restart();

} // ChangeShift.


template<class ARFLOAT>
inline void ARbdNonSymStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARbdNonSymMatrix<ARFLOAT, ARFLOAT> >::
    SetRegularMode( this->objOP, &ARbdNonSymMatrix<ARFLOAT, ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARbdNonSymStdEig<ARFLOAT>::SetShiftInvertMode(ARFLOAT sigmap)
{

  ARStdEig<ARFLOAT, ARFLOAT, ARbdNonSymMatrix<ARFLOAT, ARFLOAT> >::
    SetShiftInvertMode(sigmap,  this->objOP,
                       &ARbdNonSymMatrix<ARFLOAT, ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARbdNonSymStdEig<ARFLOAT>::
ARbdNonSymStdEig(int nevp, ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 const char* whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)

{

   this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &A,
                   &(ARbdNonSymMatrix<ARFLOAT, ARFLOAT>::MultMv),
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARbdNonSymStdEig<ARFLOAT>::
ARbdNonSymStdEig(int nevp, ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 ARFLOAT sigmap, const char* whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)

{

  this->DefineParameters(A.ncols(), nevp, &A,
                   &ARbdNonSymMatrix<ARFLOAT, ARFLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  this->ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARbdNonSymStdEig<ARFLOAT>& ARbdNonSymStdEig<ARFLOAT>::
operator=(const ARbdNonSymStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
     this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBSNSYM_H
