/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARUSComp.h.
   Arpack++ class ARumCompStdEig definition
   (umfpack version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARUSCOMP_H
#define ARUSCOMP_H

#include <cstddef>

#include "arch.h"
#include "arscomp.h"
#include "arunsmat.h"
#include "arrseig.h"


template<class ARFLOAT>
class ARumCompStdEig:
  public virtual ARCompStdEig<ARFLOAT,
                              ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> > {

 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(arcomplex<ARFLOAT> sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(arcomplex<ARFLOAT> sigmap);

 // a.2) Constructors and destructor.

  ARumCompStdEig() { }
  // Short constructor.

  ARumCompStdEig(int nevp, ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 const char* whichp = "LM", int ncvp = 0,
                 ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARumCompStdEig(int nevp, ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 arcomplex<ARFLOAT> sigma, const char* whichp = "LM",
                 int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARumCompStdEig(const ARumCompStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumCompStdEig() { }
  // Destructor.


 // b) Operators.

  ARumCompStdEig& operator=(const ARumCompStdEig& other);
  // Assignment operator.

}; // class ARumCompStdEig.


// ------------------------------------------------------------------------ //
// ARumCompStdEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARumCompStdEig<ARFLOAT>::
ChangeShift(arcomplex<ARFLOAT> sigmaRp)
{

  this->objOP->FactorAsI(sigmaRp);
  ARrcStdEig<ARFLOAT, arcomplex<ARFLOAT> >::ChangeShift(sigmaRp);

} // ChangeShift.


template<class ARFLOAT>
inline void ARumCompStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>,
           ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetRegularMode(this->objOP,
                   &ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARumCompStdEig<ARFLOAT>::
SetShiftInvertMode(arcomplex<ARFLOAT> sigmap)
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>,
           ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetShiftInvertMode(sigmap, this->objOP,
                       &ARumNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARumCompStdEig<ARFLOAT>::
ARumCompStdEig(int nevp, ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               const char* whichp, int ncvp, ARFLOAT tolp,
               int maxitp, arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &A,
                   &ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARumCompStdEig<ARFLOAT>::
ARumCompStdEig(int nevp, ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               arcomplex<ARFLOAT> sigmap, const char* whichp, int ncvp,
               ARFLOAT tolp, int maxitp, arcomplex<ARFLOAT>* residp,
               bool ishiftp)

{

  this->DefineParameters(A.ncols(), nevp, &A,
                   &ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  this->ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARumCompStdEig<ARFLOAT>& ARumCompStdEig<ARFLOAT>::
operator=(const ARumCompStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUSCOMP_H
