/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARBSComp.h.
   Arpack++ class ARbdCompStdEig definition
   (band matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARBSCOMP_H
#define ARBSCOMP_H

#include <cstddef>
#include "arch.h"
#include "arscomp.h"
#include "arbnsmat.h"
#include "arrseig.h"


template<class ARFLOAT>
class ARbdCompStdEig:
  public virtual ARCompStdEig<ARFLOAT,
                              ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> > {

 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(arcomplex<ARFLOAT> sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(arcomplex<ARFLOAT> sigmap);

 // a.2) Constructors and destructor.

  ARbdCompStdEig() { }
  // Short constructor.

  ARbdCompStdEig(int nevp, ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 const char* whichp = "LM", int ncvp = 0,
                 ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARbdCompStdEig(int nevp, ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 arcomplex<ARFLOAT> sigma, const char* whichp = "LM",
                 int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARbdCompStdEig(const ARbdCompStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARbdCompStdEig() { }
  // Destructor.


 // b) Operators.

  ARbdCompStdEig& operator=(const ARbdCompStdEig& other);
  // Assignment operator.

}; // class ARbdCompStdEig.


// ------------------------------------------------------------------------ //
// ARbdCompStdEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARbdCompStdEig<ARFLOAT>::
ChangeShift(arcomplex<ARFLOAT> sigmaRp)
{

  this->objOP->FactorAsI(sigmaRp);
  ARrcStdEig<ARFLOAT, arcomplex<ARFLOAT> >::ChangeShift(sigmaRp);

} // ChangeShift.


template<class ARFLOAT>
inline void ARbdCompStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>,
           ARbdNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> >::
    SetRegularMode(this->objOP,&ARbdNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARbdCompStdEig<ARFLOAT>::
SetShiftInvertMode(arcomplex<ARFLOAT> sigmap)
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>,
           ARbdNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT> >::
    SetShiftInvertMode(sigmap, this->objOP,
                       &ARbdNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARbdCompStdEig<ARFLOAT>::
ARbdCompStdEig(int nevp, ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               const char* whichp, int ncvp, ARFLOAT tolp,
               int maxitp, arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &A,
                   &ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARbdCompStdEig<ARFLOAT>::
ARbdCompStdEig(int nevp, ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               arcomplex<ARFLOAT> sigmap, const char* whichp, int ncvp,
               ARFLOAT tolp, int maxitp, arcomplex<ARFLOAT>* residp,
               bool ishiftp)

{

  this->DefineParameters(A.ncols(), nevp, &A,
                   &ARbdNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  this->ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARbdCompStdEig<ARFLOAT>& ARbdCompStdEig<ARFLOAT>::
operator=(const ARbdCompStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBSCOMP_H
