/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARUGComp.h.
   Arpack++ class ARumCompGenEig definition
   (umfpack version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARUGCOMP_H
#define ARUGCOMP_H

#include <cstddef>

#include "arch.h"
#include "arunsmat.h"
#include "arunspen.h"
#include "arrseig.h"
#include "argcomp.h"


template<class ARFLOAT>
class ARumCompGenEig:
  public virtual
    ARCompGenEig<ARFLOAT, ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT >,
                 ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT > > {

 private:

 // a) Data structure used to store matrices.

  ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT > Pencil;

 // b) Protected functions:

  virtual void Copy(const ARumCompGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(arcomplex<ARFLOAT> sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(arcomplex<ARFLOAT> sigmap);

 // c.2) Constructors and destructor.

  ARumCompGenEig() { }
  // Short constructor.

  ARumCompGenEig(int nevp, ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B,
                 char* whichp = "LM", int ncvp = 0,
                 ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARumCompGenEig(int nevp, ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B,
                 arcomplex<ARFLOAT> sigma, char* whichp = "LM",
                 int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARumCompGenEig(const ARumCompGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumCompGenEig() { }

 // d) Operators.

  ARumCompGenEig& operator=(const ARumCompGenEig& other);
  // Assignment operator.

}; // class ARumCompGenEig.


// ------------------------------------------------------------------------ //
// ARumCompGenEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARumCompGenEig<ARFLOAT>::
Copy(const ARumCompGenEig<ARFLOAT>& other)
{

  ARCompGenEig<ARFLOAT, ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT >,
               ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  this->objOP  = &Pencil;
  this->objB   = &Pencil;

} // Copy.


template<class ARFLOAT>
inline void ARumCompGenEig<ARFLOAT>::
ChangeShift(arcomplex<ARFLOAT> sigmaRp)
{

  this->objOP->FactorAsB(sigmaRp);
  ARrcStdEig<ARFLOAT, arcomplex<ARFLOAT> >::ChangeShift(sigmaRp);

} // ChangeShift.


template<class ARFLOAT>
inline void ARumCompGenEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>,
           ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetRegularMode(&Pencil,
                   &ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultInvBAv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARumCompGenEig<ARFLOAT>::
SetShiftInvertMode(arcomplex<ARFLOAT> sigmap)
{

  ARCompGenEig<ARFLOAT, ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>,
               ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil,
                       &ARumNonSymPencil<arcomplex<ARFLOAT>,ARFLOAT>::MultInvAsBv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARumCompGenEig<ARFLOAT>::
ARumCompGenEig(int nevp, ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B, char* whichp,
               int ncvp, ARFLOAT tolp, int maxitp,
               arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->NoShift();
  DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultInvBAv,
                   &Pencil, 
                   &ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultBv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARumCompGenEig<ARFLOAT>::
ARumCompGenEig(int nevp, ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B,
               arcomplex<ARFLOAT> sigmap, char* whichp, int ncvp,
               ARFLOAT tolp, int maxitp, arcomplex<ARFLOAT>* residp,
               bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultInvAsBv,
                   &Pencil, 
                   &ARumNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultBv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  SetShiftInvertMode(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARumCompGenEig<ARFLOAT>& ARumCompGenEig<ARFLOAT>::
operator=(const ARumCompGenEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUGCOMP_H
