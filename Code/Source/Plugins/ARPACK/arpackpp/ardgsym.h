/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARDGSym.h.
   Arpack++ class ARdsSymGenEig definition
   (dense matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARDGSYM_H
#define ARDGSYM_H

#include <cstddef>

#include "arch.h"
#include "ardsmat.h"
#include "ardspen.h"
#include "argsym.h"


template<class ARFLOAT>
class ARdsSymGenEig:
  public virtual ARSymGenEig<ARFLOAT, ARdsSymPencil<ARFLOAT>,
                             ARdsSymPencil<ARFLOAT> > {

 private:

 // a) Data structure used to store matrices.

  ARdsSymPencil<ARFLOAT> Pencil;

 // b) Protected functions:

  virtual void Copy(const ARdsSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmap);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

  virtual void SetBucklingMode(ARFLOAT sigmap);

  virtual void SetCayleyMode(ARFLOAT sigmap);

 // c.2) Constructors and destructor.

  ARdsSymGenEig() { }
  // Short constructor.

  ARdsSymGenEig(int nevp, ARdsSymMatrix<ARFLOAT>& A,
                ARdsSymMatrix<ARFLOAT>& B, const char* whichp = "LM",
                int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARdsSymGenEig(char InvertModep, int nevp, ARdsSymMatrix<ARFLOAT>& A,
                ARdsSymMatrix<ARFLOAT>& B, ARFLOAT sigma, const char* whichp = "LM",
                int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert, buckling and Cayley modes).

  ARdsSymGenEig(const ARdsSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARdsSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARdsSymGenEig& operator=(const ARdsSymGenEig& other);
  // Assignment operator.

}; // class ARdsSymGenEig.


// ------------------------------------------------------------------------ //
// ARdsSymGenEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARdsSymGenEig<ARFLOAT>::
Copy(const ARdsSymGenEig<ARFLOAT>& other)
{

  ARSymGenEig<ARFLOAT, ARdsSymPencil<ARFLOAT>,
              ARdsSymPencil<ARFLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  this->objOP  = &Pencil;
  this->objB   = &Pencil;
  this->objA   = &Pencil;

} // Copy.


template<class ARFLOAT>
inline void ARdsSymGenEig<ARFLOAT>::ChangeShift(ARFLOAT sigmap)
{

  this->objOP->FactorAsB(sigmap);
  ARrcSymGenEig<ARFLOAT>::ChangeShift(sigmap);

} // ChangeShift.


template<class ARFLOAT>
inline void ARdsSymGenEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARdsSymPencil<ARFLOAT> >::
    SetRegularMode(&Pencil, &ARdsSymPencil<ARFLOAT>::MultInvBAv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARdsSymGenEig<ARFLOAT>::
SetShiftInvertMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARdsSymPencil<ARFLOAT>, ARdsSymPencil<ARFLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil, &ARdsSymPencil<ARFLOAT>::MultInvAsBv);
  ChangeMultBx(&Pencil, &ARdsSymPencil<ARFLOAT>::MultBv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline void ARdsSymGenEig<ARFLOAT>::
SetBucklingMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARdsSymPencil<ARFLOAT>, ARdsSymPencil<ARFLOAT> >::
    SetBucklingMode(sigmap, &Pencil, &ARdsSymPencil<ARFLOAT>::MultInvAsBv);
  ChangeMultBx(&Pencil, &ARdsSymPencil<ARFLOAT>::MultAv);

} // SetBucklingMode.


template<class ARFLOAT>
inline void ARdsSymGenEig<ARFLOAT>::
SetCayleyMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARdsSymPencil<ARFLOAT>, ARdsSymPencil<ARFLOAT> >::
    SetCayleyMode(sigmap, &Pencil, &ARdsSymPencil<ARFLOAT>::MultInvAsBv,
                  &Pencil, &ARdsSymPencil<ARFLOAT>::MultAv);
  ChangeMultBx(&Pencil, &ARdsSymPencil<ARFLOAT>::MultBv);

} // SetCayleyMode.


template<class ARFLOAT>
inline ARdsSymGenEig<ARFLOAT>::
ARdsSymGenEig(int nevp, ARdsSymMatrix<ARFLOAT>& A,
              ARdsSymMatrix<ARFLOAT>& B, const char* whichp, int ncvp,
              ARFLOAT tolp, int maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->InvertMode = 'S';
  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARdsSymPencil<ARFLOAT>::MultInvBAv, &Pencil,
                   &ARdsSymPencil<ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARdsSymGenEig<ARFLOAT>::
ARdsSymGenEig(char InvertModep, int nevp, ARdsSymMatrix<ARFLOAT>& A,
              ARdsSymMatrix<ARFLOAT>& B, ARFLOAT sigmap,
              const char* whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARdsSymPencil<ARFLOAT>::MultInvAsBv, &Pencil,
                   &ARdsSymPencil<ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);
  this->InvertMode = this->CheckInvertMode(InvertModep);
  switch (this->InvertMode) {
  case 'B':
    ChangeMultBx(&Pencil, &ARdsSymPencil<ARFLOAT>::MultAv);
  case 'S':
    this->ChangeShift(sigmap);
    break;
  case 'C':
    SetCayleyMode(sigmap);
  }

} // Long constructor (shift and invert, buckling and Cayley modes).


template<class ARFLOAT>
ARdsSymGenEig<ARFLOAT>& ARdsSymGenEig<ARFLOAT>::
operator=(const ARdsSymGenEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDGSYM_H
