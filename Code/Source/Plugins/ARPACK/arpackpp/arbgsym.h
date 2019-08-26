/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARBGSym.h.
   Arpack++ class ARbdSymGenEig definition
   (band matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARBGSYM_H
#define ARBGSYM_H

#include <cstddef>
#include "arch.h"
#include "arbsmat.h"
#include "arbspen.h"
#include "argsym.h"


template<class ARFLOAT>
class ARbdSymGenEig:
  public virtual ARSymGenEig<ARFLOAT, ARbdSymPencil<ARFLOAT>,
                             ARbdSymPencil<ARFLOAT> > {

 private:

 // a) Data structure used to store matrices.

  ARbdSymPencil<ARFLOAT> Pencil;

 // b) Protected functions:

  virtual void Copy(const ARbdSymGenEig& other);
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

  ARbdSymGenEig() { }
  // Short constructor.

  ARbdSymGenEig(int nevp, ARbdSymMatrix<ARFLOAT>& A,
                ARbdSymMatrix<ARFLOAT>& B, const char* whichp = "LM",
                int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARbdSymGenEig(char InvertModep, int nevp, ARbdSymMatrix<ARFLOAT>& A,
                ARbdSymMatrix<ARFLOAT>& B, ARFLOAT sigma, const char* whichp = "LM",
                int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert, buckling and Cayley modes).

  ARbdSymGenEig(const ARbdSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARbdSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARbdSymGenEig& operator=(const ARbdSymGenEig& other);
  // Assignment operator.

}; // class ARbdSymGenEig.


// ------------------------------------------------------------------------ //
// ARbdSymGenEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARbdSymGenEig<ARFLOAT>::
Copy(const ARbdSymGenEig<ARFLOAT>& other)
{

  ARSymGenEig<ARFLOAT, ARbdSymPencil<ARFLOAT>,
              ARbdSymPencil<ARFLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  this->objOP  = &Pencil;
  this->objB   = &Pencil;
  this->objA   = &Pencil;

} // Copy.


template<class ARFLOAT>
inline void ARbdSymGenEig<ARFLOAT>::ChangeShift(ARFLOAT sigmap)
{

  this->objOP->FactorAsB(sigmap);
  ARrcSymGenEig<ARFLOAT>::ChangeShift(sigmap);

} // ChangeShift.


template<class ARFLOAT>
inline void ARbdSymGenEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARbdSymPencil<ARFLOAT> >::
    SetRegularMode(&Pencil, &ARbdSymPencil<ARFLOAT>::MultInvBAv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARbdSymGenEig<ARFLOAT>::
SetShiftInvertMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARbdSymPencil<ARFLOAT>, ARbdSymPencil<ARFLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil, &ARbdSymPencil<ARFLOAT>::MultInvAsBv);
  ChangeMultBx(&Pencil, &ARbdSymPencil<ARFLOAT>::MultBv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline void ARbdSymGenEig<ARFLOAT>::
SetBucklingMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARbdSymPencil<ARFLOAT>, ARbdSymPencil<ARFLOAT> >::
    SetBucklingMode(sigmap, &Pencil, &ARbdSymPencil<ARFLOAT>::MultInvAsBv);
  ChangeMultBx(&Pencil, &ARbdSymPencil<ARFLOAT>::MultAv);

} // SetBucklingMode.


template<class ARFLOAT>
inline void ARbdSymGenEig<ARFLOAT>::
SetCayleyMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARbdSymPencil<ARFLOAT>, ARbdSymPencil<ARFLOAT> >::
    SetCayleyMode(sigmap, &Pencil, &ARbdSymPencil<ARFLOAT>::MultInvAsBv,
                  &Pencil, &ARbdSymPencil<ARFLOAT>::MultAv);
  ChangeMultBx(&Pencil, &ARbdSymPencil<ARFLOAT>::MultBv);

} // SetCayleyMode.


template<class ARFLOAT>
inline ARbdSymGenEig<ARFLOAT>::
ARbdSymGenEig(int nevp, ARbdSymMatrix<ARFLOAT>& A,
              ARbdSymMatrix<ARFLOAT>& B, const char* whichp, int ncvp,
              ARFLOAT tolp, int maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->InvertMode = 'S';
  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARbdSymPencil<ARFLOAT>::MultInvBAv, &Pencil,
                   &ARbdSymPencil<ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARbdSymGenEig<ARFLOAT>::
ARbdSymGenEig(char InvertModep, int nevp, ARbdSymMatrix<ARFLOAT>& A,
              ARbdSymMatrix<ARFLOAT>& B, ARFLOAT sigmap,
              const char* whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARbdSymPencil<ARFLOAT>::MultInvAsBv, &Pencil,
                   &ARbdSymPencil<ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);
  this->InvertMode = this->CheckInvertMode(InvertModep);
  switch (this->InvertMode) {
  case 'B':
    ChangeMultBx(&Pencil, &ARbdSymPencil<ARFLOAT>::MultAv);
  case 'S':
    this->ChangeShift(sigmap);
    break;
  case 'C':
    SetCayleyMode(sigmap);
  }

} // Long constructor (shift and invert, buckling and Cayley modes).


template<class ARFLOAT>
ARbdSymGenEig<ARFLOAT>& ARbdSymGenEig<ARFLOAT>::
operator=(const ARbdSymGenEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBGSYM_H
