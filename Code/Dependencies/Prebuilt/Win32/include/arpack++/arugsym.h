/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARUGSym.h.
   Arpack++ class ARumSymGenEig definition
   (UMFPACK version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Kristi Maschhoff
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARUGSYM_H
#define ARUGSYM_H

#include <cstddef>

#include "arch.h"
#include "arusmat.h"
#include "aruspen.h"
#include "argsym.h"


template<class ARFLOAT>
class ARumSymGenEig:
  public virtual ARSymGenEig<ARFLOAT, ARumSymPencil<ARFLOAT>,
                             ARumSymPencil<ARFLOAT> > {

 private:

 // a) Data structure used to store matrices.

  ARumSymPencil<ARFLOAT> Pencil;

 // b) Protected functions:

  virtual void Copy(const ARumSymGenEig& other);
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

  ARumSymGenEig() { }
  // Short constructor.

  ARumSymGenEig(int nevp, ARumSymMatrix<ARFLOAT>& A,
                ARumSymMatrix<ARFLOAT>& B, char* whichp = "LM",
                int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARumSymGenEig(char InvertModep, int nevp, ARumSymMatrix<ARFLOAT>& A,
                ARumSymMatrix<ARFLOAT>& B, ARFLOAT sigma, char* whichp = "LM", 
                int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert, buckling and Cayley modes).

  ARumSymGenEig(const ARumSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARumSymGenEig& operator=(const ARumSymGenEig& other);
  // Assignment operator.

}; // class ARumSymGenEig.


// ------------------------------------------------------------------------ //
// ARumSymGenEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARumSymGenEig<ARFLOAT>::
Copy(const ARumSymGenEig<ARFLOAT>& other)
{

  ARSymGenEig<ARFLOAT, ARumSymPencil<ARFLOAT>,
              ARumSymPencil<ARFLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  this->objOP  = &Pencil;
  this->objB   = &Pencil;
  this->objA   = &Pencil;

} // Copy.


template<class ARFLOAT>
inline void ARumSymGenEig<ARFLOAT>::ChangeShift(ARFLOAT sigmap)
{

  this->objOP->FactorAsB(sigmap);
  ARrcSymGenEig<ARFLOAT>::ChangeShift(sigmap);

} // ChangeShift.


template<class ARFLOAT>
inline void ARumSymGenEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARumSymPencil<ARFLOAT> >::
    SetRegularMode(&Pencil, &ARumSymPencil<ARFLOAT>::MultInvBAv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARumSymGenEig<ARFLOAT>::
SetShiftInvertMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARumSymPencil<ARFLOAT>, ARumSymPencil<ARFLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil, &ARumSymPencil<ARFLOAT>::MultInvAsBv);
  ChangeMultBx(&Pencil, &ARumSymPencil<ARFLOAT>::MultBv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline void ARumSymGenEig<ARFLOAT>::
SetBucklingMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARumSymPencil<ARFLOAT>, ARumSymPencil<ARFLOAT> >::
    SetBucklingMode(sigmap, &Pencil, &ARumSymPencil<ARFLOAT>::MultInvAsBv);
  ChangeMultBx(&Pencil, &ARumSymPencil<ARFLOAT>::MultAv);

} // SetBucklingMode.


template<class ARFLOAT>
inline void ARumSymGenEig<ARFLOAT>::
SetCayleyMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARumSymPencil<ARFLOAT>, ARumSymPencil<ARFLOAT> >::
    SetCayleyMode(sigmap, &Pencil, &ARumSymPencil<ARFLOAT>::MultInvAsBv,
                  &Pencil, &ARumSymPencil<ARFLOAT>::MultAv);
  ChangeMultBx(&Pencil, &ARumSymPencil<ARFLOAT>::MultBv);

} // SetCayleyMode.


template<class ARFLOAT>
inline ARumSymGenEig<ARFLOAT>::
ARumSymGenEig(int nevp, ARumSymMatrix<ARFLOAT>& A,
              ARumSymMatrix<ARFLOAT>& B, char* whichp, int ncvp,
              ARFLOAT tolp, int maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->InvertMode = 'S';
  this->NoShift();
  DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARumSymPencil<ARFLOAT>::MultInvBAv, &Pencil,
                   &ARumSymPencil<ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARumSymGenEig<ARFLOAT>::
ARumSymGenEig(char InvertModep, int nevp, ARumSymMatrix<ARFLOAT>& A,
              ARumSymMatrix<ARFLOAT>& B, ARFLOAT sigmap,
              char* whichp, int ncvp, ARFLOAT tolp,
              int maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARumSymPencil<ARFLOAT>::MultInvAsBv, &Pencil,
                   &ARumSymPencil<ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);
  this->InvertMode = this->CheckInvertMode(InvertModep);
  switch (this->InvertMode) {
  case 'B':  // Buckling mode.
    ChangeMultBx(&Pencil, &ARumSymPencil<ARFLOAT>::MultAv);
  case 'S':  // Shift and invert mode.
    ChangeShift(sigmap);
    break;  
  case 'C':  // Cayley mode.
    SetCayleyMode(sigmap);
  }

} // Long constructor (shift and invert, buckling and Cayley modes).


template<class ARFLOAT>
ARumSymGenEig<ARFLOAT>& ARumSymGenEig<ARFLOAT>::
operator=(const ARumSymGenEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
   this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUGSYM_H
