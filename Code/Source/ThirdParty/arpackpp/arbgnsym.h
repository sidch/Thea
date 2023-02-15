/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARBGNSym.h.
   Arpack++ class ARbdNonSymGenEig definition
   (band matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARBGNSYM_H
#define ARBGNSYM_H

#include <cstddef>
#include "arch.h"
#include "arbnsmat.h"
#include "arbnspen.h"
#include "argnsym.h"


template<class ARFLOAT>
class ARbdNonSymGenEig:
  public virtual ARNonSymGenEig<ARFLOAT, ARbdNonSymPencil<ARFLOAT, ARFLOAT>,
                                ARbdNonSymPencil<ARFLOAT, ARFLOAT> > {

 private:

 // a) Data structure used to store matrices.

  ARbdNonSymPencil<ARFLOAT, ARFLOAT> Pencil;

 // b) Protected functions:

  virtual void Copy(const ARbdNonSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmaRp, ARFLOAT sigmaIp = 0.0);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

  virtual void SetComplexShiftMode(char partp, ARFLOAT sigmaRp, ARFLOAT sigmaIp);

 // c.2) Constructors and destructor.

  ARbdNonSymGenEig() { }
  // Short constructor.

  ARbdNonSymGenEig(int nevp, ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& B, const char* whichp = "LM",
                   int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARbdNonSymGenEig(int nevp, ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& B, ARFLOAT sigma,
                   const char* whichp = "LM", int ncvp = 0,
                   ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (real shift and invert mode).

  ARbdNonSymGenEig(int nevp, ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& B, char partp,
                   ARFLOAT sigmaRp, ARFLOAT sigmaIp, const char* whichp = "LM",
                   int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (complex shift and invert mode).

  ARbdNonSymGenEig(const ARbdNonSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARbdNonSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARbdNonSymGenEig& operator=(const ARbdNonSymGenEig& other);
  // Assignment operator.

}; // class ARbdNonSymGenEig.


// ------------------------------------------------------------------------ //
// ARbdNonSymGenEig member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARbdNonSymGenEig<ARFLOAT>::
Copy(const ARbdNonSymGenEig<ARFLOAT>& other)
{

  ARNonSymGenEig<ARFLOAT, ARbdNonSymPencil<ARFLOAT, ARFLOAT>,
                 ARbdNonSymPencil<ARFLOAT, ARFLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  this->objOP  = &Pencil;
  this->objB   = &Pencil;
  this->objA   = &Pencil;

} // Copy.


template<class ARFLOAT>
inline void ARbdNonSymGenEig<ARFLOAT>::
ChangeShift(ARFLOAT sigmaRp, ARFLOAT sigmaIp)
{

  if (sigmaIp == 0.0) {
    this->objOP->FactorAsB(sigmaRp);
  }
  else {
    this->objOP->FactorAsB(sigmaRp, sigmaIp, this->part);
  }
  ARrcNonSymGenEig<ARFLOAT>::ChangeShift(sigmaRp, sigmaIp);

} // ChangeShift.


template<class ARFLOAT>
inline void ARbdNonSymGenEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARbdNonSymPencil<ARFLOAT, ARFLOAT> >::
    SetRegularMode(&Pencil, &ARbdNonSymPencil<ARFLOAT, ARFLOAT>::MultInvBAv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARbdNonSymGenEig<ARFLOAT>::SetShiftInvertMode(ARFLOAT sigmap)
{

  ARNonSymGenEig<ARFLOAT, ARbdNonSymPencil<ARFLOAT, ARFLOAT>,
                 ARbdNonSymPencil<ARFLOAT, ARFLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil,
                       &ARbdNonSymPencil<ARFLOAT, ARFLOAT>::MultInvAsBv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline void ARbdNonSymGenEig<ARFLOAT>::
SetComplexShiftMode(char partp, ARFLOAT sigmaRp, ARFLOAT sigmaIp)
{

  ARNonSymGenEig<ARFLOAT, ARbdNonSymPencil<ARFLOAT, ARFLOAT>,
                 ARbdNonSymPencil<ARFLOAT, ARFLOAT> >::
    SetComplexShiftMode(partp, sigmaRp, sigmaIp, &Pencil,
                        &ARbdNonSymPencil<ARFLOAT, ARFLOAT>::MultInvAsBv,
                        &Pencil, &ARbdNonSymPencil<ARFLOAT, ARFLOAT>::MultAv);

} // SetComplexShiftMode.


template<class ARFLOAT>
inline ARbdNonSymGenEig<ARFLOAT>::
ARbdNonSymGenEig(int nevp, ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& B, const char* whichp, int ncvp,
                 ARFLOAT tolp, int maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARbdNonSymPencil<ARFLOAT, ARFLOAT>::MultInvBAv, &Pencil,
                   &ARbdNonSymPencil<ARFLOAT, ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARbdNonSymGenEig<ARFLOAT>::
ARbdNonSymGenEig(int nevp, ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& B, ARFLOAT sigmap,
                 const char* whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARbdNonSymPencil<ARFLOAT, ARFLOAT>::MultInvAsBv, &Pencil,
                   &ARbdNonSymPencil<ARFLOAT, ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);
  SetShiftInvertMode(sigmap);

} // Long constructor (real shift and invert mode).


template<class ARFLOAT>
inline ARbdNonSymGenEig<ARFLOAT>::
ARbdNonSymGenEig(int nevp, ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 ARbdNonSymMatrix<ARFLOAT, ARFLOAT>& B, char partp,
                 ARFLOAT sigmaRp, ARFLOAT sigmaIp, const char* whichp, int ncvp,
                 ARFLOAT tolp, int maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARbdNonSymPencil<ARFLOAT, ARFLOAT>::MultInvAsBv, &Pencil,
                   &ARbdNonSymPencil<ARFLOAT, ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);
  SetComplexShiftMode(partp, sigmaRp, sigmaIp);

} // Long constructor (complex shift and invert mode).


template<class ARFLOAT>
ARbdNonSymGenEig<ARFLOAT>& ARbdNonSymGenEig<ARFLOAT>::
operator=(const ARbdNonSymGenEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBGNSYM_H
