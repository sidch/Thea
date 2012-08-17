/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARDGComp.h.
   Arpack++ class ARdsCompGenEig definition
   (dense matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARDGCOMP_H
#define ARDGCOMP_H

#include <cstddef>

#include "arch.h"
#include "ardnsmat.h"
#include "ardnspen.h"
#include "arrseig.h"
#include "argcomp.h"


template<class ARFLOAT>
class ARdsCompGenEig:
  public virtual
    ARCompGenEig<ARFLOAT, ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT >,
                 ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT > > {

 private:

 // a) Data structure used to store matrices.

  ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT > Pencil;

 // b) Protected functions:

  virtual void Copy(const ARdsCompGenEig& other);
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

  ARdsCompGenEig() { }
  // Short constructor.

  ARdsCompGenEig(int nevp, ARdsNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 ARdsNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B,
                 char* whichp = "LM", int ncvp = 0,
                 ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARdsCompGenEig(int nevp, ARdsNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 ARdsNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B,
                 arcomplex<ARFLOAT> sigma, char* whichp = "LM",
                 int ncvp = 0, ARFLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARdsCompGenEig(const ARdsCompGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARdsCompGenEig() { }

 // d) Operators.

  ARdsCompGenEig& operator=(const ARdsCompGenEig& other);
  // Assignment operator.

}; // class ARdsCompGenEig.


// ------------------------------------------------------------------------ //
// ARdsCompGenEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARdsCompGenEig<ARFLOAT>::
Copy(const ARdsCompGenEig<ARFLOAT>& other)
{

  ARCompGenEig<ARFLOAT, ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT >,
               ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  this->objOP  = &Pencil;
  this->objB   = &Pencil;

} // Copy.


template<class ARFLOAT>
inline void ARdsCompGenEig<ARFLOAT>::
ChangeShift(arcomplex<ARFLOAT> sigmaRp)
{

  this->objOP->FactorAsB(sigmaRp);
  ARrcStdEig<ARFLOAT, arcomplex<ARFLOAT> >::ChangeShift(sigmaRp);

} // ChangeShift.


template<class ARFLOAT>
inline void ARdsCompGenEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>,
           ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetRegularMode(&Pencil,
                   &ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultInvBAv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARdsCompGenEig<ARFLOAT>::
SetShiftInvertMode(arcomplex<ARFLOAT> sigmap)
{

  ARCompGenEig<ARFLOAT, ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>,
               ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil,
                       &ARdsNonSymPencil<arcomplex<ARFLOAT>,ARFLOAT>::MultInvAsBv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARdsCompGenEig<ARFLOAT>::
ARdsCompGenEig(int nevp, ARdsNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               ARdsNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B, char* whichp,
               int ncvp, ARFLOAT tolp, int maxitp,
               arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->NoShift();
  DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultInvBAv,
                   &Pencil, 
                   &ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultBv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARdsCompGenEig<ARFLOAT>::
ARdsCompGenEig(int nevp, ARdsNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               ARdsNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& B,
               arcomplex<ARFLOAT> sigmap, char* whichp, int ncvp,
               ARFLOAT tolp, int maxitp, arcomplex<ARFLOAT>* residp,
               bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultInvAsBv,
                   &Pencil, 
                   &ARdsNonSymPencil<arcomplex<ARFLOAT>, ARFLOAT>::MultBv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  SetShiftInvertMode(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARdsCompGenEig<ARFLOAT>& ARdsCompGenEig<ARFLOAT>::
operator=(const ARdsCompGenEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDGCOMP_H
