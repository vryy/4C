/*----------------------------------------------------------------------*/
/*!
\file material.H

\brief Interface class for complex materials at Gauss points

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef MATERIAL_H
#define MATERIAL_H

#ifdef CCADISCRET

#include <Teuchos_RefCountPtr.hpp>
#include "../drt_lib/drt_parobject.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C"
{
#include "../headers/standardtypes.h"
}


/// MAT: materials
namespace MAT
{

/// Interface class for materials in elements
/*!

   The Material class provides a general way to store material history
   variables in elements. Each element owns one object derived from
   this class. Each derived class implements a specific material and
   might even store history data.

   \note There is one material object to each element, so one material
   object must store the history data of all Gauss points.

   Different material classes will (in general) not share a common set
   of functions, so the element needs to cast to the right
   material implementation. The element code will typically perform a
   test on MaterialType() and do the appropriate cast afterwards.

   <h3>Storage</h3>

   In order to support storage of material data (for restart, post
   processing, ...) this class does implement DRT::ParObject. The
   elements pack and unpack its materials along with its other data.

   \author u. kue and l.w.
   \date 06/07
*/
class Material : public DRT::ParObject
{
public:

  /// return type of this material
  virtual MATERIAL_TYP MaterialType() const = 0;

  /// create element material object given the number of a material definition
  static Teuchos::RefCountPtr<Material> Factory(int matnum);
};

  /*!
  \brief Multiply two 2nd order tensors A x B and add the result to a
   4th order symmetric material tensor in matrix notation

  \param C (in/out): Material tangent matrix to be modified
  \param AB (in)   : Scalar to multiply with A x B
  \param A (in)    : Dense matrix (3 by 3) as 2nd order tensor A
  \param B (in)    : Dense matrix (3 by 3) as 2nd order tensor B
  \param ScalarThis(in): Scalar to multiply with C before adding A x B
  */
  void ElastSymTensorMultiply(Epetra_SerialDenseMatrix& C,
                                   const double ScalarAB,
                                   const Epetra_SerialDenseMatrix& A,
                                   const Epetra_SerialDenseMatrix& B,
                                   const double ScalarThis);

  /*!
  \brief Multiply two 2nd order tensors (A x B + B x A) and add the result to a
   4th order symmetric material tensor in matrix notation

  \param C (in/out): Material tangent matrix to be modified
  \param AB (in)   : Scalar to multiply with (A x B + B x A)
  \param A (in)    : Dense matrix (3 by 3) as 2nd order tensor A
  \param B (in)    : Dense matrix (3 by 3) as 2nd order tensor B
  \param ScalarThis(in): Scalar to multiply with C before adding (A x B + B x A)
  */
  void ElastSymTensorMultiplyAddSym(Epetra_SerialDenseMatrix& C,
                                   const double ScalarAB,
                                   const Epetra_SerialDenseMatrix& A,
                                   const Epetra_SerialDenseMatrix& B,
                                   const double ScalarThis);

  /*!
  \brief Multiply two 2nd order tensors A o B and add the result to a
   4th order symmetric material tensor in matrix notation

  \param C (in/out): Material tangent matrix to be modified
  \param AB (in)   : Scalar to multiply with A o B
  \param A (in)    : Dense matrix (3 by 3) as 2nd order tensor A
  \param B (in)    : Dense matrix (3 by 3) as 2nd order tensor B
  \param ScalarThis(in): Scalar to multiply with C before adding A o B
  */
  void ElastSymTensor_o_Multiply(Epetra_SerialDenseMatrix& C,
                                   const double ScalarAB,
                                   const Epetra_SerialDenseMatrix& A,
                                   const Epetra_SerialDenseMatrix& B,
                                   const double ScalarThis);



} // namespace MAT

#endif
#endif
