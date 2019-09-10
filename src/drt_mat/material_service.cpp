/*----------------------------------------------------------------------*/
/*! \file
\brief Interface class for complex materials at Gauss points

\level 1

\maintainer Amadeus Gebauer
*/
/*----------------------------------------------------------------------*/


//#include "../drt_lib/drt_globalproblem.H"
#include "material_service.H"
#include "matpar_parameter.H"
#include "matpar_bundle.H"
#include "../linalg/linalg_utils.H"

#include <Sacado.hpp>

typedef Sacado::Fad::DFad<double> FAD;


/*----------------------------------------------------------------------*
 |  Add 'Holzapfel product' contribution to constitutive tensor         |
 |  using Voigt notation                                                |
 |                                                (public)  chfoe 04/08 |
 *----------------------------------------------------------------------*

 This function adds the following contribution to the given constitutive
 matrix cmat(6,6) based on the inverse of the right Cauchy-Green vector
 invc(6):

 scalar * ( Cinv boeppel Cinv )

 For that purpose we need the derivative

  \partial tensor(C)^-1
 -----------------------
   \partial tensor(C)

 which yields the following product

  - ( Cinv boeppel Cinv )_{abcd} = 1/2 * ( Cinv_{ac} Cinv_{bd} + Cinv_{ad} Cinv_{bc} )

 For more details see Holzapfel p. 254

 */
void MAT::AddtoCmatHolzapfelProduct(
    Epetra_SerialDenseMatrix& cmat, const Epetra_SerialDenseVector& invc, const double scalar)
{
#ifdef DEBUG
  if (cmat.M() != 6 or cmat.N() != 6 or invc.Length() != 6)
    dserror("Wrong dimensions in function AddtoCmatHolzapfelProduct");
#endif

  // and the 'boeppel-product' for the expression d(invc)/dc (see Holzapfel p. 254)
  cmat(0, 0) += scalar * invc(0) * invc(0);
  cmat(0, 1) += scalar * invc(3) * invc(3);
  cmat(0, 2) += scalar * invc(5) * invc(5);
  cmat(0, 3) += scalar * invc(0) * invc(3);
  cmat(0, 4) += scalar * invc(3) * invc(5);
  cmat(0, 5) += scalar * invc(0) * invc(5);

  cmat(1, 0) += scalar * invc(3) * invc(3);
  cmat(1, 1) += scalar * invc(1) * invc(1);
  cmat(1, 2) += scalar * invc(4) * invc(4);
  cmat(1, 3) += scalar * invc(3) * invc(1);
  cmat(1, 4) += scalar * invc(1) * invc(4);
  cmat(1, 5) += scalar * invc(3) * invc(4);

  cmat(2, 0) += scalar * invc(5) * invc(5);
  cmat(2, 1) += scalar * invc(4) * invc(4);
  cmat(2, 2) += scalar * invc(2) * invc(2);
  cmat(2, 3) += scalar * invc(5) * invc(4);
  cmat(2, 4) += scalar * invc(4) * invc(2);
  cmat(2, 5) += scalar * invc(5) * invc(2);

  cmat(3, 0) += scalar * invc(0) * invc(3);
  cmat(3, 1) += scalar * invc(3) * invc(1);
  cmat(3, 2) += scalar * invc(5) * invc(4);
  cmat(3, 3) += scalar * 0.5 * (invc(0) * invc(1) + invc(3) * invc(3));
  cmat(3, 4) += scalar * 0.5 * (invc(3) * invc(4) + invc(5) * invc(1));
  cmat(3, 5) += scalar * 0.5 * (invc(0) * invc(4) + invc(5) * invc(3));

  cmat(4, 0) += scalar * invc(3) * invc(5);
  cmat(4, 1) += scalar * invc(1) * invc(4);
  cmat(4, 2) += scalar * invc(4) * invc(2);
  cmat(4, 3) += scalar * 0.5 * (invc(3) * invc(4) + invc(5) * invc(1));
  cmat(4, 4) += scalar * 0.5 * (invc(1) * invc(2) + invc(4) * invc(4));
  cmat(4, 5) += scalar * 0.5 * (invc(3) * invc(2) + invc(4) * invc(5));

  cmat(5, 0) += scalar * invc(0) * invc(5);
  cmat(5, 1) += scalar * invc(3) * invc(4);
  cmat(5, 2) += scalar * invc(5) * invc(2);
  cmat(5, 3) += scalar * 0.5 * (invc(0) * invc(4) + invc(5) * invc(3));
  cmat(5, 4) += scalar * 0.5 * (invc(3) * invc(2) + invc(4) * invc(5));
  cmat(5, 5) += scalar * 0.5 * (invc(0) * invc(2) + invc(5) * invc(5));

  return;
}

/*----------------------------------------------------------------------*
 |  Add 'Holzapfel product' contribution to constitutive tensor         |
 |  using Voigt notation                                                |
 | This is a plain copy of the Epetra version of this method            |
 | with different parameter types
 |                                                (public)  mgee  10/08 |
 *----------------------------------------------------------------------*

 This function adds the following contribution to the given constitutive
 matrix cmat(6,6) based on the inverse of the right Cauchy-Green vector
 invc(6):

 scalar * ( Cinv boeppel Cinv )

 For that purpose we need the derivative

  \partial tensor(C)^-1
 -----------------------
   \partial tensor(C)

 which yields the following product

  - ( Cinv boeppel Cinv )_{abcd} = - 1/2 * ( Cinv_{ac} Cinv_{bd} + Cinv_{ad} Cinv_{bc} )

 For more details see Holzapfel p. 254

 */
template <typename T>
void MAT::AddtoCmatHolzapfelProduct(
    LINALG::Matrix<6, 6, T>& cmat, const LINALG::Matrix<6, 1, T>& invc, const T scalar)
{
#ifdef DEBUG
  if (cmat.M() != 6 or cmat.N() != 6 or invc.M() != 6)
    dserror("Wrong dimensions in function AddtoCmatHolzapfelProduct");
#endif

  // and the 'boeppel-product' for the expression d(invc)/dc (see Holzapfel p. 254)
  cmat(0, 0) += scalar * invc(0) * invc(0);
  cmat(0, 1) += scalar * invc(3) * invc(3);
  cmat(0, 2) += scalar * invc(5) * invc(5);
  cmat(0, 3) += scalar * invc(0) * invc(3);
  cmat(0, 4) += scalar * invc(3) * invc(5);
  cmat(0, 5) += scalar * invc(0) * invc(5);

  cmat(1, 0) += scalar * invc(3) * invc(3);
  cmat(1, 1) += scalar * invc(1) * invc(1);
  cmat(1, 2) += scalar * invc(4) * invc(4);
  cmat(1, 3) += scalar * invc(3) * invc(1);
  cmat(1, 4) += scalar * invc(1) * invc(4);
  cmat(1, 5) += scalar * invc(3) * invc(4);

  cmat(2, 0) += scalar * invc(5) * invc(5);
  cmat(2, 1) += scalar * invc(4) * invc(4);
  cmat(2, 2) += scalar * invc(2) * invc(2);
  cmat(2, 3) += scalar * invc(5) * invc(4);
  cmat(2, 4) += scalar * invc(4) * invc(2);
  cmat(2, 5) += scalar * invc(5) * invc(2);

  cmat(3, 0) += scalar * invc(0) * invc(3);
  cmat(3, 1) += scalar * invc(3) * invc(1);
  cmat(3, 2) += scalar * invc(5) * invc(4);
  cmat(3, 3) += scalar * 0.5 * (invc(0) * invc(1) + invc(3) * invc(3));
  cmat(3, 4) += scalar * 0.5 * (invc(3) * invc(4) + invc(5) * invc(1));
  cmat(3, 5) += scalar * 0.5 * (invc(0) * invc(4) + invc(5) * invc(3));

  cmat(4, 0) += scalar * invc(3) * invc(5);
  cmat(4, 1) += scalar * invc(1) * invc(4);
  cmat(4, 2) += scalar * invc(4) * invc(2);
  cmat(4, 3) += scalar * 0.5 * (invc(3) * invc(4) + invc(5) * invc(1));
  cmat(4, 4) += scalar * 0.5 * (invc(1) * invc(2) + invc(4) * invc(4));
  cmat(4, 5) += scalar * 0.5 * (invc(3) * invc(2) + invc(4) * invc(5));

  cmat(5, 0) += scalar * invc(0) * invc(5);
  cmat(5, 1) += scalar * invc(3) * invc(4);
  cmat(5, 2) += scalar * invc(5) * invc(2);
  cmat(5, 3) += scalar * 0.5 * (invc(0) * invc(4) + invc(5) * invc(3));
  cmat(5, 4) += scalar * 0.5 * (invc(3) * invc(2) + invc(4) * invc(5));
  cmat(5, 5) += scalar * 0.5 * (invc(0) * invc(2) + invc(5) * invc(5));

  return;
}


/*----------------------------------------------------------------------*
 | compute the "elasticity tensor product" A x B of                     |
 | two 2nd order tensors (in matrix notation) and add the result to     |
 | a 4th order tensor (in Voigt matrix notation!) using the             |
 | symmetry-conditions inherent to elasticity tensors                   |
 | The implementation is based on the Epetra-Method Matrix.Multiply.    |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
void MAT::ElastSymTensorMultiply(Epetra_SerialDenseMatrix& C, const double ScalarAB,
    const Epetra_SerialDenseMatrix& A, const Epetra_SerialDenseMatrix& B, const double ScalarThis)
{
#ifdef DEBUG
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3)
  {
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");
#endif

  // everything in Voigt-Notation
  LINALG::SerialDenseMatrix AVoigt(6, 1);
  LINALG::SerialDenseMatrix BVoigt(6, 1);

  AVoigt(0, 0) = A(0, 0);
  AVoigt(1, 0) = A(1, 1);
  AVoigt(2, 0) = A(2, 2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3, 0) = A(1, 0);
  AVoigt(4, 0) = A(2, 1);
  AVoigt(5, 0) = A(2, 0);

  BVoigt(0, 0) = B(0, 0);
  BVoigt(1, 0) = B(1, 1);
  BVoigt(2, 0) = B(2, 2);
  BVoigt(3, 0) = B(1, 0);
  BVoigt(4, 0) = B(2, 1);
  BVoigt(5, 0) = B(2, 0);

  C.Multiply('N', 'T', ScalarAB, AVoigt, BVoigt, ScalarThis);

  // this is explicitly what the former .Multiply does:
  //  C(0,0)= ScalarThis*C(0,0) + ScalarAB * A(0,0)*B(0,0);
  //  C(0,1)= ScalarThis*C(0,1) + ScalarAB * A(0,0)*B(1,1);
  //  C(0,2)= ScalarThis*C(0,2) + ScalarAB * A(0,0)*B(2,2);
  //  C(0,3)= ScalarThis*C(0,3) + ScalarAB * A(0,0)*B(1,0);
  //  C(0,4)= ScalarThis*C(0,4) + ScalarAB * A(0,0)*B(2,1);
  //  C(0,5)= ScalarThis*C(0,5) + ScalarAB * A(0,0)*B(2,0);
  //
  //  C(1,0)= ScalarThis*C(1,0) + ScalarAB * A(1,1)*B(0,0);
  //  C(1,1)= ScalarThis*C(1,1) + ScalarAB * A(1,1)*B(1,1);
  //  C(1,2)= ScalarThis*C(1,2) + ScalarAB * A(1,1)*B(2,2);
  //  C(1,3)= ScalarThis*C(1,3) + ScalarAB * A(1,1)*B(1,0);
  //  C(1,4)= ScalarThis*C(1,4) + ScalarAB * A(1,1)*B(2,1);
  //  C(1,5)= ScalarThis*C(1,5) + ScalarAB * A(1,1)*B(2,0);
  //
  //  C(2,0)= ScalarThis*C(2,0) + ScalarAB * A(2,2)*B(0,0);
  //  C(2,1)= ScalarThis*C(2,1) + ScalarAB * A(2,2)*B(1,1);
  //  C(2,2)= ScalarThis*C(2,2) + ScalarAB * A(2,2)*B(2,2);
  //  C(2,3)= ScalarThis*C(2,3) + ScalarAB * A(2,2)*B(1,0);
  //  C(2,4)= ScalarThis*C(2,4) + ScalarAB * A(2,2)*B(2,1);
  //  C(2,5)= ScalarThis*C(2,5) + ScalarAB * A(2,2)*B(2,0);
  //
  //  C(3,0)= ScalarThis*C(3,0) + ScalarAB * A(1,0)*B(0,0);
  //  C(3,1)= ScalarThis*C(3,1) + ScalarAB * A(1,0)*B(1,1);
  //  C(3,2)= ScalarThis*C(3,2) + ScalarAB * A(1,0)*B(2,2);
  //  C(3,3)= ScalarThis*C(3,3) + ScalarAB * A(1,0)*B(1,0);
  //  C(3,4)= ScalarThis*C(3,4) + ScalarAB * A(1,0)*B(2,1);
  //  C(3,5)= ScalarThis*C(3,5) + ScalarAB * A(1,0)*B(2,0);
  //
  //  C(4,0)= ScalarThis*C(4,0) + ScalarAB * A(2,1)*B(0,0);
  //  C(4,1)= ScalarThis*C(4,1) + ScalarAB * A(2,1)*B(1,1);
  //  C(4,2)= ScalarThis*C(4,2) + ScalarAB * A(2,1)*B(2,2);
  //  C(4,3)= ScalarThis*C(4,3) + ScalarAB * A(2,1)*B(1,0);
  //  C(4,4)= ScalarThis*C(4,4) + ScalarAB * A(2,1)*B(2,1);
  //  C(4,5)= ScalarThis*C(4,5) + ScalarAB * A(2,1)*B(2,0);
  //
  //  C(5,0)= ScalarThis*C(5,0) + ScalarAB * A(2,0)*B(0,0);
  //  C(5,1)= ScalarThis*C(5,1) + ScalarAB * A(2,0)*B(1,1);
  //  C(5,2)= ScalarThis*C(5,2) + ScalarAB * A(2,0)*B(2,2);
  //  C(5,3)= ScalarThis*C(5,3) + ScalarAB * A(2,0)*B(1,0);
  //  C(5,4)= ScalarThis*C(5,4) + ScalarAB * A(2,0)*B(2,1);
  //  C(5,5)= ScalarThis*C(5,5) + ScalarAB * A(2,0)*B(2,0);

  return;
}


/*----------------------------------------------------------------------*
 | compute the "elasticity tensor product" A x B of                     |
 | two 2nd order tensors (in matrix notation) and add the result to     |
 | a 4th order tensor (in Voigt matrix notation!) using the             |
 | symmetry-conditions inherent to elasticity tensors                   |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
// This is a copy of the above function, using the fixed size matrix.
void MAT::ElastSymTensorMultiply(LINALG::Matrix<6, 6>& C, const double ScalarAB,
    const LINALG::Matrix<3, 3>& A, const LINALG::Matrix<3, 3>& B, const double ScalarThis)
{
  // everything in Voigt-Notation
  LINALG::Matrix<6, 1> AVoigt;
  LINALG::Matrix<6, 1> BVoigt;

  AVoigt(0, 0) = A(0, 0);
  AVoigt(1, 0) = A(1, 1);
  AVoigt(2, 0) = A(2, 2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3, 0) = A(1, 0);
  AVoigt(4, 0) = A(2, 1);
  AVoigt(5, 0) = A(2, 0);

  BVoigt(0, 0) = B(0, 0);
  BVoigt(1, 0) = B(1, 1);
  BVoigt(2, 0) = B(2, 2);
  BVoigt(3, 0) = B(1, 0);
  BVoigt(4, 0) = B(2, 1);
  BVoigt(5, 0) = B(2, 0);

  C.MultiplyNT(ScalarAB, AVoigt, BVoigt, ScalarThis);

  return;
}


/*----------------------------------------------------------------------*
 | compute the "elasticity tensor product" (A x B + B x A) of           |
 | two 2nd order tensors (in matrix notation) and add the result to     |
 | a 4th order tensor (in Voigt matrix notation!) using the             |
 | symmetry-conditions inherent to elasticity tensors                   |
 | The implementation is based on the Epetra-Method Matrix.Multiply.    |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
void MAT::ElastSymTensorMultiplyAddSym(Epetra_SerialDenseMatrix& C, const double ScalarAB,
    const Epetra_SerialDenseMatrix& A, const Epetra_SerialDenseMatrix& B, const double ScalarThis)
{
#ifdef DEBUG
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3)
  {
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");
#endif

  // everything in Voigt-Notation
  LINALG::SerialDenseMatrix AVoigt(6, 1);
  LINALG::SerialDenseMatrix BVoigt(6, 1);

  AVoigt(0, 0) = A(0, 0);
  AVoigt(1, 0) = A(1, 1);
  AVoigt(2, 0) = A(2, 2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3, 0) = A(1, 0);
  AVoigt(4, 0) = A(2, 1);
  AVoigt(5, 0) = A(2, 0);

  BVoigt(0, 0) = B(0, 0);
  BVoigt(1, 0) = B(1, 1);
  BVoigt(2, 0) = B(2, 2);
  BVoigt(3, 0) = B(1, 0);
  BVoigt(4, 0) = B(2, 1);
  BVoigt(5, 0) = B(2, 0);

  C.Multiply('N', 'T', ScalarAB, AVoigt, BVoigt, ScalarThis);
  C.Multiply('N', 'T', ScalarAB, BVoigt, AVoigt, 1.0);

  // this is explicitly what the former .Multiplies do:
  //  C(0,0)= ScalarThis*C(0,0) + ScalarAB * (A(0,0)*B(0,0) + B(0,0)*A(0,0));
  //  C(0,1)= ScalarThis*C(0,1) + ScalarAB * (A(0,0)*B(1,1) + B(0,0)*A(1,1));
  //  C(0,2)= ScalarThis*C(0,2) + ScalarAB * (A(0,0)*B(2,2) + B(0,0)*A(2,2));
  //  C(0,3)= ScalarThis*C(0,3) + ScalarAB * (A(0,0)*B(1,0) + B(0,0)*A(1,0));
  //  C(0,4)= ScalarThis*C(0,4) + ScalarAB * (A(0,0)*B(2,1) + B(0,0)*A(2,1));
  //  C(0,5)= ScalarThis*C(0,5) + ScalarAB * (A(0,0)*B(2,0) + B(0,0)*A(2,0));
  //
  //  C(1,0)= ScalarThis*C(1,0) + ScalarAB * (A(1,1)*B(0,0) + B(1,1)*A(0,0));
  //  C(1,1)= ScalarThis*C(1,1) + ScalarAB * (A(1,1)*B(1,1) + B(1,1)*A(1,1));
  //  C(1,2)= ScalarThis*C(1,2) + ScalarAB * (A(1,1)*B(2,2) + B(1,1)*A(2,2));
  //  C(1,3)= ScalarThis*C(1,3) + ScalarAB * (A(1,1)*B(1,0) + B(1,1)*A(1,0));
  //  C(1,4)= ScalarThis*C(1,4) + ScalarAB * (A(1,1)*B(2,1) + B(1,1)*A(2,1));
  //  C(1,5)= ScalarThis*C(1,5) + ScalarAB * (A(1,1)*B(2,0) + B(1,1)*A(2,0));
  //
  //  C(2,0)= ScalarThis*C(2,0) + ScalarAB * (A(2,2)*B(0,0) + B(2,2)*A(0,0));
  //  C(2,1)= ScalarThis*C(2,1) + ScalarAB * (A(2,2)*B(1,1) + B(2,2)*A(1,1));
  //  C(2,2)= ScalarThis*C(2,2) + ScalarAB * (A(2,2)*B(2,2) + B(2,2)*A(2,2));
  //  C(2,3)= ScalarThis*C(2,3) + ScalarAB * (A(2,2)*B(1,0) + B(2,2)*A(1,0));
  //  C(2,4)= ScalarThis*C(2,4) + ScalarAB * (A(2,2)*B(2,1) + B(2,2)*A(2,1));
  //  C(2,5)= ScalarThis*C(2,5) + ScalarAB * (A(2,2)*B(2,0) + B(2,2)*A(2,0));
  //
  //  C(3,0)= ScalarThis*C(3,0) + ScalarAB * (A(1,0)*B(0,0) + B(1,0)*A(0,0));
  //  C(3,1)= ScalarThis*C(3,1) + ScalarAB * (A(1,0)*B(1,1) + B(1,0)*A(1,1));
  //  C(3,2)= ScalarThis*C(3,2) + ScalarAB * (A(1,0)*B(2,2) + B(1,0)*A(2,2));
  //  C(3,3)= ScalarThis*C(3,3) + ScalarAB * (A(1,0)*B(1,0) + B(1,0)*A(1,0));
  //  C(3,4)= ScalarThis*C(3,4) + ScalarAB * (A(1,0)*B(2,1) + B(1,0)*A(2,1));
  //  C(3,5)= ScalarThis*C(3,5) + ScalarAB * (A(1,0)*B(2,0) + B(1,0)*A(2,0));
  //
  //  C(4,0)= ScalarThis*C(4,0) + ScalarAB * (A(2,1)*B(0,0) + B(2,1)*A(0,0));
  //  C(4,1)= ScalarThis*C(4,1) + ScalarAB * (A(2,1)*B(1,1) + B(2,1)*A(1,1));
  //  C(4,2)= ScalarThis*C(4,2) + ScalarAB * (A(2,1)*B(2,2) + B(2,1)*A(2,2));
  //  C(4,3)= ScalarThis*C(4,3) + ScalarAB * (A(2,1)*B(1,0) + B(2,1)*A(1,0));
  //  C(4,4)= ScalarThis*C(4,4) + ScalarAB * (A(2,1)*B(2,1) + B(2,1)*A(2,1));
  //  C(4,5)= ScalarThis*C(4,5) + ScalarAB * (A(2,1)*B(2,0) + B(2,1)*A(2,0));
  //
  //  C(5,0)= ScalarThis*C(5,0) + ScalarAB * (A(2,0)*B(0,0) + B(2,0)*A(0,0));
  //  C(5,1)= ScalarThis*C(5,1) + ScalarAB * (A(2,0)*B(1,1) + B(2,0)*A(1,1));
  //  C(5,2)= ScalarThis*C(5,2) + ScalarAB * (A(2,0)*B(2,2) + B(2,0)*A(2,2));
  //  C(5,3)= ScalarThis*C(5,3) + ScalarAB * (A(2,0)*B(1,0) + B(2,0)*A(1,0));
  //  C(5,4)= ScalarThis*C(5,4) + ScalarAB * (A(2,0)*B(2,1) + B(2,0)*A(2,1));
  //  C(5,5)= ScalarThis*C(5,5) + ScalarAB * (A(2,0)*B(2,0) + B(2,0)*A(2,0));

  return;
}

/*----------------------------------------------------------------------*
 | compute the "elasticity tensor product" (A x B + B x A) of           |
 | two 2nd order tensors (in matrix notation) and add the result to     |
 | a 4th order tensor (in Voigt matrix notation!) using the             |
 | symmetry-conditions inherent to elasticity tensors                   |
 | The implementation is based on the Epetra-Method Matrix.Multiply.    |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
// This is a copy of the above function, using the fixed size matrix.
void MAT::ElastSymTensorMultiplyAddSym(LINALG::Matrix<6, 6>& C, const double ScalarAB,
    const LINALG::Matrix<3, 3>& A, const LINALG::Matrix<3, 3>& B, const double ScalarThis)
{
#ifdef DEBUG
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3)
  {
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");
#endif

  // everything in Voigt-Notation
  LINALG::Matrix<6, 1> AVoigt;
  LINALG::Matrix<6, 1> BVoigt;

  AVoigt(0, 0) = A(0, 0);
  AVoigt(1, 0) = A(1, 1);
  AVoigt(2, 0) = A(2, 2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3, 0) = A(1, 0);
  AVoigt(4, 0) = A(2, 1);
  AVoigt(5, 0) = A(2, 0);

  BVoigt(0, 0) = B(0, 0);
  BVoigt(1, 0) = B(1, 1);
  BVoigt(2, 0) = B(2, 2);
  BVoigt(3, 0) = B(1, 0);
  BVoigt(4, 0) = B(2, 1);
  BVoigt(5, 0) = B(2, 0);

  C.MultiplyNT(ScalarAB, AVoigt, BVoigt, ScalarThis);
  C.MultiplyNT(ScalarAB, BVoigt, AVoigt, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 | compute the "material tensor product" A o B (also known as           |
 | kronecker-tensor-product) of two 2nd order tensors                   |
 | (in matrix notation) and add the result to a 4th order tensor        |
 | (also in matrix notation) using the symmetry-conditions inherent to  |
 | material tensors, or tangent matrices, respectively                  |
 | AND the Voigt notation of E,S, and C with the famous factor 2!       |
 | The implementation is based on the Epetra-Method Matrix.Multiply.    |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
void MAT::ElastSymTensor_o_Multiply(Epetra_SerialDenseMatrix& C, const double ScalarAB,
    const Epetra_SerialDenseMatrix& A, const Epetra_SerialDenseMatrix& B, const double ScalarThis)
{
#ifdef DEBUG
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3)
  {
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");
#endif

  /* the kronecker-product in matrix notation is:
   * A11*B11 A11*B12 A11*B13   A12*B11 A12*B12 A12*B13   A13*B11 A13*B12 A13*B13
   * A11*B21 ...
   * A11*B31 ...
   *
   * A21*B11
   * A21*B21
   * A21*B31
   * ...                                                 A33*B11 A33*B12 A33*B13
   *                                                     A33*B21 A33*B22 A33*B23
   *                                                     A33*B31 A33*B32 A33*B33
   */
  /* to reduce the resulting 9by9 matrix to 6by6 we refer to the
   * Diss. from Balzani, Anhang D, BUT
   * we consider a factor 2 for colums/rows 4-6 :
   *  C(1)               2* 1/2*(C(2)+C(3))
   *  2* 1/2*(C(2)+C(3)  2* 1/4*(C(4)+2*C(5)+C(6))
   * which is repaired later due to the "voigt-matrix":
   *    1                 1/2
   *   1/2                1/2
   */

  //  C(0,0)= ScalarThis*C(0,0) + ScalarAB * A(0,0)*B(0,0);
  //  C(0,1)= ScalarThis*C(0,1) + ScalarAB * A(0,0)*B(0,1);
  //  C(0,2)= ScalarThis*C(0,2) + ScalarAB * A(0,0)*B(0,2);
  //  C(0,3)= ScalarThis*C(0,3) + ScalarAB * (A(0,1)*B(0,0) + A(0,2)*B(0,0));
  //  C(0,4)= ScalarThis*C(0,4) + ScalarAB * (A(0,1)*B(0,1) + A(0,2)*B(0,1));
  //  C(0,5)= ScalarThis*C(0,5) + ScalarAB * (A(0,1)*B(0,2) + A(0,2)*B(0,2));
  //
  //  C(1,0)= ScalarThis*C(1,0) + ScalarAB * A(0,0)*B(1,0);
  //  C(1,1)= ScalarThis*C(1,1) + ScalarAB * A(0,0)*B(1,1);
  //  C(1,2)= ScalarThis*C(1,2) + ScalarAB * A(0,0)*B(1,2);
  //  C(1,3)= ScalarThis*C(1,3) + ScalarAB * (A(0,1)*B(1,0) + A(0,2)*B(1,0));
  //  C(1,4)= ScalarThis*C(1,4) + ScalarAB * (A(0,1)*B(1,1) + A(0,2)*B(1,1));
  //  C(1,5)= ScalarThis*C(1,5) + ScalarAB * (A(0,1)*B(1,2) + A(0,2)*B(1,2));
  //
  //  C(2,0)= ScalarThis*C(2,0) + ScalarAB * A(0,0)*B(2,0);
  //  C(2,1)= ScalarThis*C(2,1) + ScalarAB * A(0,0)*B(2,1);
  //  C(2,2)= ScalarThis*C(2,2) + ScalarAB * A(0,0)*B(2,2);
  //  C(2,3)= ScalarThis*C(2,3) + ScalarAB * (A(0,1)*B(2,0) + A(0,2)*B(2,0));
  //  C(2,4)= ScalarThis*C(2,4) + ScalarAB * (A(0,1)*B(2,1) + A(0,2)*B(2,1));
  //  C(2,5)= ScalarThis*C(2,5) + ScalarAB * (A(0,1)*B(2,2) + A(0,2)*B(2,2));
  //
  //  C(3,0)= ScalarThis*C(3,0) + ScalarAB * (A(1,0)*B(0,0) + A(2,0)*B(0,0));
  //  C(3,1)= ScalarThis*C(3,1) + ScalarAB * (A(1,0)*B(0,1) + A(2,0)*B(0,1));
  //  C(3,2)= ScalarThis*C(3,2) + ScalarAB * (A(1,0)*B(0,2) + A(2,0)*B(0,2));
  //  C(3,3)= ScalarThis*C(3,3) + ScalarAB * 0.5*(A(1,1)*B(0,0) + 2.0*A(1,2)*B(0,0) +
  //  A(2,2)*B(0,0)); C(3,4)= ScalarThis*C(3,4) + ScalarAB * 0.5*(A(1,1)*B(0,1) + 2.0*A(1,2)*B(0,1)
  //  + A(2,2)*B(0,1)); C(3,5)= ScalarThis*C(3,5) + ScalarAB * 0.5*(A(1,1)*B(0,2)
  //  + 2.0*A(1,2)*B(0,2) + A(2,2)*B(0,2));
  //
  //  C(4,0)= ScalarThis*C(4,0) + ScalarAB * (A(1,0)*B(1,0) + A(2,0)*B(1,0));
  //  C(4,1)= ScalarThis*C(4,1) + ScalarAB * (A(1,0)*B(1,1) + A(2,0)*B(1,1));
  //  C(4,2)= ScalarThis*C(4,2) + ScalarAB * (A(1,0)*B(1,2) + A(2,0)*B(1,2));
  //  C(4,3)= ScalarThis*C(4,3) + ScalarAB * 0.5*(A(1,1)*B(1,0) + 2.0*A(1,2)*B(1,0) +
  //  A(2,2)*B(1,0)); C(4,4)= ScalarThis*C(4,4) + ScalarAB * 0.5*(A(1,1)*B(1,1) + 2.0*A(1,2)*B(1,1)
  //  + A(2,2)*B(1,1)); C(4,5)= ScalarThis*C(4,5) + ScalarAB * 0.5*(A(1,1)*B(1,2)
  //  + 2.0*A(1,2)*B(1,2) + A(2,2)*B(1,2));
  //
  //  C(5,0)= ScalarThis*C(5,0) + ScalarAB * (A(1,0)*B(2,0) + A(2,0)*B(2,0));
  //  C(5,1)= ScalarThis*C(5,1) + ScalarAB * (A(1,0)*B(2,1) + A(2,0)*B(2,1));
  //  C(5,2)= ScalarThis*C(5,2) + ScalarAB * (A(1,0)*B(2,2) + A(2,0)*B(2,2));
  //  C(5,3)= ScalarThis*C(5,3) + ScalarAB * 0.5*(A(1,1)*B(2,0) + 2.0*A(1,2)*B(2,0) +
  //  A(2,2)*B(2,0)); C(5,4)= ScalarThis*C(5,4) + ScalarAB * 0.5*(A(1,1)*B(2,1) + 2.0*A(1,2)*B(2,1)
  //  + A(2,2)*B(2,1)); C(5,5)= ScalarThis*C(5,5) + ScalarAB * 0.5*(A(1,1)*B(2,2)
  //  + 2.0*A(1,2)*B(2,2) + A(2,2)*B(2,2));

  C(0, 0) =
      ScalarThis * C(0, 0) + ScalarAB * 0.5 * (A(0, 0) * B(0, 0) + A(0, 0) * B(0, 0));  // C1111
  C(0, 1) =
      ScalarThis * C(0, 1) + ScalarAB * 0.5 * (A(0, 1) * B(0, 1) + A(0, 1) * B(0, 1));  // C1122
  C(0, 2) =
      ScalarThis * C(0, 2) + ScalarAB * 0.5 * (A(0, 2) * B(0, 2) + A(0, 2) * B(0, 2));  // C1133
  C(0, 3) =
      ScalarThis * C(0, 3) + ScalarAB * 0.5 * (A(0, 0) * B(0, 1) + A(0, 1) * B(0, 0));  // C1112
  C(0, 4) =
      ScalarThis * C(0, 4) + ScalarAB * 0.5 * (A(0, 1) * B(0, 2) + A(0, 2) * B(0, 1));  // C1123
  C(0, 5) =
      ScalarThis * C(0, 5) + ScalarAB * 0.5 * (A(0, 0) * B(0, 2) + A(0, 2) * B(0, 0));  // C1113

  C(1, 0) =
      ScalarThis * C(1, 0) + ScalarAB * 0.5 * (A(1, 0) * B(1, 0) + A(1, 0) * B(1, 0));  // C2211
  C(1, 1) =
      ScalarThis * C(1, 1) + ScalarAB * 0.5 * (A(1, 1) * B(1, 1) + A(1, 1) * B(1, 1));  // C2222
  C(1, 2) =
      ScalarThis * C(1, 2) + ScalarAB * 0.5 * (A(1, 2) * B(1, 2) + A(1, 2) * B(1, 2));  // C2233
  C(1, 3) =
      ScalarThis * C(1, 3) + ScalarAB * 0.5 * (A(1, 0) * B(1, 1) + A(1, 1) * B(1, 0));  // C2212
  C(1, 4) =
      ScalarThis * C(1, 4) + ScalarAB * 0.5 * (A(1, 1) * B(1, 2) + A(1, 2) * B(1, 1));  // C2223
  C(1, 5) =
      ScalarThis * C(1, 5) + ScalarAB * 0.5 * (A(1, 0) * B(1, 2) + A(1, 2) * B(1, 0));  // C2213

  C(2, 0) =
      ScalarThis * C(2, 0) + ScalarAB * 0.5 * (A(2, 0) * B(2, 0) + A(2, 0) * B(2, 0));  // C3311
  C(2, 1) =
      ScalarThis * C(2, 1) + ScalarAB * 0.5 * (A(2, 1) * B(2, 1) + A(2, 1) * B(2, 1));  // C3322
  C(2, 2) =
      ScalarThis * C(2, 2) + ScalarAB * 0.5 * (A(2, 2) * B(2, 2) + A(2, 2) * B(2, 2));  // C3333
  C(2, 3) =
      ScalarThis * C(2, 3) + ScalarAB * 0.5 * (A(2, 1) * B(2, 1) + A(2, 1) * B(2, 0));  // C3312
  C(2, 4) =
      ScalarThis * C(2, 4) + ScalarAB * 0.5 * (A(2, 1) * B(2, 2) + A(2, 2) * B(2, 1));  // C3323
  C(2, 5) =
      ScalarThis * C(2, 5) + ScalarAB * 0.5 * (A(2, 0) * B(2, 2) + A(2, 2) * B(2, 0));  // C3313

  C(3, 0) =
      ScalarThis * C(3, 0) + ScalarAB * 0.5 * (A(0, 0) * B(1, 0) + A(0, 0) * B(1, 0));  // C1211
  C(3, 1) =
      ScalarThis * C(3, 1) + ScalarAB * 0.5 * (A(0, 1) * B(1, 1) + A(0, 1) * B(1, 1));  // C1222
  C(3, 2) =
      ScalarThis * C(3, 2) + ScalarAB * 0.5 * (A(0, 2) * B(1, 2) + A(0, 2) * B(1, 2));  // C1233
  C(3, 3) =
      ScalarThis * C(3, 3) + ScalarAB * 0.5 * (A(0, 0) * B(1, 1) + A(0, 1) * B(1, 0));  // C1212
  C(3, 4) =
      ScalarThis * C(3, 4) + ScalarAB * 0.5 * (A(0, 1) * B(1, 2) + A(0, 2) * B(1, 1));  // C1223
  C(3, 5) =
      ScalarThis * C(3, 5) + ScalarAB * 0.5 * (A(0, 0) * B(1, 2) + A(0, 2) * B(1, 0));  // C1213

  C(4, 0) =
      ScalarThis * C(4, 0) + ScalarAB * 0.5 * (A(1, 0) * B(2, 0) + A(1, 0) * B(2, 0));  // C2311
  C(4, 1) =
      ScalarThis * C(4, 1) + ScalarAB * 0.5 * (A(1, 1) * B(2, 1) + A(1, 1) * B(2, 1));  // C2322
  C(4, 2) =
      ScalarThis * C(4, 2) + ScalarAB * 0.5 * (A(1, 2) * B(2, 2) + A(1, 2) * B(2, 2));  // C2333
  C(4, 3) =
      ScalarThis * C(4, 3) + ScalarAB * 0.5 * (A(1, 0) * B(2, 1) + A(1, 1) * B(2, 0));  // C2312
  C(4, 4) =
      ScalarThis * C(4, 4) + ScalarAB * 0.5 * (A(1, 1) * B(2, 2) + A(1, 2) * B(2, 1));  // C2323
  C(4, 5) =
      ScalarThis * C(4, 5) + ScalarAB * 0.5 * (A(1, 0) * B(2, 2) + A(1, 2) * B(2, 0));  // C2313

  C(5, 0) =
      ScalarThis * C(5, 0) + ScalarAB * 0.5 * (A(0, 0) * B(2, 0) + A(0, 0) * B(2, 0));  // C1311
  C(5, 1) =
      ScalarThis * C(5, 1) + ScalarAB * 0.5 * (A(0, 1) * B(2, 1) + A(0, 1) * B(2, 1));  // C1322
  C(5, 2) =
      ScalarThis * C(5, 2) + ScalarAB * 0.5 * (A(0, 2) * B(2, 2) + A(0, 2) * B(2, 2));  // C1333
  C(5, 3) =
      ScalarThis * C(5, 3) + ScalarAB * 0.5 * (A(0, 0) * B(2, 1) + A(0, 1) * B(2, 0));  // C1312
  C(5, 4) =
      ScalarThis * C(5, 4) + ScalarAB * 0.5 * (A(0, 1) * B(2, 2) + A(0, 2) * B(2, 1));  // C1323
  C(5, 5) =
      ScalarThis * C(5, 5) + ScalarAB * 0.5 * (A(0, 0) * B(2, 2) + A(0, 2) * B(2, 0));  // C1313

  return;
}

/*----------------------------------------------------------------------*
 | compute the "material tensor product" A o B (also known as           |
 | kronecker-tensor-product) of two 2nd order tensors                   |
 | (in matrix notation) and add the result to a 4th order tensor        |
 | (also in matrix notation) using the symmetry-conditions inherent to  |
 | material tensors, or tangent matrices, respectively                  |
 | AND the Voigt notation of E,S, and C with the famous factor 2!       |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
// This is a copy of the above function, using the fixed size matrix.ss
void MAT::ElastSymTensor_o_Multiply(LINALG::Matrix<6, 6>& C, const double ScalarAB,
    const LINALG::Matrix<3, 3>& A, const LINALG::Matrix<3, 3>& B, const double ScalarThis)
{
  // To keep the code somewhat shorter I removed the explanation comment here,
  // but they are still in the Epetra Version
  const double ScalarABhalf = ScalarAB * 0.5;
  C(0, 0) = ScalarThis * C(0, 0) + ScalarAB * (A(0, 0) * B(0, 0));                          // C1111
  C(0, 1) = ScalarThis * C(0, 1) + ScalarAB * (A(0, 1) * B(0, 1));                          // C1122
  C(0, 2) = ScalarThis * C(0, 2) + ScalarAB * (A(0, 2) * B(0, 2));                          // C1133
  C(0, 3) = ScalarThis * C(0, 3) + ScalarABhalf * (A(0, 0) * B(0, 1) + A(0, 1) * B(0, 0));  // C1112
  C(0, 4) = ScalarThis * C(0, 4) + ScalarABhalf * (A(0, 1) * B(0, 2) + A(0, 2) * B(0, 1));  // C1123
  C(0, 5) = ScalarThis * C(0, 5) + ScalarABhalf * (A(0, 0) * B(0, 2) + A(0, 2) * B(0, 0));  // C1113

  C(1, 0) = ScalarThis * C(1, 0) + ScalarAB * (A(1, 0) * B(1, 0));                          // C2211
  C(1, 1) = ScalarThis * C(1, 1) + ScalarAB * (A(1, 1) * B(1, 1));                          // C2222
  C(1, 2) = ScalarThis * C(1, 2) + ScalarAB * (A(1, 2) * B(1, 2));                          // C2233
  C(1, 3) = ScalarThis * C(1, 3) + ScalarABhalf * (A(1, 0) * B(1, 1) + A(1, 1) * B(1, 0));  // C2212
  C(1, 4) = ScalarThis * C(1, 4) + ScalarABhalf * (A(1, 1) * B(1, 2) + A(1, 2) * B(1, 1));  // C2223
  C(1, 5) = ScalarThis * C(1, 5) + ScalarABhalf * (A(1, 0) * B(1, 2) + A(1, 2) * B(1, 0));  // C2213

  C(2, 0) = ScalarThis * C(2, 0) + ScalarAB * (A(2, 0) * B(2, 0));                          // C3311
  C(2, 1) = ScalarThis * C(2, 1) + ScalarAB * (A(2, 1) * B(2, 1));                          // C3322
  C(2, 2) = ScalarThis * C(2, 2) + ScalarAB * (A(2, 2) * B(2, 2));                          // C3333
  C(2, 3) = ScalarThis * C(2, 3) + ScalarABhalf * (A(2, 1) * B(2, 1) + A(2, 1) * B(2, 0));  // C3312
  C(2, 4) = ScalarThis * C(2, 4) + ScalarABhalf * (A(2, 1) * B(2, 2) + A(2, 2) * B(2, 1));  // C3323
  C(2, 5) = ScalarThis * C(2, 5) + ScalarABhalf * (A(2, 0) * B(2, 2) + A(2, 2) * B(2, 0));  // C3313

  C(3, 0) = ScalarThis * C(3, 0) + ScalarAB * (A(0, 0) * B(1, 0));                          // C1211
  C(3, 1) = ScalarThis * C(3, 1) + ScalarAB * (A(0, 1) * B(1, 1));                          // C1222
  C(3, 2) = ScalarThis * C(3, 2) + ScalarAB * (A(0, 2) * B(1, 2));                          // C1233
  C(3, 3) = ScalarThis * C(3, 3) + ScalarABhalf * (A(0, 0) * B(1, 1) + A(0, 1) * B(1, 0));  // C1212
  C(3, 4) = ScalarThis * C(3, 4) + ScalarABhalf * (A(0, 1) * B(1, 2) + A(0, 2) * B(1, 1));  // C1223
  C(3, 5) = ScalarThis * C(3, 5) + ScalarABhalf * (A(0, 0) * B(1, 2) + A(0, 2) * B(1, 0));  // C1213

  C(4, 0) = ScalarThis * C(4, 0) + ScalarAB * (A(1, 0) * B(2, 0));                          // C2311
  C(4, 1) = ScalarThis * C(4, 1) + ScalarAB * (A(1, 1) * B(2, 1));                          // C2322
  C(4, 2) = ScalarThis * C(4, 2) + ScalarAB * (A(1, 2) * B(2, 2));                          // C2333
  C(4, 3) = ScalarThis * C(4, 3) + ScalarABhalf * (A(1, 0) * B(2, 1) + A(1, 1) * B(2, 0));  // C2312
  C(4, 4) = ScalarThis * C(4, 4) + ScalarABhalf * (A(1, 1) * B(2, 2) + A(1, 2) * B(2, 1));  // C2323
  C(4, 5) = ScalarThis * C(4, 5) + ScalarABhalf * (A(1, 0) * B(2, 2) + A(1, 2) * B(2, 0));  // C2313

  C(5, 0) = ScalarThis * C(5, 0) + ScalarAB * (A(0, 0) * B(2, 0));                          // C1311
  C(5, 1) = ScalarThis * C(5, 1) + ScalarAB * (A(0, 1) * B(2, 1));                          // C1322
  C(5, 2) = ScalarThis * C(5, 2) + ScalarAB * (A(0, 2) * B(2, 2));                          // C1333
  C(5, 3) = ScalarThis * C(5, 3) + ScalarABhalf * (A(0, 0) * B(2, 1) + A(0, 1) * B(2, 0));  // C1312
  C(5, 4) = ScalarThis * C(5, 4) + ScalarABhalf * (A(0, 1) * B(2, 2) + A(0, 2) * B(2, 1));  // C1323
  C(5, 5) = ScalarThis * C(5, 5) + ScalarABhalf * (A(0, 0) * B(2, 2) + A(0, 2) * B(2, 0));  // C1313

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::VolumetrifyAndIsochorify(LINALG::Matrix<6, 1>* pk2vol, LINALG::Matrix<6, 6>* cvol,
    LINALG::Matrix<6, 1>* pk2iso, LINALG::Matrix<6, 6>* ciso, const LINALG::Matrix<6, 1>& gl,
    const LINALG::Matrix<6, 1>& pk2, const LINALG::Matrix<6, 6>& cmat)
{
  // useful call?
#ifdef DEBUG
  if ((pk2vol == NULL) and (cvol == NULL) and (pk2iso == NULL) and (ciso == NULL))
    dserror("Useful call? Apparently you do not want to compute anything");
#endif

  // right Cauchy--Green tensor
  // REMARK: stored in _strain_-like 6-Voigt vector
  LINALG::Matrix<6, 1> rcg(gl);
  rcg.Scale(2.0);
  for (int i = 0; i < 3; i++) rcg(i) += 1.0;

  // third invariant (determinant) of right Cauchy--Green strains
  const double rcg3rd = rcg(0) * rcg(1) * rcg(2) + 0.25 * rcg(3) * rcg(4) * rcg(5) -
                        0.25 * rcg(1) * rcg(5) * rcg(5) - 0.25 * rcg(2) * rcg(3) * rcg(3) -
                        0.25 * rcg(0) * rcg(4) * rcg(4);

  // inverse right Cauchy--Green tensor C^{-1}
  // REMARK: stored in as _stress_ 6-Voigt vector
  LINALG::Matrix<6, 1> icg(false);
  icg(0) = (rcg(1) * rcg(2) - 0.25 * rcg(4) * rcg(4)) / rcg3rd;        // (C^{-1})^{11}
  icg(1) = (rcg(0) * rcg(2) - 0.25 * rcg(5) * rcg(5)) / rcg3rd;        // (C^{-1})^{22}
  icg(2) = (rcg(0) * rcg(1) - 0.25 * rcg(3) * rcg(3)) / rcg3rd;        // (C^{-1})^{33}
  icg(3) = (0.25 * rcg(5) * rcg(4) - 0.5 * rcg(3) * rcg(2)) / rcg3rd;  // (C^{-1})^{12}
  icg(4) = (0.25 * rcg(3) * rcg(5) - 0.5 * rcg(0) * rcg(4)) / rcg3rd;  // (C^{-1})^{23}
  icg(5) = (0.25 * rcg(3) * rcg(4) - 0.5 * rcg(5) * rcg(1)) / rcg3rd;  // (C^{-1})^{31}

  // double contraction of 2nd Piola--Kirchhoff stress and right Cauchy--Green strain,
  // i.e. in index notation S^{AB} C_{AB}
  // REMARK: equal to S^T C, because S is stress-like and C is strain-like 6-Voigt vector
  const double pk2rcg = pk2.Dot(rcg);

  // stress splitting
  {
    // volumetric 2nd Piola--Kirchhoff stress
    LINALG::Matrix<6, 1> pk2vol_(false);
    if (pk2vol != NULL) pk2vol_.SetView(*pk2vol);
    pk2vol_.Update(pk2rcg / 3.0, icg);

    // isochoric 2nd Piola--Kirchhoff stress
    // S^{AB}_iso = S^{AB} - S^{AB}_{vol}
    if (pk2iso != NULL) pk2iso->Update(1.0, pk2, -1.0, pk2vol_);
  }

  // elasticity tensor splitting
  {
    // 'linearised' 2nd Piola--Kirchhoff stress
    // S^{CD}_lin = S^{CD} + 1/2 C_{AB} C^{ABCD}
    LINALG::Matrix<6, 1> pk2lin(pk2);
    pk2lin.MultiplyTN(0.5, cmat, rcg, 1.0);  // transpose on purpose

    // volumetric part of constitutive tensor
    // C^{ABCD}_vol = 2/3 (C^{-1})^{AB} S^{CD}_lin
    //              - 2/3 (S^{EF} C_{EF}) ( 1/2 (
    //                (C^{-1})^{AC} (C^{-1})^{BD} + (C^{-1})^{AD} (C^{-1})^{BC}
    //              ) )
    LINALG::Matrix<6, 6> cvol_(false);
    if (cvol != NULL) cvol_.SetView(*cvol);
    cvol_.MultiplyNT(2.0 / 3.0, icg, pk2lin);
    AddtoCmatHolzapfelProduct(cvol_, icg, -2.0 / 3.0 * pk2rcg);

    // isochoric part of constitutive tensor
    // C^{ABCD}_iso = C^{ABCD} - C^{ABCD}_vol
    if (ciso != NULL) ciso->Update(1.0, cmat, -1.0, cvol_);
  }

  //
  return;
}

/*----------------------------------------------------------------------*
 |  add dX^2/dX to cmat                                     seitz 07/13 |
 *----------------------------------------------------------------------*/
void MAT::AddToCmatDerivTensorSquare(
    LINALG::Matrix<6, 6>& C, double ScalarDX2, LINALG::Matrix<3, 3> X, double ScalarThis)
{
  C(0, 0) = ScalarThis * C(0, 0) + ScalarDX2 * 2. * X(0, 0);  // C1111
  C(0, 1) = ScalarThis * C(0, 1);                             // C1122
  C(0, 2) = ScalarThis * C(0, 2);                             // C1133
  C(0, 3) = ScalarThis * C(0, 3) + ScalarDX2 * X(0, 1);       // C1112
  C(0, 4) = ScalarThis * C(0, 4);                             // C1123
  C(0, 5) = ScalarThis * C(0, 5) + ScalarDX2 * X(0, 2);       // C1113

  C(1, 0) = ScalarThis * C(1, 0);                             // C2211
  C(1, 1) = ScalarThis * C(1, 1) + ScalarDX2 * 2. * X(1, 1);  // C2222
  C(1, 2) = ScalarThis * C(1, 2);                             // C2233
  C(1, 3) = ScalarThis * C(1, 3) + ScalarDX2 * X(0, 1);       // C2212
  C(1, 4) = ScalarThis * C(1, 4) + ScalarDX2 * X(1, 2);       // C2223
  C(1, 5) = ScalarThis * C(1, 5);                             // C2213

  C(2, 0) = ScalarThis * C(2, 0);                             // C3311
  C(2, 1) = ScalarThis * C(2, 1);                             // C3322
  C(2, 2) = ScalarThis * C(2, 2) + ScalarDX2 * 2. * X(2, 2);  // C3333
  C(2, 3) = ScalarThis * C(2, 3);                             // C3312
  C(2, 4) = ScalarThis * C(2, 4) + ScalarDX2 * X(1, 2);       // C3323
  C(2, 5) = ScalarThis * C(2, 5) + ScalarDX2 * X(0, 2);       // C3313

  C(3, 0) = ScalarThis * C(3, 0) + ScalarDX2 * X(0, 1);                    // C1211
  C(3, 1) = ScalarThis * C(3, 1) + ScalarDX2 * X(0, 1);                    // C1222
  C(3, 2) = ScalarThis * C(3, 2);                                          // C1233
  C(3, 3) = ScalarThis * C(3, 3) + ScalarDX2 * 0.5 * (X(0, 0) + X(1, 1));  // C1212
  C(3, 4) = ScalarThis * C(3, 4) + ScalarDX2 * 0.5 * X(0, 2);              // C1223
  C(3, 5) = ScalarThis * C(3, 5) + ScalarDX2 * 0.5 * X(1, 2);              // C1213

  C(4, 0) = ScalarThis * C(4, 0);                                          // C2311
  C(4, 1) = ScalarThis * C(4, 1) + ScalarDX2 * X(1, 2);                    // C2322
  C(4, 2) = ScalarThis * C(4, 2) + ScalarDX2 * X(1, 2);                    // C2333
  C(4, 3) = ScalarThis * C(4, 3) + ScalarDX2 * 0.5 * X(0, 2);              // C2312
  C(4, 4) = ScalarThis * C(4, 4) + ScalarDX2 * 0.5 * (X(1, 1) + X(2, 2));  // C2323
  C(4, 5) = ScalarThis * C(4, 5) + ScalarDX2 * 0.5 * X(0, 1);              // C2313

  C(5, 0) = ScalarThis * C(5, 0) + ScalarDX2 * X(0, 2);                    // C1311
  C(5, 1) = ScalarThis * C(5, 1);                                          // C1322
  C(5, 2) = ScalarThis * C(5, 2) + ScalarDX2 * X(0, 2);                    // C1333
  C(5, 3) = ScalarThis * C(5, 3) + ScalarDX2 * 0.5 * X(1, 2);              // C1312
  C(5, 4) = ScalarThis * C(5, 4) + ScalarDX2 * 0.5 * X(0, 1);              // C1323
  C(5, 5) = ScalarThis * C(5, 5) + ScalarDX2 * 0.5 * (X(2, 2) + X(0, 0));  // C1313

  return;
}


/*----------------------------------------------------------------------*
 |                                                          seitz 09/13 |
 *----------------------------------------------------------------------*/
void MAT::AddSymmetricHolzapfelProduct(LINALG::Matrix<6, 6>& X, const LINALG::Matrix<3, 3>& A,
    const LINALG::Matrix<3, 3>& B, const double fac)
{
  X(0, 0) += 4 * fac * A(0, 0) * B(0, 0);
  X(0, 3) += fac * (2 * A(0, 0) * B(1, 0) + 2 * A(1, 0) * B(0, 0));
  X(0, 5) += fac * (2 * A(0, 0) * B(2, 0) + 2 * A(2, 0) * B(0, 0));
  X(0, 1) += 4 * fac * A(1, 0) * B(1, 0);
  X(0, 4) += fac * (2 * A(1, 0) * B(2, 0) + 2 * A(2, 0) * B(1, 0));
  X(0, 2) += 4 * fac * A(2, 0) * B(2, 0);

  X(3, 0) += fac * (2 * A(0, 0) * B(0, 1) + 2 * A(0, 1) * B(0, 0));
  X(3, 3) += fac * (A(0, 0) * B(1, 1) + A(1, 0) * B(0, 1) + A(1, 1) * B(0, 0) + A(0, 1) * B(1, 0));
  X(3, 5) += fac * (A(0, 0) * B(2, 1) + A(2, 0) * B(0, 1) + A(2, 1) * B(0, 0) + A(0, 1) * B(2, 0));
  X(3, 1) += fac * (2 * A(1, 0) * B(1, 1) + 2 * A(1, 1) * B(1, 0));
  X(3, 4) += fac * (A(1, 0) * B(2, 1) + A(2, 0) * B(1, 1) + A(2, 1) * B(1, 0) + A(1, 1) * B(2, 0));
  X(3, 2) += fac * (2 * A(2, 0) * B(2, 1) + 2 * A(2, 1) * B(2, 0));

  X(5, 0) += fac * (2 * A(0, 0) * B(0, 2) + 2 * A(0, 2) * B(0, 0));
  X(5, 3) += fac * (A(0, 0) * B(1, 2) + A(1, 0) * B(0, 2) + A(1, 2) * B(0, 0) + A(0, 2) * B(1, 0));
  X(5, 5) += fac * (A(0, 0) * B(2, 2) + A(2, 0) * B(0, 2) + A(2, 2) * B(0, 0) + A(0, 2) * B(2, 0));
  X(5, 1) += fac * (2 * A(1, 0) * B(1, 2) + 2 * A(1, 2) * B(1, 0));
  X(5, 4) += fac * (A(1, 0) * B(2, 2) + A(2, 0) * B(1, 2) + A(2, 2) * B(1, 0) + A(1, 2) * B(2, 0));
  X(5, 2) += fac * (2 * A(2, 0) * B(2, 2) + 2 * A(2, 2) * B(2, 0));

  X(1, 0) += 4 * fac * A(0, 1) * B(0, 1);
  X(1, 3) += fac * (2 * A(0, 1) * B(1, 1) + 2 * A(1, 1) * B(0, 1));
  X(1, 5) += fac * (2 * A(0, 1) * B(2, 1) + 2 * A(2, 1) * B(0, 1));
  X(1, 1) += 4 * fac * A(1, 1) * B(1, 1);
  X(1, 4) += fac * (2 * A(1, 1) * B(2, 1) + 2 * A(2, 1) * B(1, 1));
  X(1, 2) += 4 * fac * A(2, 1) * B(2, 1);

  X(4, 0) += fac * (2 * A(0, 1) * B(0, 2) + 2 * A(0, 2) * B(0, 1));
  X(4, 3) += fac * (A(0, 1) * B(1, 2) + A(1, 1) * B(0, 2) + A(1, 2) * B(0, 1) + A(0, 2) * B(1, 1));
  X(4, 5) += fac * (A(0, 1) * B(2, 2) + A(2, 1) * B(0, 2) + A(2, 2) * B(0, 1) + A(0, 2) * B(2, 1));
  X(4, 1) += fac * (2 * A(1, 1) * B(1, 2) + 2 * A(1, 2) * B(1, 1));
  X(4, 4) += fac * (A(1, 1) * B(2, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(1, 1) + A(1, 2) * B(2, 1));
  X(4, 2) += fac * (2 * A(2, 1) * B(2, 2) + 2 * A(2, 2) * B(2, 1));

  X(2, 0) += 4 * fac * A(0, 2) * B(0, 2);
  X(2, 3) += fac * (2 * A(0, 2) * B(1, 2) + 2 * A(1, 2) * B(0, 2));
  X(2, 5) += fac * (2 * A(0, 2) * B(2, 2) + 2 * A(2, 2) * B(0, 2));
  X(2, 1) += 4 * fac * A(1, 2) * B(1, 2);
  X(2, 4) += fac * (2 * A(1, 2) * B(2, 2) + 2 * A(2, 2) * B(1, 2));
  X(2, 2) += 4 * fac * A(2, 2) * B(2, 2);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
void MAT::AddRightNonSymmetricHolzapfelProduct(LINALG::Matrix<6, 9, T>& out,
    LINALG::Matrix<3, 3, T> const& A, LINALG::Matrix<3, 3, T> const& B, T const fac)
{
  out(0, 0) += 2 * fac * A(0, 0) * B(0, 0);
  out(0, 3) += 2 * fac * A(0, 0) * B(0, 1);
  out(0, 5) += 2 * fac * A(0, 0) * B(0, 2);
  out(0, 6) += 2 * fac * A(0, 1) * B(0, 0);
  out(0, 1) += 2 * fac * A(0, 1) * B(0, 1);
  out(0, 4) += 2 * fac * A(0, 1) * B(0, 2);
  out(0, 8) += 2 * fac * A(0, 2) * B(0, 0);
  out(0, 7) += 2 * fac * A(0, 2) * B(0, 1);
  out(0, 2) += 2 * fac * A(0, 2) * B(0, 2);

  out(1, 0) += 2 * fac * A(1, 0) * B(1, 0);
  out(1, 3) += 2 * fac * A(1, 0) * B(1, 1);
  out(1, 5) += 2 * fac * A(1, 0) * B(1, 2);
  out(1, 6) += 2 * fac * A(1, 1) * B(1, 0);
  out(1, 1) += 2 * fac * A(1, 1) * B(1, 1);
  out(1, 4) += 2 * fac * A(1, 1) * B(1, 2);
  out(1, 8) += 2 * fac * A(1, 2) * B(1, 0);
  out(1, 7) += 2 * fac * A(1, 2) * B(1, 1);
  out(1, 2) += 2 * fac * A(1, 2) * B(1, 2);

  out(2, 0) += 2 * fac * A(2, 0) * B(2, 0);
  out(2, 3) += 2 * fac * A(2, 0) * B(2, 1);
  out(2, 5) += 2 * fac * A(2, 0) * B(2, 2);
  out(2, 6) += 2 * fac * A(2, 1) * B(2, 0);
  out(2, 1) += 2 * fac * A(2, 1) * B(2, 1);
  out(2, 4) += 2 * fac * A(2, 1) * B(2, 2);
  out(2, 8) += 2 * fac * A(2, 2) * B(2, 0);
  out(2, 7) += 2 * fac * A(2, 2) * B(2, 1);
  out(2, 2) += 2 * fac * A(2, 2) * B(2, 2);

  out(3, 0) += fac * (A(0, 0) * B(1, 0) + A(1, 0) * B(0, 0));
  out(3, 3) += fac * (A(0, 0) * B(1, 1) + A(1, 0) * B(0, 1));
  out(3, 5) += fac * (A(0, 0) * B(1, 2) + A(1, 0) * B(0, 2));
  out(3, 6) += fac * (A(0, 1) * B(1, 0) + A(1, 1) * B(0, 0));
  out(3, 1) += fac * (A(0, 1) * B(1, 1) + A(1, 1) * B(0, 1));
  out(3, 4) += fac * (A(0, 1) * B(1, 2) + A(1, 1) * B(0, 2));
  out(3, 8) += fac * (A(0, 2) * B(1, 0) + A(1, 2) * B(0, 0));
  out(3, 7) += fac * (A(0, 2) * B(1, 1) + A(1, 2) * B(0, 1));
  out(3, 2) += fac * (A(0, 2) * B(1, 2) + A(1, 2) * B(0, 2));

  out(4, 0) += fac * (A(1, 0) * B(2, 0) + A(2, 0) * B(1, 0));
  out(4, 3) += fac * (A(1, 0) * B(2, 1) + A(2, 0) * B(1, 1));
  out(4, 5) += fac * (A(1, 0) * B(2, 2) + A(2, 0) * B(1, 2));
  out(4, 6) += fac * (A(1, 1) * B(2, 0) + A(2, 1) * B(1, 0));
  out(4, 1) += fac * (A(1, 1) * B(2, 1) + A(2, 1) * B(1, 1));
  out(4, 4) += fac * (A(1, 1) * B(2, 2) + A(2, 1) * B(1, 2));
  out(4, 8) += fac * (A(1, 2) * B(2, 0) + A(2, 2) * B(1, 0));
  out(4, 7) += fac * (A(1, 2) * B(2, 1) + A(2, 2) * B(1, 1));
  out(4, 2) += fac * (A(1, 2) * B(2, 2) + A(2, 2) * B(1, 2));

  out(5, 0) += fac * (A(0, 0) * B(2, 0) + A(2, 0) * B(0, 0));
  out(5, 3) += fac * (A(0, 0) * B(2, 1) + A(2, 0) * B(0, 1));
  out(5, 5) += fac * (A(0, 0) * B(2, 2) + A(2, 0) * B(0, 2));
  out(5, 6) += fac * (A(0, 1) * B(2, 0) + A(2, 1) * B(0, 0));
  out(5, 1) += fac * (A(0, 1) * B(2, 1) + A(2, 1) * B(0, 1));
  out(5, 4) += fac * (A(0, 1) * B(2, 2) + A(2, 1) * B(0, 2));
  out(5, 8) += fac * (A(0, 2) * B(2, 0) + A(2, 2) * B(0, 0));
  out(5, 7) += fac * (A(0, 2) * B(2, 1) + A(2, 2) * B(0, 1));
  out(5, 2) += fac * (A(0, 2) * B(2, 2) + A(2, 2) * B(0, 2));

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
void MAT::AddRightNonSymmetricHolzapfelProductStrainLike(LINALG::Matrix<6, 9, T>& out,
    LINALG::Matrix<3, 3, T> const& A, LINALG::Matrix<3, 3, T> const& B, T const fac)
{
  out(0, 0) += 2 * fac * A(0, 0) * B(0, 0);
  out(0, 3) += 2 * fac * A(0, 0) * B(0, 1);
  out(0, 5) += 2 * fac * A(0, 0) * B(0, 2);
  out(0, 6) += 2 * fac * A(0, 1) * B(0, 0);
  out(0, 1) += 2 * fac * A(0, 1) * B(0, 1);
  out(0, 4) += 2 * fac * A(0, 1) * B(0, 2);
  out(0, 8) += 2 * fac * A(0, 2) * B(0, 0);
  out(0, 7) += 2 * fac * A(0, 2) * B(0, 1);
  out(0, 2) += 2 * fac * A(0, 2) * B(0, 2);

  out(1, 0) += 2 * fac * A(1, 0) * B(1, 0);
  out(1, 3) += 2 * fac * A(1, 0) * B(1, 1);
  out(1, 5) += 2 * fac * A(1, 0) * B(1, 2);
  out(1, 6) += 2 * fac * A(1, 1) * B(1, 0);
  out(1, 1) += 2 * fac * A(1, 1) * B(1, 1);
  out(1, 4) += 2 * fac * A(1, 1) * B(1, 2);
  out(1, 8) += 2 * fac * A(1, 2) * B(1, 0);
  out(1, 7) += 2 * fac * A(1, 2) * B(1, 1);
  out(1, 2) += 2 * fac * A(1, 2) * B(1, 2);

  out(2, 0) += 2 * fac * A(2, 0) * B(2, 0);
  out(2, 3) += 2 * fac * A(2, 0) * B(2, 1);
  out(2, 5) += 2 * fac * A(2, 0) * B(2, 2);
  out(2, 6) += 2 * fac * A(2, 1) * B(2, 0);
  out(2, 1) += 2 * fac * A(2, 1) * B(2, 1);
  out(2, 4) += 2 * fac * A(2, 1) * B(2, 2);
  out(2, 8) += 2 * fac * A(2, 2) * B(2, 0);
  out(2, 7) += 2 * fac * A(2, 2) * B(2, 1);
  out(2, 2) += 2 * fac * A(2, 2) * B(2, 2);

  out(3, 0) += 2 * fac * (A(0, 0) * B(1, 0) + A(1, 0) * B(0, 0));
  out(3, 3) += 2 * fac * (A(0, 0) * B(1, 1) + A(1, 0) * B(0, 1));
  out(3, 5) += 2 * fac * (A(0, 0) * B(1, 2) + A(1, 0) * B(0, 2));
  out(3, 6) += 2 * fac * (A(0, 1) * B(1, 0) + A(1, 1) * B(0, 0));
  out(3, 1) += 2 * fac * (A(0, 1) * B(1, 1) + A(1, 1) * B(0, 1));
  out(3, 4) += 2 * fac * (A(0, 1) * B(1, 2) + A(1, 1) * B(0, 2));
  out(3, 8) += 2 * fac * (A(0, 2) * B(1, 0) + A(1, 2) * B(0, 0));
  out(3, 7) += 2 * fac * (A(0, 2) * B(1, 1) + A(1, 2) * B(0, 1));
  out(3, 2) += 2 * fac * (A(0, 2) * B(1, 2) + A(1, 2) * B(0, 2));

  out(4, 0) += 2 * fac * (A(1, 0) * B(2, 0) + A(2, 0) * B(1, 0));
  out(4, 3) += 2 * fac * (A(1, 0) * B(2, 1) + A(2, 0) * B(1, 1));
  out(4, 5) += 2 * fac * (A(1, 0) * B(2, 2) + A(2, 0) * B(1, 2));
  out(4, 6) += 2 * fac * (A(1, 1) * B(2, 0) + A(2, 1) * B(1, 0));
  out(4, 1) += 2 * fac * (A(1, 1) * B(2, 1) + A(2, 1) * B(1, 1));
  out(4, 4) += 2 * fac * (A(1, 1) * B(2, 2) + A(2, 1) * B(1, 2));
  out(4, 8) += 2 * fac * (A(1, 2) * B(2, 0) + A(2, 2) * B(1, 0));
  out(4, 7) += 2 * fac * (A(1, 2) * B(2, 1) + A(2, 2) * B(1, 1));
  out(4, 2) += 2 * fac * (A(1, 2) * B(2, 2) + A(2, 2) * B(1, 2));

  out(5, 0) += 2 * fac * (A(0, 0) * B(2, 0) + A(2, 0) * B(0, 0));
  out(5, 3) += 2 * fac * (A(0, 0) * B(2, 1) + A(2, 0) * B(0, 1));
  out(5, 5) += 2 * fac * (A(0, 0) * B(2, 2) + A(2, 0) * B(0, 2));
  out(5, 6) += 2 * fac * (A(0, 1) * B(2, 0) + A(2, 1) * B(0, 0));
  out(5, 1) += 2 * fac * (A(0, 1) * B(2, 1) + A(2, 1) * B(0, 1));
  out(5, 4) += 2 * fac * (A(0, 1) * B(2, 2) + A(2, 1) * B(0, 2));
  out(5, 8) += 2 * fac * (A(0, 2) * B(2, 0) + A(2, 2) * B(0, 0));
  out(5, 7) += 2 * fac * (A(0, 2) * B(2, 1) + A(2, 2) * B(0, 1));
  out(5, 2) += 2 * fac * (A(0, 2) * B(2, 2) + A(2, 2) * B(0, 2));

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::AddLeftNonSymmetricHolzapfelProduct(LINALG::Matrix<9, 6>& out,
    LINALG::Matrix<3, 3> const& A, LINALG::Matrix<3, 3> const& B, double const fac)
{
  out(0, 0) += fac * 2 * A(0, 0) * B(0, 0);
  out(0, 3) += fac * (A(0, 0) * B(0, 1) + A(0, 1) * B(0, 0));
  out(0, 5) += fac * (A(0, 0) * B(0, 2) + A(0, 2) * B(0, 0));
  out(0, 1) += fac * 2 * A(0, 1) * B(0, 1);
  out(0, 4) += fac * (A(0, 1) * B(0, 2) + A(0, 2) * B(0, 1));
  out(0, 2) += fac * 2 * A(0, 2) * B(0, 2);

  out(3, 0) += fac * 2 * A(0, 0) * B(1, 0);
  out(3, 3) += fac * (A(0, 0) * B(1, 1) + A(0, 1) * B(1, 0));
  out(3, 5) += fac * (A(0, 0) * B(1, 2) + A(0, 2) * B(1, 0));
  out(3, 1) += fac * 2 * A(0, 1) * B(1, 1);
  out(3, 4) += fac * (A(0, 1) * B(1, 2) + A(0, 2) * B(1, 1));
  out(3, 2) += fac * 2 * A(0, 2) * B(1, 2);

  out(5, 0) += fac * 2 * A(0, 0) * B(2, 0);
  out(5, 3) += fac * (A(0, 0) * B(2, 1) + A(0, 1) * B(2, 0));
  out(5, 5) += fac * (A(0, 0) * B(2, 2) + A(0, 2) * B(2, 0));
  out(5, 1) += fac * 2 * A(0, 1) * B(2, 1);
  out(5, 4) += fac * (A(0, 1) * B(2, 2) + A(0, 2) * B(2, 1));
  out(5, 2) += fac * 2 * A(0, 2) * B(2, 2);

  out(6, 0) += fac * 2 * A(1, 0) * B(0, 0);
  out(6, 3) += fac * (A(1, 0) * B(0, 1) + A(1, 1) * B(0, 0));
  out(6, 5) += fac * (A(1, 0) * B(0, 2) + A(1, 2) * B(0, 0));
  out(6, 1) += fac * 2 * A(1, 1) * B(0, 1);
  out(6, 4) += fac * (A(1, 1) * B(0, 2) + A(1, 2) * B(0, 1));
  out(6, 2) += fac * 2 * A(1, 2) * B(0, 2);

  out(1, 0) += fac * 2 * A(1, 0) * B(1, 0);
  out(1, 3) += fac * (A(1, 0) * B(1, 1) + A(1, 1) * B(1, 0));
  out(1, 5) += fac * (A(1, 0) * B(1, 2) + A(1, 2) * B(1, 0));
  out(1, 1) += fac * 2 * A(1, 1) * B(1, 1);
  out(1, 4) += fac * (A(1, 1) * B(1, 2) + A(1, 2) * B(1, 1));
  out(1, 2) += fac * 2 * A(1, 2) * B(1, 2);

  out(4, 0) += fac * 2 * A(1, 0) * B(2, 0);
  out(4, 3) += fac * (A(1, 0) * B(2, 1) + A(1, 1) * B(2, 0));
  out(4, 5) += fac * (A(1, 0) * B(2, 2) + A(1, 2) * B(2, 0));
  out(4, 1) += fac * 2 * A(1, 1) * B(2, 1);
  out(4, 4) += fac * (A(1, 1) * B(2, 2) + A(1, 2) * B(2, 1));
  out(4, 2) += fac * 2 * A(1, 2) * B(2, 2);

  out(8, 0) += fac * 2 * A(2, 0) * B(0, 0);
  out(8, 3) += fac * (A(2, 0) * B(0, 1) + A(2, 1) * B(0, 0));
  out(8, 5) += fac * (A(2, 0) * B(0, 2) + A(2, 2) * B(0, 0));
  out(8, 1) += fac * 2 * A(2, 1) * B(0, 1);
  out(8, 4) += fac * (A(2, 1) * B(0, 2) + A(2, 2) * B(0, 1));
  out(8, 2) += fac * 2 * A(2, 2) * B(0, 2);

  out(7, 0) += fac * 2 * A(2, 0) * B(1, 0);
  out(7, 3) += fac * (A(2, 0) * B(1, 1) + A(2, 1) * B(1, 0));
  out(7, 5) += fac * (A(2, 0) * B(1, 2) + A(2, 2) * B(1, 0));
  out(7, 1) += fac * 2 * A(2, 1) * B(1, 1);
  out(7, 4) += fac * (A(2, 1) * B(1, 2) + A(2, 2) * B(1, 1));
  out(7, 2) += fac * 2 * A(2, 2) * B(1, 2);

  out(2, 0) += fac * 2 * A(2, 0) * B(2, 0);
  out(2, 3) += fac * (A(2, 0) * B(2, 1) + A(2, 1) * B(2, 0));
  out(2, 5) += fac * (A(2, 0) * B(2, 2) + A(2, 2) * B(2, 0));
  out(2, 1) += fac * 2 * A(2, 1) * B(2, 1);
  out(2, 4) += fac * (A(2, 1) * B(2, 2) + A(2, 2) * B(2, 1));
  out(2, 2) += fac * 2 * A(2, 2) * B(2, 2);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::AddNonSymmetricProduct(double const& fac, LINALG::Matrix<3, 3> const& A,
    LINALG::Matrix<3, 3> const& B, LINALG::Matrix<9, 9>& out)
{
  out(0, 0) += fac * A(0, 0) * B(0, 0);
  out(0, 3) += fac * A(0, 0) * B(1, 0);
  out(0, 5) += fac * A(0, 0) * B(2, 0);
  out(0, 6) += fac * A(0, 1) * B(0, 0);
  out(0, 1) += fac * A(0, 1) * B(1, 0);
  out(0, 4) += fac * A(0, 1) * B(2, 0);
  out(0, 8) += fac * A(0, 2) * B(0, 0);
  out(0, 7) += fac * A(0, 2) * B(1, 0);
  out(0, 2) += fac * A(0, 2) * B(2, 0);

  out(3, 0) += fac * A(0, 0) * B(0, 1);
  out(3, 3) += fac * A(0, 0) * B(1, 1);
  out(3, 5) += fac * A(0, 0) * B(2, 1);
  out(3, 6) += fac * A(0, 1) * B(0, 1);
  out(3, 1) += fac * A(0, 1) * B(1, 1);
  out(3, 4) += fac * A(0, 1) * B(2, 1);
  out(3, 8) += fac * A(0, 2) * B(0, 1);
  out(3, 7) += fac * A(0, 2) * B(1, 1);
  out(3, 2) += fac * A(0, 2) * B(2, 1);

  out(5, 0) += fac * A(0, 0) * B(0, 2);
  out(5, 3) += fac * A(0, 0) * B(1, 2);
  out(5, 5) += fac * A(0, 0) * B(2, 2);
  out(5, 6) += fac * A(0, 1) * B(0, 2);
  out(5, 1) += fac * A(0, 1) * B(1, 2);
  out(5, 4) += fac * A(0, 1) * B(2, 2);
  out(5, 8) += fac * A(0, 2) * B(0, 2);
  out(5, 7) += fac * A(0, 2) * B(1, 2);
  out(5, 2) += fac * A(0, 2) * B(2, 2);

  out(6, 0) += fac * A(1, 0) * B(0, 0);
  out(6, 3) += fac * A(1, 0) * B(1, 0);
  out(6, 5) += fac * A(1, 0) * B(2, 0);
  out(6, 6) += fac * A(1, 1) * B(0, 0);
  out(6, 1) += fac * A(1, 1) * B(1, 0);
  out(6, 4) += fac * A(1, 1) * B(2, 0);
  out(6, 8) += fac * A(1, 2) * B(0, 0);
  out(6, 7) += fac * A(1, 2) * B(1, 0);
  out(6, 2) += fac * A(1, 2) * B(2, 0);

  out(1, 0) += fac * A(1, 0) * B(0, 1);
  out(1, 3) += fac * A(1, 0) * B(1, 1);
  out(1, 5) += fac * A(1, 0) * B(2, 1);
  out(1, 6) += fac * A(1, 1) * B(0, 1);
  out(1, 1) += fac * A(1, 1) * B(1, 1);
  out(1, 4) += fac * A(1, 1) * B(2, 1);
  out(1, 8) += fac * A(1, 2) * B(0, 1);
  out(1, 7) += fac * A(1, 2) * B(1, 1);
  out(1, 2) += fac * A(1, 2) * B(2, 1);

  out(4, 0) += fac * A(1, 0) * B(0, 2);
  out(4, 3) += fac * A(1, 0) * B(1, 2);
  out(4, 5) += fac * A(1, 0) * B(2, 2);
  out(4, 6) += fac * A(1, 1) * B(0, 2);
  out(4, 1) += fac * A(1, 1) * B(1, 2);
  out(4, 4) += fac * A(1, 1) * B(2, 2);
  out(4, 8) += fac * A(1, 2) * B(0, 2);
  out(4, 7) += fac * A(1, 2) * B(1, 2);
  out(4, 2) += fac * A(1, 2) * B(2, 2);

  out(8, 0) += fac * A(2, 0) * B(0, 0);
  out(8, 3) += fac * A(2, 0) * B(1, 0);
  out(8, 5) += fac * A(2, 0) * B(2, 0);
  out(8, 6) += fac * A(2, 1) * B(0, 0);
  out(8, 1) += fac * A(2, 1) * B(1, 0);
  out(8, 4) += fac * A(2, 1) * B(2, 0);
  out(8, 8) += fac * A(2, 2) * B(0, 0);
  out(8, 7) += fac * A(2, 2) * B(1, 0);
  out(8, 2) += fac * A(2, 2) * B(2, 0);

  out(7, 0) += fac * A(2, 0) * B(0, 1);
  out(7, 3) += fac * A(2, 0) * B(1, 1);
  out(7, 5) += fac * A(2, 0) * B(2, 1);
  out(7, 6) += fac * A(2, 1) * B(0, 1);
  out(7, 1) += fac * A(2, 1) * B(1, 1);
  out(7, 4) += fac * A(2, 1) * B(2, 1);
  out(7, 8) += fac * A(2, 2) * B(0, 1);
  out(7, 7) += fac * A(2, 2) * B(1, 1);
  out(7, 2) += fac * A(2, 2) * B(2, 1);

  out(2, 0) += fac * A(2, 0) * B(0, 2);
  out(2, 3) += fac * A(2, 0) * B(1, 2);
  out(2, 5) += fac * A(2, 0) * B(2, 2);
  out(2, 6) += fac * A(2, 1) * B(0, 2);
  out(2, 1) += fac * A(2, 1) * B(1, 2);
  out(2, 4) += fac * A(2, 1) * B(2, 2);
  out(2, 8) += fac * A(2, 2) * B(0, 2);
  out(2, 7) += fac * A(2, 2) * B(1, 2);
  out(2, 2) += fac * A(2, 2) * B(2, 2);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
void MAT::MatrixtoStressLikeVoigtNotation(
    LINALG::Matrix<3, 3, T> const& in, LINALG::Matrix<6, 1, T>& out)
{
  for (int i = 0; i < 3; ++i) out(i) = in(i, i);
  out(3) = 0.5 * (in(0, 1) + in(1, 0));
  out(4) = 0.5 * (in(1, 2) + in(2, 1));
  out(5) = 0.5 * (in(0, 2) + in(2, 0));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::MatrixtoStrainLikeVoigtNotation(LINALG::Matrix<3, 3> const& in, LINALG::Matrix<6, 1>& out)
{
  for (int i = 0; i < 3; ++i) out(i) = in(i, i);
  out(3) = in(0, 1) + in(1, 0);
  out(4) = in(1, 2) + in(2, 1);
  out(5) = in(0, 2) + in(2, 0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StressLikeVoigtNotationtoMatrix(LINALG::Matrix<6, 1> const& in, LINALG::Matrix<3, 3>& out)
{
  for (int i = 0; i < 3; ++i) out(i, i) = in(i);
  out(0, 1) = out(1, 0) = in(3);
  out(1, 2) = out(2, 1) = in(4);
  out(0, 2) = out(2, 0) = in(5);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Matrix3x3to9x1(LINALG::Matrix<3, 3> const& in, LINALG::Matrix<9, 1>& out)
{
  for (int i = 0; i < 3; ++i) out(i) = in(i, i);
  out(3) = in(0, 1);
  out(4) = in(1, 2);
  out(5) = in(0, 2);
  out(6) = in(1, 0);
  out(7) = in(2, 1);
  out(8) = in(2, 0);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::InvariantsModified(LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<3, 1>& prinv)
{
  // 1st invariant, trace
  modinv(0) = prinv(0) * std::pow(prinv(2), -1. / 3.);
  // 2nd invariant
  modinv(1) = prinv(1) * std::pow(prinv(2), -2. / 3.);
  // J
  modinv(2) = std::pow(prinv(2), 1. / 2.);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StretchesPrincipal(
    LINALG::Matrix<3, 1>& prstr, LINALG::Matrix<3, 3>& prdir, const LINALG::Matrix<6, 1>& rcg)
{
  // create right Cauchy-Green 2-tensor
  LINALG::Matrix<3, 3> rcgt(false);
  rcgt(0, 0) = rcg(0);
  rcgt(1, 1) = rcg(1);
  rcgt(2, 2) = rcg(2);
  rcgt(0, 1) = rcgt(1, 0) = 0.5 * rcg(3);
  rcgt(1, 2) = rcgt(2, 1) = 0.5 * rcg(4);
  rcgt(2, 0) = rcgt(0, 2) = 0.5 * rcg(5);

  // eigenvalue decomposition
  LINALG::Matrix<3, 3> prstr2;  // squared principal stretches
  LINALG::SYEV(rcgt, prstr2, prdir);

  // THE principal stretches
  for (int al = 0; al < 3; ++al) prstr(al) = std::sqrt(prstr2(al, al));

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StretchesModified(LINALG::Matrix<3, 1>& modstr, const LINALG::Matrix<3, 1>& prstr)
{
  // determinant of deformation gradient
  const double detdefgrad = prstr(0) * prstr(1) * prstr(2);

  // determine modified principal stretches
  modstr.Update(std::pow(detdefgrad, -1.0 / 3.0), prstr);

  return;
}

void MAT::IdentityMatrix(LINALG::Matrix<6, 1>& id)
{
  id.Clear();
  for (unsigned i = 0; i < 3; ++i) id(i) = 1.0;
}

template <unsigned int N>
void MAT::IdentityMatrix(LINALG::Matrix<N, N>& id)
{
  id.Clear();
  for (unsigned i = 0; i < N; ++i) id(i, i) = 1.0;
}

template <MAT::VoigtNotation rows_notation, MAT::VoigtNotation cols_notation>
void MAT::FourthOrderIdentityMatrix(LINALG::Matrix<6, 6>& id)
{
  id.Clear();

  for (unsigned int i = 0; i < 3; ++i) id(i, i) = 1.0;

  for (unsigned int i = 3; i < 6; ++i)
    id(i, i) = 0.5 * ScaleFactor<rows_notation>(i) * ScaleFactor<cols_notation>(i);
}

/*----------------------------------------------------------------------------*/
// initialization of const static members has to be done in the cpp file since
// GCC complains otherwise
const unsigned MAT::IMap::second_[3][3] = {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}};
const unsigned MAT::IMap::fourth_[6][6][4] = {
    {{0, 0, 0, 0}, {0, 0, 1, 1}, {0, 0, 2, 2}, {0, 0, 0, 1}, {0, 0, 1, 2}, {0, 0, 0, 2}},
    {{1, 1, 0, 0}, {1, 1, 1, 1}, {1, 1, 2, 2}, {1, 1, 0, 1}, {1, 1, 1, 2}, {1, 1, 0, 2}},
    {{2, 2, 0, 0}, {2, 2, 1, 1}, {2, 2, 2, 2}, {2, 2, 0, 1}, {2, 2, 1, 2}, {2, 2, 0, 2}},
    {{0, 1, 0, 0}, {0, 1, 1, 1}, {0, 1, 2, 2}, {0, 1, 0, 1}, {0, 1, 1, 2}, {0, 1, 0, 2}},
    {{1, 2, 0, 0}, {1, 2, 1, 1}, {1, 2, 2, 2}, {1, 2, 0, 1}, {1, 2, 1, 2}, {1, 2, 0, 2}},
    {{0, 2, 0, 0}, {0, 2, 1, 1}, {0, 2, 2, 2}, {0, 2, 0, 1}, {0, 2, 1, 2}, {0, 2, 0, 2}}};
template <MAT::VoigtNotation type>
const double MAT::VoigtUtils<type>::unscale_fac_[6] = {1.0, 1.0, 1.0, 0.5, 0.5, 0.5};
template <MAT::VoigtNotation type>
const double MAT::VoigtUtils<type>::scale_fac_[6] = {1.0, 1.0, 1.0, 2.0, 2.0, 2.0};

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <MAT::VoigtNotation type>
void MAT::VoigtUtils<type>::SymmetricOuterProduct(const LINALG::Matrix<3, 1>& vec_a,
    const LINALG::Matrix<3, 1>& vec_b, LINALG::Matrix<6, 1>& ab_ba)
{
  std::fill(ab_ba.A(), ab_ba.A() + 6, 0.0);

  LINALG::Matrix<3, 3> outer_product;
  outer_product.MultiplyNT(vec_a, vec_b);

  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = i; j < 3; ++j)
      ab_ba(IMap::second_[i][j]) += outer_product(i, j) + outer_product(j, i);

  // scale off-diagonal values
  ScaleOffDiagonalVals(ab_ba);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <MAT::VoigtNotation type>
void MAT::VoigtUtils<type>::MultiplyTensorVector(
    const LINALG::Matrix<6, 1>& strain, const LINALG::Matrix<3, 1>& vec, LINALG::Matrix<3, 1>& res)
{
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
    {
      const double fac = UnscaleFactor<type>(IMap::second_[i][j]);
      res(i, 0) += strain(IMap::second_[i][j]) * fac * vec(j, 0);
    }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <MAT::VoigtNotation type>
void MAT::VoigtUtils<type>::PowerOfSymmetricTensor(
    const unsigned pow, const LINALG::Matrix<6, 1>& strain, LINALG::Matrix<6, 1>& strain_pow)
{
  std::copy(strain.A(), strain.A() + 6, strain_pow.A());

  if (pow > 1)
  {
    // unscale the off-diagonal values
    UnscaleOffDiagonalVals(strain_pow);

    LINALG::Matrix<6, 1> prod(false);
    for (unsigned p = 1; p < pow; ++p)
    {
      std::fill(prod.A(), prod.A() + 6, 0.0);

      for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = i; j < 3; ++j)
          for (unsigned k = 0; k < 3; ++k)
            prod(IMap::second_[i][j], 0) += strain_pow(IMap::second_[i][k], 0) *
                                            unscale_fac_[IMap::second_[k][j]] *
                                            strain(IMap::second_[k][j], 0);

      std::copy(prod.A(), prod.A() + 6, strain_pow.A());
    }

    // scale the off-diagonal values again
    ScaleOffDiagonalVals(strain_pow);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <MAT::VoigtNotation type>
void MAT::VoigtUtils<type>::InverseTensor(
    const LINALG::Matrix<6, 1>& tens, LINALG::Matrix<6, 1>& tens_inv)
{
  double det = MAT::VoigtUtils<type>::Determinant(tens);
  tens_inv(0) = (tens(1) * tens(2) -
                    MAT::UnscaleFactor<type>(4) * MAT::UnscaleFactor<type>(4) * tens(4) * tens(4)) /
                det * MAT::ScaleFactor<type>(0);
  tens_inv(1) = (tens(0) * tens(2) -
                    MAT::UnscaleFactor<type>(5) * MAT::UnscaleFactor<type>(5) * tens(5) * tens(5)) /
                det * MAT::ScaleFactor<type>(1);
  tens_inv(2) = (tens(0) * tens(1) -
                    MAT::UnscaleFactor<type>(3) * MAT::UnscaleFactor<type>(3) * tens(3) * tens(3)) /
                det * MAT::ScaleFactor<type>(2);
  tens_inv(3) = (MAT::UnscaleFactor<type>(5) * MAT::UnscaleFactor<type>(4) * tens(5) * tens(4) -
                    MAT::UnscaleFactor<type>(3) * MAT::UnscaleFactor<type>(2) * tens(3) * tens(2)) /
                det * MAT::ScaleFactor<type>(3);
  tens_inv(4) = (MAT::UnscaleFactor<type>(3) * MAT::UnscaleFactor<type>(5) * tens(3) * tens(5) -
                    MAT::UnscaleFactor<type>(0) * MAT::UnscaleFactor<type>(4) * tens(0) * tens(4)) /
                det * MAT::ScaleFactor<type>(4);
  tens_inv(5) = (MAT::UnscaleFactor<type>(3) * MAT::UnscaleFactor<type>(4) * tens(3) * tens(4) -
                    MAT::UnscaleFactor<type>(5) * MAT::UnscaleFactor<type>(1) * tens(5) * tens(1)) /
                det * MAT::ScaleFactor<type>(5);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <MAT::VoigtNotation type>
void MAT::VoigtUtils<type>::ScaleOffDiagonalVals(LINALG::Matrix<6, 1>& strain)
{
  for (unsigned i = 3; i < 6; ++i) strain(i, 0) *= ScaleFactor<type>(i);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <MAT::VoigtNotation type>
void MAT::VoigtUtils<type>::UnscaleOffDiagonalVals(LINALG::Matrix<6, 1>& strain)
{
  for (unsigned i = 3; i < 6; ++i) strain(i, 0) *= UnscaleFactor<type>(i);
}

template <MAT::VoigtNotation type>
void MAT::VoigtUtils<type>::ToStressLike(
    const LINALG::Matrix<6, 1>& vtensor_in, LINALG::Matrix<6, 1>& vtensor_out)
{
  for (unsigned i = 0; i < 6; ++i) vtensor_out(i) = MAT::UnscaleFactor<type>(i) * vtensor_in(i);
}

template <MAT::VoigtNotation type>
void MAT::VoigtUtils<type>::ToStrainLike(
    const LINALG::Matrix<6, 1>& vtensor_in, LINALG::Matrix<6, 1>& vtensor_out)
{
  for (unsigned i = 0; i < 6; ++i)
    vtensor_out(i) =
        MAT::UnscaleFactor<type>(i) * vtensor_in(i) * MAT::ScaleFactor<VoigtNotation::strain>(i);
}

/*----------------------------------------------------------------------------*/
// explicit instantiation of template functions
template void MAT::AddRightNonSymmetricHolzapfelProduct<double>(LINALG::Matrix<6, 9, double>&,
    LINALG::Matrix<3, 3, double> const&, LINALG::Matrix<3, 3, double> const&, double const);
template void MAT::AddRightNonSymmetricHolzapfelProduct<FAD>(LINALG::Matrix<6, 9, FAD>&,
    LINALG::Matrix<3, 3, FAD> const&, LINALG::Matrix<3, 3, FAD> const&, FAD const);
template void MAT::AddRightNonSymmetricHolzapfelProductStrainLike<double>(
    LINALG::Matrix<6, 9, double>& out, LINALG::Matrix<3, 3, double> const& A,
    LINALG::Matrix<3, 3, double> const& B, double const fac);
template void MAT::AddRightNonSymmetricHolzapfelProductStrainLike<FAD>(
    LINALG::Matrix<6, 9, FAD>& out, LINALG::Matrix<3, 3, FAD> const& A,
    LINALG::Matrix<3, 3, FAD> const& B, FAD const fac);
template void MAT::AddtoCmatHolzapfelProduct<double>(
    LINALG::Matrix<6, 6, double>&, const LINALG::Matrix<6, 1, double>&, const double scalar);
template void MAT::AddtoCmatHolzapfelProduct<FAD>(
    LINALG::Matrix<6, 6, FAD>&, const LINALG::Matrix<6, 1, FAD>&, const FAD scalar);
template void MAT::MatrixtoStressLikeVoigtNotation<double>(
    LINALG::Matrix<3, 3, double> const& in, LINALG::Matrix<6, 1, double>& out);
template void MAT::MatrixtoStressLikeVoigtNotation<FAD>(
    LINALG::Matrix<3, 3, FAD> const& in, LINALG::Matrix<6, 1, FAD>& out);

template class MAT::VoigtUtils<MAT::VoigtNotation::strain>;
template class MAT::VoigtUtils<MAT::VoigtNotation::stress>;

template void MAT::IdentityMatrix<6>(LINALG::Matrix<6, 6>& id);
template void MAT::FourthOrderIdentityMatrix<MAT::VoigtNotation::stress,
    MAT::VoigtNotation::stress>(LINALG::Matrix<6, 6>& id);
template void MAT::FourthOrderIdentityMatrix<MAT::VoigtNotation::stress,
    MAT::VoigtNotation::strain>(LINALG::Matrix<6, 6>& id);