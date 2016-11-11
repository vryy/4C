/*----------------------------------------------------------------------*/
/*!
\file material_service.cpp

\brief Interface class for complex materials at Gauss points

\level 1

<pre>
\maintainer Lena Wiechert
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/


//#include "../drt_lib/drt_globalproblem.H"
#include "material_service.H"
#include "matpar_parameter.H"
#include "matpar_bundle.H"


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
void MAT::AddtoCmatHolzapfelProduct(Epetra_SerialDenseMatrix& cmat,
                                    const Epetra_SerialDenseVector& invc,
                                    const double scalar)
{
#ifdef DEBUG
  if (cmat.M()!=6 or cmat.N()!=6 or invc.Length()!=6)
    dserror("Wrong dimensions in function AddtoCmatHolzapfelProduct");
#endif

  // and the 'boeppel-product' for the expression d(invc)/dc (see Holzapfel p. 254)
  cmat(0,0) += scalar * invc(0)*invc(0);
  cmat(0,1) += scalar * invc(3)*invc(3);
  cmat(0,2) += scalar * invc(5)*invc(5);
  cmat(0,3) += scalar * invc(0)*invc(3);
  cmat(0,4) += scalar * invc(3)*invc(5);
  cmat(0,5) += scalar * invc(0)*invc(5);

  cmat(1,0) += scalar * invc(3)*invc(3);
  cmat(1,1) += scalar * invc(1)*invc(1);
  cmat(1,2) += scalar * invc(4)*invc(4);
  cmat(1,3) += scalar * invc(3)*invc(1);
  cmat(1,4) += scalar * invc(1)*invc(4);
  cmat(1,5) += scalar * invc(3)*invc(4);

  cmat(2,0) += scalar * invc(5)*invc(5);
  cmat(2,1) += scalar * invc(4)*invc(4);
  cmat(2,2) += scalar * invc(2)*invc(2);
  cmat(2,3) += scalar * invc(5)*invc(4);
  cmat(2,4) += scalar * invc(4)*invc(2);
  cmat(2,5) += scalar * invc(5)*invc(2);

  cmat(3,0) += scalar * invc(0)*invc(3);
  cmat(3,1) += scalar * invc(3)*invc(1);
  cmat(3,2) += scalar * invc(5)*invc(4);
  cmat(3,3) += scalar * 0.5*( invc(0)*invc(1) + invc(3)*invc(3) );
  cmat(3,4) += scalar * 0.5*( invc(3)*invc(4) + invc(5)*invc(1) );
  cmat(3,5) += scalar * 0.5*( invc(0)*invc(4) + invc(5)*invc(3) );

  cmat(4,0) += scalar * invc(3)*invc(5);
  cmat(4,1) += scalar * invc(1)*invc(4);
  cmat(4,2) += scalar * invc(4)*invc(2);
  cmat(4,3) += scalar * 0.5*( invc(3)*invc(4) + invc(5)*invc(1) );
  cmat(4,4) += scalar * 0.5*( invc(1)*invc(2) + invc(4)*invc(4) );
  cmat(4,5) += scalar * 0.5*( invc(3)*invc(2) + invc(4)*invc(5) );

  cmat(5,0) += scalar * invc(0)*invc(5);
  cmat(5,1) += scalar * invc(3)*invc(4);
  cmat(5,2) += scalar * invc(5)*invc(2);
  cmat(5,3) += scalar * 0.5*( invc(0)*invc(4) + invc(5)*invc(3) );
  cmat(5,4) += scalar * 0.5*( invc(3)*invc(2) + invc(4)*invc(5) );
  cmat(5,5) += scalar * 0.5*( invc(0)*invc(2) + invc(5)*invc(5) );

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
void MAT::AddtoCmatHolzapfelProduct(LINALG::Matrix<6,6>& cmat,
                                    const LINALG::Matrix<6,1>& invc,
                                    const double scalar)
{
#ifdef DEBUG
  if (cmat.M()!=6 or cmat.N()!=6 or invc.M()!=6)
    dserror("Wrong dimensions in function AddtoCmatHolzapfelProduct");
#endif

  // and the 'boeppel-product' for the expression d(invc)/dc (see Holzapfel p. 254)
  cmat(0,0) += scalar * invc(0)*invc(0);
  cmat(0,1) += scalar * invc(3)*invc(3);
  cmat(0,2) += scalar * invc(5)*invc(5);
  cmat(0,3) += scalar * invc(0)*invc(3);
  cmat(0,4) += scalar * invc(3)*invc(5);
  cmat(0,5) += scalar * invc(0)*invc(5);

  cmat(1,0) += scalar * invc(3)*invc(3);
  cmat(1,1) += scalar * invc(1)*invc(1);
  cmat(1,2) += scalar * invc(4)*invc(4);
  cmat(1,3) += scalar * invc(3)*invc(1);
  cmat(1,4) += scalar * invc(1)*invc(4);
  cmat(1,5) += scalar * invc(3)*invc(4);

  cmat(2,0) += scalar * invc(5)*invc(5);
  cmat(2,1) += scalar * invc(4)*invc(4);
  cmat(2,2) += scalar * invc(2)*invc(2);
  cmat(2,3) += scalar * invc(5)*invc(4);
  cmat(2,4) += scalar * invc(4)*invc(2);
  cmat(2,5) += scalar * invc(5)*invc(2);

  cmat(3,0) += scalar * invc(0)*invc(3);
  cmat(3,1) += scalar * invc(3)*invc(1);
  cmat(3,2) += scalar * invc(5)*invc(4);
  cmat(3,3) += scalar * 0.5*( invc(0)*invc(1) + invc(3)*invc(3) );
  cmat(3,4) += scalar * 0.5*( invc(3)*invc(4) + invc(5)*invc(1) );
  cmat(3,5) += scalar * 0.5*( invc(0)*invc(4) + invc(5)*invc(3) );

  cmat(4,0) += scalar * invc(3)*invc(5);
  cmat(4,1) += scalar * invc(1)*invc(4);
  cmat(4,2) += scalar * invc(4)*invc(2);
  cmat(4,3) += scalar * 0.5*( invc(3)*invc(4) + invc(5)*invc(1) );
  cmat(4,4) += scalar * 0.5*( invc(1)*invc(2) + invc(4)*invc(4) );
  cmat(4,5) += scalar * 0.5*( invc(3)*invc(2) + invc(4)*invc(5) );

  cmat(5,0) += scalar * invc(0)*invc(5);
  cmat(5,1) += scalar * invc(3)*invc(4);
  cmat(5,2) += scalar * invc(5)*invc(2);
  cmat(5,3) += scalar * 0.5*( invc(0)*invc(4) + invc(5)*invc(3) );
  cmat(5,4) += scalar * 0.5*( invc(3)*invc(2) + invc(4)*invc(5) );
  cmat(5,5) += scalar * 0.5*( invc(0)*invc(2) + invc(5)*invc(5) );

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
void MAT::ElastSymTensorMultiply(Epetra_SerialDenseMatrix& C,
                                 const double ScalarAB,
                                 const Epetra_SerialDenseMatrix& A,
                                 const Epetra_SerialDenseMatrix& B,
                                 const double ScalarThis)
{
#ifdef DEBUG
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3){
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");
#endif

  // everything in Voigt-Notation
  LINALG::SerialDenseMatrix AVoigt(6,1);
  LINALG::SerialDenseMatrix BVoigt(6,1);

  AVoigt(0,0) = A(0,0); AVoigt(1,0) = A(1,1); AVoigt(2,0) = A(2,2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3,0) = A(1,0); AVoigt(4,0) = A(2,1); AVoigt(5,0) = A(2,0);

  BVoigt(0,0) = B(0,0); BVoigt(1,0) = B(1,1); BVoigt(2,0) = B(2,2);
  BVoigt(3,0) = B(1,0); BVoigt(4,0) = B(2,1); BVoigt(5,0) = B(2,0);

  C.Multiply('N','T',ScalarAB,AVoigt,BVoigt,ScalarThis);

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
void MAT::ElastSymTensorMultiply(LINALG::Matrix<6,6>& C,
                                 const double ScalarAB,
                                 const LINALG::Matrix<3,3>& A,
                                 const LINALG::Matrix<3,3>& B,
                                 const double ScalarThis)
{
  // everything in Voigt-Notation
  LINALG::Matrix<6,1> AVoigt;
  LINALG::Matrix<6,1> BVoigt;

  AVoigt(0,0) = A(0,0); AVoigt(1,0) = A(1,1); AVoigt(2,0) = A(2,2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3,0) = A(1,0); AVoigt(4,0) = A(2,1); AVoigt(5,0) = A(2,0);

  BVoigt(0,0) = B(0,0); BVoigt(1,0) = B(1,1); BVoigt(2,0) = B(2,2);
  BVoigt(3,0) = B(1,0); BVoigt(4,0) = B(2,1); BVoigt(5,0) = B(2,0);

  C.MultiplyNT(ScalarAB,AVoigt,BVoigt,ScalarThis);

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
void MAT::ElastSymTensorMultiplyAddSym(Epetra_SerialDenseMatrix& C,
                                 const double ScalarAB,
                                 const Epetra_SerialDenseMatrix& A,
                                 const Epetra_SerialDenseMatrix& B,
                                 const double ScalarThis)
{
#ifdef DEBUG
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3){
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");
#endif

  // everything in Voigt-Notation
  LINALG::SerialDenseMatrix AVoigt(6,1);
  LINALG::SerialDenseMatrix BVoigt(6,1);

  AVoigt(0,0) = A(0,0); AVoigt(1,0) = A(1,1); AVoigt(2,0) = A(2,2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3,0) = A(1,0); AVoigt(4,0) = A(2,1); AVoigt(5,0) = A(2,0);

  BVoigt(0,0) = B(0,0); BVoigt(1,0) = B(1,1); BVoigt(2,0) = B(2,2);
  BVoigt(3,0) = B(1,0); BVoigt(4,0) = B(2,1); BVoigt(5,0) = B(2,0);

  C.Multiply('N','T',ScalarAB,AVoigt,BVoigt,ScalarThis);
  C.Multiply('N','T',ScalarAB,BVoigt,AVoigt,1.0);

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
void MAT::ElastSymTensorMultiplyAddSym(LINALG::Matrix<6,6>& C,
        const double ScalarAB,
        const LINALG::Matrix<3,3>& A,
        const LINALG::Matrix<3,3>& B,
        const double ScalarThis)
{
#ifdef DEBUG
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3){
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");
#endif

  // everything in Voigt-Notation
  LINALG::Matrix<6,1> AVoigt;
  LINALG::Matrix<6,1> BVoigt;

  AVoigt(0,0) = A(0,0); AVoigt(1,0) = A(1,1); AVoigt(2,0) = A(2,2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3,0) = A(1,0); AVoigt(4,0) = A(2,1); AVoigt(5,0) = A(2,0);

  BVoigt(0,0) = B(0,0); BVoigt(1,0) = B(1,1); BVoigt(2,0) = B(2,2);
  BVoigt(3,0) = B(1,0); BVoigt(4,0) = B(2,1); BVoigt(5,0) = B(2,0);

  C.MultiplyNT(ScalarAB,AVoigt,BVoigt,ScalarThis);
  C.MultiplyNT(ScalarAB,BVoigt,AVoigt,1.0);

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
void MAT::ElastSymTensor_o_Multiply(Epetra_SerialDenseMatrix& C,
                                 const double ScalarAB,
                                 const Epetra_SerialDenseMatrix& A,
                                 const Epetra_SerialDenseMatrix& B,
                                 const double ScalarThis)
{
#ifdef DEBUG
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3){
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
//  C(3,3)= ScalarThis*C(3,3) + ScalarAB * 0.5*(A(1,1)*B(0,0) + 2.0*A(1,2)*B(0,0) + A(2,2)*B(0,0));
//  C(3,4)= ScalarThis*C(3,4) + ScalarAB * 0.5*(A(1,1)*B(0,1) + 2.0*A(1,2)*B(0,1) + A(2,2)*B(0,1));
//  C(3,5)= ScalarThis*C(3,5) + ScalarAB * 0.5*(A(1,1)*B(0,2) + 2.0*A(1,2)*B(0,2) + A(2,2)*B(0,2));
//
//  C(4,0)= ScalarThis*C(4,0) + ScalarAB * (A(1,0)*B(1,0) + A(2,0)*B(1,0));
//  C(4,1)= ScalarThis*C(4,1) + ScalarAB * (A(1,0)*B(1,1) + A(2,0)*B(1,1));
//  C(4,2)= ScalarThis*C(4,2) + ScalarAB * (A(1,0)*B(1,2) + A(2,0)*B(1,2));
//  C(4,3)= ScalarThis*C(4,3) + ScalarAB * 0.5*(A(1,1)*B(1,0) + 2.0*A(1,2)*B(1,0) + A(2,2)*B(1,0));
//  C(4,4)= ScalarThis*C(4,4) + ScalarAB * 0.5*(A(1,1)*B(1,1) + 2.0*A(1,2)*B(1,1) + A(2,2)*B(1,1));
//  C(4,5)= ScalarThis*C(4,5) + ScalarAB * 0.5*(A(1,1)*B(1,2) + 2.0*A(1,2)*B(1,2) + A(2,2)*B(1,2));
//
//  C(5,0)= ScalarThis*C(5,0) + ScalarAB * (A(1,0)*B(2,0) + A(2,0)*B(2,0));
//  C(5,1)= ScalarThis*C(5,1) + ScalarAB * (A(1,0)*B(2,1) + A(2,0)*B(2,1));
//  C(5,2)= ScalarThis*C(5,2) + ScalarAB * (A(1,0)*B(2,2) + A(2,0)*B(2,2));
//  C(5,3)= ScalarThis*C(5,3) + ScalarAB * 0.5*(A(1,1)*B(2,0) + 2.0*A(1,2)*B(2,0) + A(2,2)*B(2,0));
//  C(5,4)= ScalarThis*C(5,4) + ScalarAB * 0.5*(A(1,1)*B(2,1) + 2.0*A(1,2)*B(2,1) + A(2,2)*B(2,1));
//  C(5,5)= ScalarThis*C(5,5) + ScalarAB * 0.5*(A(1,1)*B(2,2) + 2.0*A(1,2)*B(2,2) + A(2,2)*B(2,2));

  C(0,0)= ScalarThis*C(0,0) + ScalarAB * 0.5 * (A(0,0)*B(0,0) + A(0,0)*B(0,0));//C1111
  C(0,1)= ScalarThis*C(0,1) + ScalarAB * 0.5 * (A(0,1)*B(0,1) + A(0,1)*B(0,1));//C1122
  C(0,2)= ScalarThis*C(0,2) + ScalarAB * 0.5 * (A(0,2)*B(0,2) + A(0,2)*B(0,2));//C1133
  C(0,3)= ScalarThis*C(0,3) + ScalarAB * 0.5 * (A(0,0)*B(0,1) + A(0,1)*B(0,0));//C1112
  C(0,4)= ScalarThis*C(0,4) + ScalarAB * 0.5 * (A(0,1)*B(0,2) + A(0,2)*B(0,1));//C1123
  C(0,5)= ScalarThis*C(0,5) + ScalarAB * 0.5 * (A(0,0)*B(0,2) + A(0,2)*B(0,0));//C1113

  C(1,0)= ScalarThis*C(1,0) + ScalarAB * 0.5 * (A(1,0)*B(1,0) + A(1,0)*B(1,0));//C2211
  C(1,1)= ScalarThis*C(1,1) + ScalarAB * 0.5 * (A(1,1)*B(1,1) + A(1,1)*B(1,1));//C2222
  C(1,2)= ScalarThis*C(1,2) + ScalarAB * 0.5 * (A(1,2)*B(1,2) + A(1,2)*B(1,2));//C2233
  C(1,3)= ScalarThis*C(1,3) + ScalarAB * 0.5 * (A(1,0)*B(1,1) + A(1,1)*B(1,0));//C2212
  C(1,4)= ScalarThis*C(1,4) + ScalarAB * 0.5 * (A(1,1)*B(1,2) + A(1,2)*B(1,1));//C2223
  C(1,5)= ScalarThis*C(1,5) + ScalarAB * 0.5 * (A(1,0)*B(1,2) + A(1,2)*B(1,0));//C2213

  C(2,0)= ScalarThis*C(2,0) + ScalarAB * 0.5 * (A(2,0)*B(2,0) + A(2,0)*B(2,0));//C3311
  C(2,1)= ScalarThis*C(2,1) + ScalarAB * 0.5 * (A(2,1)*B(2,1) + A(2,1)*B(2,1));//C3322
  C(2,2)= ScalarThis*C(2,2) + ScalarAB * 0.5 * (A(2,2)*B(2,2) + A(2,2)*B(2,2));//C3333
  C(2,3)= ScalarThis*C(2,3) + ScalarAB * 0.5 * (A(2,1)*B(2,1) + A(2,1)*B(2,0));//C3312
  C(2,4)= ScalarThis*C(2,4) + ScalarAB * 0.5 * (A(2,1)*B(2,2) + A(2,2)*B(2,1));//C3323
  C(2,5)= ScalarThis*C(2,5) + ScalarAB * 0.5 * (A(2,0)*B(2,2) + A(2,2)*B(2,0));//C3313

  C(3,0)= ScalarThis*C(3,0) + ScalarAB * 0.5 * (A(0,0)*B(1,0) + A(0,0)*B(1,0));//C1211
  C(3,1)= ScalarThis*C(3,1) + ScalarAB * 0.5 * (A(0,1)*B(1,1) + A(0,1)*B(1,1));//C1222
  C(3,2)= ScalarThis*C(3,2) + ScalarAB * 0.5 * (A(0,2)*B(1,2) + A(0,2)*B(1,2));//C1233
  C(3,3)= ScalarThis*C(3,3) + ScalarAB * 0.5 * (A(0,0)*B(1,1) + A(0,1)*B(1,0));//C1212
  C(3,4)= ScalarThis*C(3,4) + ScalarAB * 0.5 * (A(0,1)*B(1,2) + A(0,2)*B(1,1));//C1223
  C(3,5)= ScalarThis*C(3,5) + ScalarAB * 0.5 * (A(0,0)*B(1,2) + A(0,2)*B(1,0));//C1213

  C(4,0)= ScalarThis*C(4,0) + ScalarAB * 0.5 * (A(1,0)*B(2,0) + A(1,0)*B(2,0));//C2311
  C(4,1)= ScalarThis*C(4,1) + ScalarAB * 0.5 * (A(1,1)*B(2,1) + A(1,1)*B(2,1));//C2322
  C(4,2)= ScalarThis*C(4,2) + ScalarAB * 0.5 * (A(1,2)*B(2,2) + A(1,2)*B(2,2));//C2333
  C(4,3)= ScalarThis*C(4,3) + ScalarAB * 0.5 * (A(1,0)*B(2,1) + A(1,1)*B(2,0));//C2312
  C(4,4)= ScalarThis*C(4,4) + ScalarAB * 0.5 * (A(1,1)*B(2,2) + A(1,2)*B(2,1));//C2323
  C(4,5)= ScalarThis*C(4,5) + ScalarAB * 0.5 * (A(1,0)*B(2,2) + A(1,2)*B(2,0));//C2313

  C(5,0)= ScalarThis*C(5,0) + ScalarAB * 0.5 * (A(0,0)*B(2,0) + A(0,0)*B(2,0));//C1311
  C(5,1)= ScalarThis*C(5,1) + ScalarAB * 0.5 * (A(0,1)*B(2,1) + A(0,1)*B(2,1));//C1322
  C(5,2)= ScalarThis*C(5,2) + ScalarAB * 0.5 * (A(0,2)*B(2,2) + A(0,2)*B(2,2));//C1333
  C(5,3)= ScalarThis*C(5,3) + ScalarAB * 0.5 * (A(0,0)*B(2,1) + A(0,1)*B(2,0));//C1312
  C(5,4)= ScalarThis*C(5,4) + ScalarAB * 0.5 * (A(0,1)*B(2,2) + A(0,2)*B(2,1));//C1323
  C(5,5)= ScalarThis*C(5,5) + ScalarAB * 0.5 * (A(0,0)*B(2,2) + A(0,2)*B(2,0));//C1313

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
void MAT::ElastSymTensor_o_Multiply(LINALG::Matrix<6,6>& C,
                                    const double ScalarAB,
                                    const LINALG::Matrix<3,3>& A,
                                    const LINALG::Matrix<3,3>& B,
                                    const double ScalarThis)
{
  // To keep the code somewhat shorter I removed the explanation comment here,
  // but they are still in the Epetra Version
  const double ScalarABhalf = ScalarAB * 0.5;
  C(0,0)= ScalarThis*C(0,0) + ScalarAB     * (A(0,0)*B(0,0)                );//C1111
  C(0,1)= ScalarThis*C(0,1) + ScalarAB     * (A(0,1)*B(0,1)                );//C1122
  C(0,2)= ScalarThis*C(0,2) + ScalarAB     * (A(0,2)*B(0,2)                );//C1133
  C(0,3)= ScalarThis*C(0,3) + ScalarABhalf * (A(0,0)*B(0,1) + A(0,1)*B(0,0));//C1112
  C(0,4)= ScalarThis*C(0,4) + ScalarABhalf * (A(0,1)*B(0,2) + A(0,2)*B(0,1));//C1123
  C(0,5)= ScalarThis*C(0,5) + ScalarABhalf * (A(0,0)*B(0,2) + A(0,2)*B(0,0));//C1113

  C(1,0)= ScalarThis*C(1,0) + ScalarAB     * (A(1,0)*B(1,0)                );//C2211
  C(1,1)= ScalarThis*C(1,1) + ScalarAB     * (A(1,1)*B(1,1)                );//C2222
  C(1,2)= ScalarThis*C(1,2) + ScalarAB     * (A(1,2)*B(1,2)                );//C2233
  C(1,3)= ScalarThis*C(1,3) + ScalarABhalf * (A(1,0)*B(1,1) + A(1,1)*B(1,0));//C2212
  C(1,4)= ScalarThis*C(1,4) + ScalarABhalf * (A(1,1)*B(1,2) + A(1,2)*B(1,1));//C2223
  C(1,5)= ScalarThis*C(1,5) + ScalarABhalf * (A(1,0)*B(1,2) + A(1,2)*B(1,0));//C2213

  C(2,0)= ScalarThis*C(2,0) + ScalarAB     * (A(2,0)*B(2,0)                );//C3311
  C(2,1)= ScalarThis*C(2,1) + ScalarAB     * (A(2,1)*B(2,1)                );//C3322
  C(2,2)= ScalarThis*C(2,2) + ScalarAB     * (A(2,2)*B(2,2)                );//C3333
  C(2,3)= ScalarThis*C(2,3) + ScalarABhalf * (A(2,1)*B(2,1) + A(2,1)*B(2,0));//C3312
  C(2,4)= ScalarThis*C(2,4) + ScalarABhalf * (A(2,1)*B(2,2) + A(2,2)*B(2,1));//C3323
  C(2,5)= ScalarThis*C(2,5) + ScalarABhalf * (A(2,0)*B(2,2) + A(2,2)*B(2,0));//C3313

  C(3,0)= ScalarThis*C(3,0) + ScalarAB     * (A(0,0)*B(1,0)                );//C1211
  C(3,1)= ScalarThis*C(3,1) + ScalarAB     * (A(0,1)*B(1,1)                );//C1222
  C(3,2)= ScalarThis*C(3,2) + ScalarAB     * (A(0,2)*B(1,2)                );//C1233
  C(3,3)= ScalarThis*C(3,3) + ScalarABhalf * (A(0,0)*B(1,1) + A(0,1)*B(1,0));//C1212
  C(3,4)= ScalarThis*C(3,4) + ScalarABhalf * (A(0,1)*B(1,2) + A(0,2)*B(1,1));//C1223
  C(3,5)= ScalarThis*C(3,5) + ScalarABhalf * (A(0,0)*B(1,2) + A(0,2)*B(1,0));//C1213

  C(4,0)= ScalarThis*C(4,0) + ScalarAB     * (A(1,0)*B(2,0)                );//C2311
  C(4,1)= ScalarThis*C(4,1) + ScalarAB     * (A(1,1)*B(2,1)                );//C2322
  C(4,2)= ScalarThis*C(4,2) + ScalarAB     * (A(1,2)*B(2,2)                );//C2333
  C(4,3)= ScalarThis*C(4,3) + ScalarABhalf * (A(1,0)*B(2,1) + A(1,1)*B(2,0));//C2312
  C(4,4)= ScalarThis*C(4,4) + ScalarABhalf * (A(1,1)*B(2,2) + A(1,2)*B(2,1));//C2323
  C(4,5)= ScalarThis*C(4,5) + ScalarABhalf * (A(1,0)*B(2,2) + A(1,2)*B(2,0));//C2313

  C(5,0)= ScalarThis*C(5,0) + ScalarAB     * (A(0,0)*B(2,0)                );//C1311
  C(5,1)= ScalarThis*C(5,1) + ScalarAB     * (A(0,1)*B(2,1)                );//C1322
  C(5,2)= ScalarThis*C(5,2) + ScalarAB     * (A(0,2)*B(2,2)                );//C1333
  C(5,3)= ScalarThis*C(5,3) + ScalarABhalf * (A(0,0)*B(2,1) + A(0,1)*B(2,0));//C1312
  C(5,4)= ScalarThis*C(5,4) + ScalarABhalf * (A(0,1)*B(2,2) + A(0,2)*B(2,1));//C1323
  C(5,5)= ScalarThis*C(5,5) + ScalarABhalf * (A(0,0)*B(2,2) + A(0,2)*B(2,0));//C1313

  return;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::VolumetrifyAndIsochorify(LINALG::Matrix<6,1>* pk2vol,
                                   LINALG::Matrix<6,6>* cvol,
                                   LINALG::Matrix<6,1>* pk2iso,
                                   LINALG::Matrix<6,6>* ciso,
                                   const LINALG::Matrix<6,1>& gl,
                                   const LINALG::Matrix<6,1>& pk2,
                                   const LINALG::Matrix<6,6>& cmat)
{
  // useful call?
#ifdef DEBUG
  if ( (pk2vol == NULL) and (cvol == NULL) and (pk2iso == NULL) and (ciso == NULL) )
    dserror("Useful call? Apparently you do not want to compute anything");
#endif

  // right Cauchy--Green tensor
  // REMARK: stored in _strain_-like 6-Voigt vector
  LINALG::Matrix<6,1> rcg(gl);
  rcg.Scale(2.0);
  for (int i=0; i<3; i++) rcg(i) += 1.0;

  // third invariant (determinant) of right Cauchy--Green strains
  const double rcg3rd = rcg(0)*rcg(1)*rcg(2)
                      + 0.25 * rcg(3)*rcg(4)*rcg(5)
                      - 0.25 * rcg(1)*rcg(5)*rcg(5)
                      - 0.25 * rcg(2)*rcg(3)*rcg(3)
                      - 0.25 * rcg(0)*rcg(4)*rcg(4);

  // inverse right Cauchy--Green tensor C^{-1}
  // REMARK: stored in as _stress_ 6-Voigt vector
  LINALG::Matrix<6,1> icg(false);
  icg(0) = (rcg(1)*rcg(2) - 0.25*rcg(4)*rcg(4))/rcg3rd;  // (C^{-1})^{11}
  icg(1) = (rcg(0)*rcg(2) - 0.25*rcg(5)*rcg(5))/rcg3rd;  // (C^{-1})^{22}
  icg(2) = (rcg(0)*rcg(1) - 0.25*rcg(3)*rcg(3))/rcg3rd;  // (C^{-1})^{33}
  icg(3) = (0.25*rcg(5)*rcg(4) - 0.5*rcg(3)*rcg(2))/rcg3rd;  // (C^{-1})^{12}
  icg(4) = (0.25*rcg(3)*rcg(5) - 0.5*rcg(0)*rcg(4))/rcg3rd;  // (C^{-1})^{23}
  icg(5) = (0.25*rcg(3)*rcg(4) - 0.5*rcg(5)*rcg(1))/rcg3rd;  // (C^{-1})^{31}

  // double contraction of 2nd Piola--Kirchhoff stress and right Cauchy--Green strain,
  // i.e. in index notation S^{AB} C_{AB}
  // REMARK: equal to S^T C, because S is stress-like and C is strain-like 6-Voigt vector
  const double pk2rcg = pk2.Dot(rcg);

  // stress splitting
  {
    // volumetric 2nd Piola--Kirchhoff stress
    LINALG::Matrix<6,1> pk2vol_(false);
    if (pk2vol != NULL)
      pk2vol_.SetView(*pk2vol);
    pk2vol_.Update(pk2rcg/3.0, icg);

    // isochoric 2nd Piola--Kirchhoff stress
    // S^{AB}_iso = S^{AB} - S^{AB}_{vol}
    if (pk2iso != NULL)
      pk2iso->Update(1.0, pk2, -1.0, pk2vol_);
  }

  // elasticity tensor splitting
  {
    // 'linearised' 2nd Piola--Kirchhoff stress
    // S^{CD}_lin = S^{CD} + 1/2 C_{AB} C^{ABCD}
    LINALG::Matrix<6,1> pk2lin(pk2);
    pk2lin.MultiplyTN(0.5, cmat, rcg, 1.0);  // transpose on purpose

    // volumetric part of constitutive tensor
    // C^{ABCD}_vol = 2/3 (C^{-1})^{AB} S^{CD}_lin
    //              - 2/3 (S^{EF} C_{EF}) ( 1/2 (
    //                (C^{-1})^{AC} (C^{-1})^{BD} + (C^{-1})^{AD} (C^{-1})^{BC}
    //              ) )
    LINALG::Matrix<6,6> cvol_(false);
    if (cvol != NULL)
      cvol_.SetView(*cvol);
    cvol_.MultiplyNT(2.0/3.0, icg, pk2lin);
    AddtoCmatHolzapfelProduct(cvol_, icg, -2.0/3.0*pk2rcg);

    // isochoric part of constitutive tensor
    // C^{ABCD}_iso = C^{ABCD} - C^{ABCD}_vol
    if (ciso != NULL)
      ciso->Update(1.0, cmat, -1.0, cvol_);
  }

  //
  return;
}

/*----------------------------------------------------------------------*
 |  add dX^2/dX to cmat                                     seitz 07/13 |
 *----------------------------------------------------------------------*/
void MAT::AddToCmatDerivTensorSquare(LINALG::Matrix<6,6>& C, double ScalarDX2, LINALG::Matrix<3,3> X, double ScalarThis)
{
  C(0,0)= ScalarThis*C(0,0) + ScalarDX2 * 2.*X(0,0);//C1111
  C(0,1)= ScalarThis*C(0,1)  ;//C1122
  C(0,2)= ScalarThis*C(0,2)  ;//C1133
  C(0,3)= ScalarThis*C(0,3) + ScalarDX2 * X(0,1);//C1112
  C(0,4)= ScalarThis*C(0,4)  ;//C1123
  C(0,5)= ScalarThis*C(0,5) + ScalarDX2 * X(0,2);//C1113

  C(1,0)= ScalarThis*C(1,0)  ;//C2211
  C(1,1)= ScalarThis*C(1,1) + ScalarDX2 * 2.*X(1,1);//C2222
  C(1,2)= ScalarThis*C(1,2)  ;//C2233
  C(1,3)= ScalarThis*C(1,3) + ScalarDX2 * X(0,1);//C2212
  C(1,4)= ScalarThis*C(1,4) + ScalarDX2 * X(1,2) ;//C2223
  C(1,5)= ScalarThis*C(1,5)  ;//C2213

  C(2,0)= ScalarThis*C(2,0)   ;//C3311
  C(2,1)= ScalarThis*C(2,1)  ;//C3322
  C(2,2)= ScalarThis*C(2,2) + ScalarDX2 * 2.*X(2,2)  ;//C3333
  C(2,3)= ScalarThis*C(2,3) ;//C3312
  C(2,4)= ScalarThis*C(2,4) + ScalarDX2 * X(1,2) ;//C3323
  C(2,5)= ScalarThis*C(2,5) + ScalarDX2 * X(0,2);//C3313

  C(3,0)= ScalarThis*C(3,0) + ScalarDX2 * X(0,1) ;//C1211
  C(3,1)= ScalarThis*C(3,1) + ScalarDX2 * X(0,1) ;//C1222
  C(3,2)= ScalarThis*C(3,2)   ;//C1233
  C(3,3)= ScalarThis*C(3,3) + ScalarDX2 * 0.5*(X(0,0)+X(1,1));//C1212
  C(3,4)= ScalarThis*C(3,4) + ScalarDX2 * 0.5*X(0,2);//C1223
  C(3,5)= ScalarThis*C(3,5) + ScalarDX2 * 0.5 * X(1,2)  ;//C1213

  C(4,0)= ScalarThis*C(4,0)  ;//C2311
  C(4,1)= ScalarThis*C(4,1) + ScalarDX2 * X(1,2) ;//C2322
  C(4,2)= ScalarThis*C(4,2) + ScalarDX2 * X(1,2) ;//C2333
  C(4,3)= ScalarThis*C(4,3) + ScalarDX2 * 0.5*X(0,2) ;//C2312
  C(4,4)= ScalarThis*C(4,4) + ScalarDX2 * 0.5*(X(1,1)+X(2,2)) ;//C2323
  C(4,5)= ScalarThis*C(4,5) + ScalarDX2 * 0.5*X(0,1) ;//C2313

  C(5,0)= ScalarThis*C(5,0) + ScalarDX2 * X(0,2) ;//C1311
  C(5,1)= ScalarThis*C(5,1)   ;//C1322
  C(5,2)= ScalarThis*C(5,2) + ScalarDX2 * X(0,2) ;//C1333
  C(5,3)= ScalarThis*C(5,3) + ScalarDX2 * 0.5*X(1,2);//C1312
  C(5,4)= ScalarThis*C(5,4) + ScalarDX2 * 0.5*X(0,1);//C1323
  C(5,5)= ScalarThis*C(5,5) + ScalarDX2 * 0.5*(X(2,2) + X(0,0)) ;//C1313

  return;
}


/*----------------------------------------------------------------------*
 |                                                          seitz 09/13 |
 *----------------------------------------------------------------------*/
void MAT::AddSymmetricHolzapfelProduct(LINALG::Matrix<6,6>& X,
    const LINALG::Matrix<3,3>& A, const LINALG::Matrix<3,3>& B, double fac)
{
  X(0,0) += 4 * fac * A(0,0) * B(0,0);
  X(0,3) += fac * (2 * A(0,0) * B(1,0) + 2 * A(1,0) * B(0,0));
  X(0,5) += fac * (2 * A(0,0) * B(2,0) + 2 * A(2,0) * B(0,0));
  X(0,1) += 4 * fac * A(1,0) * B(1,0);
  X(0,4) += fac * (2 * A(1,0) * B(2,0) + 2 * A(2,0) * B(1,0));
  X(0,2) += 4 * fac * A(2,0) * B(2,0);

  X(3,0) += fac * (2 * A(0,0) * B(0,1) + 2 * A(0,1) * B(0,0));
  X(3,3) += fac * (A(0,0) * B(1,1) + A(1,0) * B(0,1) + A(1,1) * B(0,0) + A(0,1) * B(1,0));
  X(3,5) += fac * (A(0,0) * B(2,1) + A(2,0) * B(0,1) + A(2,1) * B(0,0) + A(0,1) * B(2,0));
  X(3,1) += fac * (2 * A(1,0) * B(1,1) + 2 * A(1,1) * B(1,0));
  X(3,4) += fac * (A(1,0) * B(2,1) + A(2,0) * B(1,1) + A(2,1) * B(1,0) + A(1,1) * B(2,0));
  X(3,2) += fac * (2 * A(2,0) * B(2,1) + 2 * A(2,1) * B(2,0));

  X(5,0) += fac * (2 * A(0,0) * B(0,2) + 2 * A(0,2) * B(0,0));
  X(5,3) += fac * (A(0,0) * B(1,2) + A(1,0) * B(0,2) + A(1,2) * B(0,0) + A(0,2) * B(1,0));
  X(5,5) += fac * (A(0,0) * B(2,2) + A(2,0) * B(0,2) + A(2,2) * B(0,0) + A(0,2) * B(2,0));
  X(5,1) += fac * (2 * A(1,0) * B(1,2) + 2 * A(1,2) * B(1,0));
  X(5,4) += fac * (A(1,0) * B(2,2) + A(2,0) * B(1,2) + A(2,2) * B(1,0) + A(1,2) * B(2,0));
  X(5,2) += fac * (2 * A(2,0) * B(2,2) + 2 * A(2,2) * B(2,0));

  X(1,0) += 4 * fac * A(0,1) * B(0,1);
  X(1,3) += fac * (2 * A(0,1) * B(1,1) + 2 * A(1,1) * B(0,1));
  X(1,5) += fac * (2 * A(0,1) * B(2,1) + 2 * A(2,1) * B(0,1));
  X(1,1) += 4 * fac * A(1,1) * B(1,1);
  X(1,4) += fac * (2 * A(1,1) * B(2,1) + 2 * A(2,1) * B(1,1));
  X(1,2) += 4 * fac * A(2,1) * B(2,1);

  X(4,0) += fac * (2 * A(0,1) * B(0,2) + 2 * A(0,2) * B(0,1));
  X(4,3) += fac * (A(0,1) * B(1,2) + A(1,1) * B(0,2) + A(1,2) * B(0,1) + A(0,2) * B(1,1));
  X(4,5) += fac * (A(0,1) * B(2,2) + A(2,1) * B(0,2) + A(2,2) * B(0,1) + A(0,2) * B(2,1));
  X(4,1) += fac * (2 * A(1,1) * B(1,2) + 2 * A(1,2) * B(1,1));
  X(4,4) += fac * (A(1,1) * B(2,2) + A(2,1) * B(1,2) + A(2,2) * B(1,1) + A(1,2) * B(2,1));
  X(4,2) += fac * (2 * A(2,1) * B(2,2) + 2 * A(2,2) * B(2,1));

  X(2,0) += 4 * fac * A(0,2) * B(0,2);
  X(2,3) += fac * (2 * A(0,2) * B(1,2) + 2 * A(1,2) * B(0,2));
  X(2,5) += fac * (2 * A(0,2) * B(2,2) + 2 * A(2,2) * B(0,2));
  X(2,1) += 4 * fac * A(1,2) * B(1,2);
  X(2,4) += fac * (2 * A(1,2) * B(2,2) + 2 * A(2,2) * B(1,2));
  X(2,2) += 4 * fac * A(2,2) * B(2,2);

  return;
}
