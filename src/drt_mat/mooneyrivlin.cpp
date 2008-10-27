/*!----------------------------------------------------------------------
\file mooneyrivlin.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "Epetra_SerialDenseSolver.h"
#include "mooneyrivlin.H"

extern struct _MATERIAL *mat; ///< C-style material struct


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)     maf 04/08|
 *----------------------------------------------------------------------*/
MAT::MooneyRivlin::MooneyRivlin()
  : matdata_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)      maf 04/08|
 *----------------------------------------------------------------------*/
MAT::MooneyRivlin::MooneyRivlin(MATERIAL* matdata)
  : matdata_(matdata)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)     maf 04/08|
 *----------------------------------------------------------------------*/
void MAT::MooneyRivlin::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;   // pointer difference to reach 0-entry
  AddtoPack(data,matdata);
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)     maf 04/08|
 *----------------------------------------------------------------------*/
void MAT::MooneyRivlin::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matdata
  int matdata;
  ExtractfromPack(position,data,matdata);
  matdata_ = &mat[matdata];     // unpack pointer to my specific matdata_

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}


/*----------------------------------------------------------------------*
 |  Return density                                (public)     maf 04/07|
 *----------------------------------------------------------------------*/
double MAT::MooneyRivlin::Density()
{
  return matdata_->m.mooneyrivlin->density;  // density, returned to evaluate mass matrix
}


/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)     maf 04/08|
 *----------------------------------------------------------------------*

Mooney-Rivlin type nearly incompressible, hyperelastic 3D material law.

The underlying strain-energy function is ('a' is 'alpha' and J = I3^(1/2)):

W = sum_(p=1)^2 ( mu_p / a_p (l1^a_p + l2^a_p + l3^a_p) - mu_p ln(J) )
  + Lambda/4 * (J^2 - 1 - 2 ln(J) )

which can be expressed in terms of invariants I1 and I2 as:

W = mu_1/a_1 (I1 - 3)  +  mu_2/a_2 (I2 - 3)  -  mu_1 ln(J)  -  mu_2 ln(J)
  + Lambda/4 * (J^2 - 1 - 2 ln(J) )

For references see Holzapfel p. 245, Simo&Miehe 1992, Klinkel et al 2007.

Parameters are mu_1, mu_2, a_1, a_2 and Lambda = kappa - 2/3 mu as penalty factor
to effectuate incompressibility (shear modulus mu = (mu_1 a_1 + mu_2 a_2)/2 )

Mind that it is not stress-free in reference configuration!
*/

void MAT::MooneyRivlin::Evaluate(const Epetra_SerialDenseVector* glstrain_e,
                                       Epetra_SerialDenseMatrix* cmat_e,
                                       Epetra_SerialDenseVector* stress_e)

{
  // this is temporary as long as the material does not have a 
  // FixedSizeSerialDenseMatrix-type interface
  const LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> glstrain(glstrain_e->A(),true);
        LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,NUM_STRESS_3D> cmat(cmat_e->A(),true);
        LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> stress(stress_e->A(),true);

        // get material parameters
  const double m1  = matdata_->m.mooneyrivlin->mu1;
  const double m2  = matdata_->m.mooneyrivlin->mu2;
  const double c1  = m1 / matdata_->m.mooneyrivlin->alpha1;  // c1 = mu1/a1
  const double c2  = m2 / matdata_->m.mooneyrivlin->alpha2;  // c2 = mu2/a2
  const double pen = matdata_->m.mooneyrivlin->penalty;
  
  // right Cauchy-Green Tensor  C = 2 * E + I
  // build identity tensor I
  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> C(glstrain);
  C.Scale(2.0);
  C += Id;

  // invariants
  const double I1 = C(0) + C(1) + C(2);  // 1st invariant, trace
  const double I3 = C(0)*C(1)*C(2)
        + 0.25 * C(3)*C(4)*C(5)
        - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3)
        - 0.25 * C(0)*C(4)*C(4);    // 3rd invariant, determinant

  // invert C
  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> Cinv(false);

  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);
  Cinv.Scale(1.0/I3);

  /* //Compute strain-energy function W
  // I2 = 1/2 (trace(C)^2 - trace(CxC))
  Epetra_SerialDenseMatrix CC(3,3);
  CC.Multiply('N','N',1.0,C,C,0.0);
  double I2 = 0.5*( C(0,0)*C(0,0) + C(1,1)*C(1,1) + C(2,2)*C(2,2) - CC(0,0) - CC(1,1) - CC(2,2));
  double J = sqrt(I3);                     // J = I3^(1/2)
  double W = c1*(I1-3.0) + c2*(I2-3.0) - m1*log(J) - m2*log(J) + 0.25*pen*(J*J -1.0 -2.0*log(J));
  */


  // ******* evaluate 2nd PK stress ********************
  double scalar1 = -2.0 * c2;     // -2 dW/dI2
  stress.Update(scalar1,C);       // S = -2 dW/dI2 C
  
  // S += 2 (dW/dI1 + I1 dW/dI2) times Identity
  double scalar2 = 2.0 * (c1 + I1*c2);
  stress.Update(scalar2,Id);

  double scalar3 = - m1 - m2 + 0.5*I3*pen - 0.5*pen;   // 2 I3 dW/dI3
  stress.Update(scalar3,Cinv);    // S += 2 I3 dW/dI3 Cinv
  // end of ******* evaluate 2nd PK stress ********************

  // ********** evaluate C-Matrix *****************************
  double delta1 = 4.0 * c2;
  double delta6 = pen * I3;
  double delta7 = 2.0*m1 + 2.0*m2 - pen*I3 + pen;
  double delta8 = - 4.0*c2;

  
  for (unsigned int i = 0; i < NUM_STRESS_3D; ++i) {
    for (unsigned int j = 0; j < NUM_STRESS_3D; ++j) {
      cmat(i,j) += delta1 * Id(i)*Id(j)         // add d1 I x I
                 + delta6 * Cinv(i)*Cinv(j);    // add d6 Cinv x Cinv
    }
  }

  AddtoCmatHolzapfelProduct(cmat,Cinv,delta7);  // add d7 Cinv o Cinv
  AddtoCmatHolzapfelProduct(cmat,Id,delta8);    // add d8 I o I
  // end of ********** evaluate C-Matrix *****************************



  return;
}


#endif
