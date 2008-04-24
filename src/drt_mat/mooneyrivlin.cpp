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

extern struct _MATERIAL *mat;


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

*/

void MAT::MooneyRivlin::Evaluate(const Epetra_SerialDenseVector* glstrain,
                                       Epetra_SerialDenseMatrix* cmat,
                                       Epetra_SerialDenseVector* stress)

{
  // get material parameters
  double m1  = matdata_->m.mooneyrivlin->mu1;
  double m2  = matdata_->m.mooneyrivlin->mu2;
  double c1  = m1 / matdata_->m.mooneyrivlin->alpha1;  // c1 = mu1/a1
  double c2  = m2 / matdata_->m.mooneyrivlin->alpha2;  // c2 = mu2/a2
  double pen = matdata_->m.mooneyrivlin->penalty;

  // Identity Matrix
  Epetra_SerialDenseMatrix I(3,3);
  for (int i = 0; i < 3; ++i) I(i,i) = 1.0;

  // Green-Lagrange Strain Tensor
  Epetra_SerialDenseMatrix E(3,3);
  E(0,0) = (*glstrain)(0);
  E(1,1) = (*glstrain)(1);
  E(2,2) = (*glstrain)(2);
  E(0,1) = 0.5 * (*glstrain)(3);  E(1,0) = 0.5 * (*glstrain)(3);
  E(1,2) = 0.5 * (*glstrain)(4);  E(2,1) = 0.5 * (*glstrain)(4);
  E(0,2) = 0.5 * (*glstrain)(5);  E(2,0) = 0.5 * (*glstrain)(5);

  // Right Cauchy-Green Tensor  C = 2 * E + I
  Epetra_SerialDenseMatrix C(E);
  C.Scale(2.0);
  C += I;
  
  // compute Invariants
  double I1 = C(0,0) + C(1,1) + C(2,2);    // I1 = trace(C)
  
  // compute determinant of C by Sarrus' rule
  double I3  = C(0,0) * C(1,1) * C(2,2)    // I3 = det(C)
             + C(0,1) * C(1,2) * C(2,0)
             + C(0,2) * C(1,0) * C(2,1)
             - C(0,0) * C(1,2) * C(2,1)
             - C(0,1) * C(1,0) * C(2,2)
             - C(0,2) * C(1,1) * C(2,0);
  
  /* //Compute strain-energy function W
  // I2 = 1/2 (trace(C)^2 - trace(CxC))
  Epetra_SerialDenseMatrix CC(3,3);
  CC.Multiply('N','N',1.0,C,C,0.0);
  double I2 = 0.5*( C(0,0)*C(0,0) + C(1,1)*C(1,1) + C(2,2)*C(2,2) - CC(0,0) - CC(1,1) - CC(2,2));
  double J = sqrt(I3);                     // J = I3^(1/2)
  double W = c1*(I1-3.0) + c2*(I2-3.0) - m1*log(J) - m2*log(J) + 0.25*pen*(J*J -1.0 -2.0*log(J));
  */

  // compute C^-1
  Epetra_SerialDenseMatrix Cinv(C);
  Epetra_SerialDenseSolver solve_for_inverseC;
  solve_for_inverseC.SetMatrix(Cinv);
  int err2 = solve_for_inverseC.Factor();
  int err = solve_for_inverseC.Invert();
  if ((err != 0) && (err2!=0)) dserror("Inversion of Cauchy-Green failed");

  // ******* evaluate 2nd PK stress ********************
  Epetra_SerialDenseMatrix S(C);    // S = C
  double scalar1 = -2.0 * c2;       // -2 dW/dI2
  S.Scale(scalar1);
  
  // S += 2 (dW/dI1 + I1 dW/dI2) times Identity 
  double scalar2 = 2.0 * (c1 + I1*c2);
  for (int i = 0; i < 3; ++i) S(i,i) += scalar2;
  
  Epetra_SerialDenseMatrix Cinvscale(Cinv);
  double scalar3 = - m1 - m2 + 0.5*I3*pen - 0.5*pen;   // 2 I3 dW/dI3
  Cinvscale.Scale(scalar3);
  S += Cinvscale;
  
  (*stress)(0) = S(0,0);
  (*stress)(1) = S(1,1);
  (*stress)(2) = S(2,2);
  (*stress)(3) = S(0,1);
  (*stress)(4) = S(1,2);
  (*stress)(5) = S(0,2);
  // end of ******* evaluate 2nd PK stress ********************
  
  // ********** evaluate C-Matrix *****************************
  double delta1 = 4.0 * c2;
  double delta6 = pen * I3;
  double delta7 = 2.0*m1 + 2.0*m2 - pen*I3 + pen;
  double delta8 = - 4.0*c2;
  
  ElastSymTensorMultiply((*cmat),delta1,I,I,0.0);  // initialize cmat and add d1 I x I
  ElastSymTensorMultiply((*cmat),delta6,Cinv,Cinv,1.0);    // add d6 Cinv x Cinv
  ElastSymTensor_o_Multiply((*cmat),delta7,Cinv,Cinv,1.0); // add d7 Cinv o Cinv
  ElastSymTensor_o_Multiply((*cmat),delta8,I,I,1.0);       // add d8 I o I
  // end of ********** evaluate C-Matrix *****************************

  
  
  return;
}


#endif
