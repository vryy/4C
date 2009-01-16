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
 *----------------------------------------------------------------------*/
void MAT::MooneyRivlin::Evaluate(
        const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
        LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
        LINALG::Matrix<NUM_STRESS_3D,1> * stress)
{
  // get material parameters
  const double c1  = matdata_->m.mooneyrivlin->c1;
  const double c2  = matdata_->m.mooneyrivlin->c2;
  const double kappa_q1 = matdata_->m.mooneyrivlin->kap; // kappa_q1*(J-1)^2
  const double lambda = matdata_->m.mooneyrivlin->lambda;
  
  // penalty param for Klinkel-formulation 
  // Watch out: It is not stress free in reference configuration!
  double kappa_q2 = 0.;  
  double kappa_ln = 0.; 

  if (lambda !=0.0)
  {
    kappa_q2 = lambda/4.;  // kappa_q2*(J^2-1)
    kappa_ln = lambda/2.-6.*c2;  //kappa_ln*ln(J)
  }
  
  // aux param d (to ensure stress free reference configuration -> Holzapfel)
  const double d = 2.*c1+4.*c2+kappa_ln;
  
  // right Cauchy-Green Tensor  C = 2 * E + I
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  for (int i = 0; i < 3; i++) C(i) = 2.*C(i);  // respect factor 2 in shear terms of glstrain
  C += Id;

  // invariants
  const double I1 = C(0) + C(1) + C(2);  // 1st invariant, trace
  const double I3 = C(0)*C(1)*C(2)
        + C(3)*C(4)*C(5)
        - C(1)*C(5)*C(5)
        - C(2)*C(3)*C(3)
        - C(0)*C(4)*C(4);    // 3rd invariant, determinant

  const double J = sqrt(I3);
  
  // invert C
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv(false);

  Cinv(0) = C(1)*C(2) - C(4)*C(4);
  Cinv(1) = C(0)*C(2) - C(5)*C(5);
  Cinv(2) = C(0)*C(1) - C(3)*C(3);
  Cinv(3) = C(5)*C(4) - C(3)*C(2);
  Cinv(4) = C(3)*C(5) - C(0)*C(4);
  Cinv(5) = C(3)*C(4) - C(5)*C(1);
  Cinv.Scale(1.0/I3);

  /* //Compute strain-energy function W
  // I2 = 1/2 (trace(C)^2 - trace(CxC))
  Epetra_SerialDenseMatrix CC(3,3);
  CC.Multiply('N','N',1.0,C,C,0.0);
  double I2 = 0.5*( C(0,0)*C(0,0) + C(1,1)*C(1,1) + C(2,2)*C(2,2) - CC(0,0) - CC(1,1) - CC(2,2));
  double J = sqrt(I3);                     // J = I3^(1/2)
  double W = c1*(I1-3.0) + c2*(I2-3.0) - d*log(J) + pen*pow(J -1.0,2)
  */


  // ******* evaluate 2nd PK stress ********************
  // gammas from Holzapfel page 248
  
  const double gamma1 = 2.0 * (c1 + I1*c2); //2 (dW/dI1 + I1 dW/dI2)
  (*stress).Update(gamma1,Id,0.0); //S +=  gamma1 times Identity

  const double gamma2 = -2.0 * c2;     // -2 dW/dI2
  (*stress).Update(gamma2,C,1.0);       // S = gamma2 times C
    
  const double gamma3 = 2.*kappa_q1*J*(J-1.0)-d+2.*kappa_q2*J*J;   // 2 I3 dW/dI3
  (*stress).Update(gamma3,Cinv,1.0);    // S += gamma3 times Cinv
  // end of ******* evaluate 2nd PK stress ********************

  // ********** evaluate C-Matrix *****************************
  // deltas from Holzapfel page 262
  const double delta1 = 4.0 * c2;
  const double delta6 = 2.*kappa_q1*J*(2.*J-1.) + 4.*kappa_q2*J*J;
  const double delta7 = -2.*(2.*kappa_q1*J*(J-1.)-d) - 4.*kappa_q2*J*J;
  const double delta8 = - 4.0*c2;

  for (unsigned int i = 0; i < NUM_STRESS_3D; ++i) {
    for (unsigned int j = 0; j < NUM_STRESS_3D; ++j) {
      (*cmat)(i,j) += delta1 * Id(i)*Id(j)         // add d1 I x I
                 + delta6 * Cinv(i)*Cinv(j);    // add d6 Cinv x Cinv
    }
  }

  AddtoCmatHolzapfelProduct((*cmat),Cinv,delta7);  // add d7 Cinv o Cinv
  AddtoCmatHolzapfelProduct((*cmat),Id,delta8);    // add d8 I o I
  
  // end of ********** evaluate C-Matrix *****************************



  return;
}


#endif
