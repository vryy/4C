/*!----------------------------------------------------------------------
\file anisotropic_balzani.cpp
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
#include "../drt_lib/linalg_utils.H"
#include "anisotropic_balzani.H"

extern struct _MATERIAL *mat;

using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)     maf 07/07|
 *----------------------------------------------------------------------*/
MAT::AnisotropicBalzani::AnisotropicBalzani()
  : matdata_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)      maf 07/07|
 *----------------------------------------------------------------------*/
MAT::AnisotropicBalzani::AnisotropicBalzani(MATERIAL* matdata)
  : matdata_(matdata)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)     maf 07/07|
 *----------------------------------------------------------------------*/
void MAT::AnisotropicBalzani::Pack(vector<char>& data) const
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
 |  Unpack                                        (public)     maf 07/07|
 *----------------------------------------------------------------------*/
void MAT::AnisotropicBalzani::Unpack(const vector<char>& data)
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
double MAT::AnisotropicBalzani::Density()
{
  return matdata_->m.anisotropic_balzani->density;  // density, returned to evaluate mass matrix
}


/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)     maf 07/07|
 *----------------------------------------------------------------------*

This routine establishes a local material law, stress-strain relationship
for hyperelastic, anisotropic material for a 3D-hex-element. 
Used for biological, soft, collagenous tissues.

Documented in the 'Diplomarbeit' by Barbara Roehrnbauer at LNM.
Based on Holzapfel [1], Ogden [2] and Balzani, Schroeder, Neff [3].

[1] G.A.Holzapfel, R.W.Ogden, A New Consitutive Framework for Arterial Wall Mechanics and
  a Comparative Study of Material Models, Journal of Elasticity 61, 1-48, 2000. 
[2] R.W.Ogden, Anisotropy and Nonlinear Elasticity in Arterial Wall Mechanics,  
  CISM Course on Biomechanical Modeling, Lectures 2,3, 2006.
[3] D.Balzani, P.Neff, J.Schroeder, G.A.Holzapfel, A Polyconvex Framework for Soft Biological Tissues
  Adjustment to Experimental Data, Report-Preprint No. 22, 2005.
*/

void MAT::AnisotropicBalzani::Evaluate(const Epetra_SerialDenseVector* glstrain,
                                    const Epetra_SerialDenseMatrix* defgrd,
                                    const int gp, const int ele_ID, const double time,
                                      Epetra_SerialDenseMatrix* cmat,
                                      Epetra_SerialDenseVector* stress)
{
  // get material parameters
  double c1 = matdata_->m.anisotropic_balzani->c1;          //parameter for ground substance
  double eps1 = matdata_->m.anisotropic_balzani->eps1;      //parameter for fiber potential
  double eps2 = matdata_->m.anisotropic_balzani->eps2;      //parameter for fiber potential
  double alpha1 = matdata_->m.anisotropic_balzani->alpha1;  //parameter for fiber potential
  double alpha2 = matdata_->m.anisotropic_balzani->alpha2;  //parameter for fiber potential
    
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
  
//  cout << "Enhanced Cauchy-Green: " << C;
//  
//  Epetra_SerialDenseMatrix C_lock(3,3);
//  C_lock.Multiply('T','N',1.0,*defgrd,*defgrd,0.0);
//  cout << "Disp-based Cauchy-Green: " << C_lock;
  
  // compute eigenvalues of C
  Epetra_SerialDenseMatrix Ccopy(C);
  Epetra_SerialDenseVector lambda(3);
  SymmetricEigen(Ccopy,lambda,3,'N');
  
  // evaluate principle Invariants of C
  //Epetra_SerialDenseVector Inv(3);
  double I1 = lambda(0) + lambda(1) + lambda(2);
  //double I2 = lambda(0) * lambda(1) + lambda(0) * lambda(2) + lambda(1) * lambda(2);
  double I3 = lambda(0) * lambda(1) * lambda(2);  // = det(C)
  
  // compute C^-1
  Epetra_SerialDenseMatrix Cinv(C);
  Epetra_SerialDenseSolver solve_for_inverseC;
  solve_for_inverseC.SetMatrix(Cinv);
  int err2 = solve_for_inverseC.Factor();        
  int err = solve_for_inverseC.Invert();
  if ((err != 0) && (err2!=0)) dserror("Inversion of Cauchy-Green failed");
  // Cinv = C^-T
  Cinv.SetUseTranspose(true);
  
  // Structural Tensor M, defined by a x a
  Epetra_SerialDenseVector a(3);
  a(0) = matdata_->m.anisotropic_balzani->a1[0];  // first fiber vector from input
  a(1) = matdata_->m.anisotropic_balzani->a1[1];
  a(2) = matdata_->m.anisotropic_balzani->a1[2];
  Epetra_SerialDenseMatrix M(3,3);
  M.Multiply('N','T',1.0,a,a,0.0);
  
  Epetra_SerialDenseMatrix CM(3,3);
  CM.Multiply('N','N',1.0,C,M,0.0);
  // compute CM + MC
  Epetra_SerialDenseMatrix CMMC(CM);  // CMMC = CM
  CMMC.Multiply('N','N',1.0,M,C,1.0); // CMMC = CM + MC
  // anisotropic Invariant J4 = tr[CM]
  double J4 = CM(0,0) + CM(1,1) + CM(2,2);

  Epetra_SerialDenseMatrix Csq(3,3);
  Csq.Multiply('N','N',1.0,C,C,0.0);
  Epetra_SerialDenseMatrix CsqM(3,3);
  CsqM.Multiply('N','N',1.0,Csq,M,0.0);
  // anisotropic Invariant J5 = tr[C^2M]
  double J5 = CsqM(0,0) + CsqM(1,1) + CsqM(2,2);
  
  // anisotropic associated reducible Invariant K3 = I1*J4 - J5
  double K3 = I1 * J4 - J5;
  double fiberscalar = alpha1*pow(K3-2.0,alpha2)*alpha2 / (K3-2.0);
  
  // cofC = det(C)*C^-T
  Epetra_SerialDenseMatrix CofC(Cinv);
  CofC.Scale(I3);
  
  /* Underlying strain-energy function
   * Wiso = c1*(\frac{I_1}{I_3^{\frac{1}{3}} -3)                  // ground substance SEF
   *      + eps1 * (I_3^eps2 + \frac{1}{I_3^eps2} - 2)            // penalty function
   *      + alpha1*(I_1*J4 - J5 - 2)^alpha2 if ... >= 2; else 0   // fiber SEF
   */
    
  // ******* evaluate 2nd PK stress ********************
  Epetra_SerialDenseMatrix S(CofC);   // S = C^{-T}
  double scalar1 = -1.0/3.0 * c1*I1 * pow(I3,-4.0/3.0)
                 + eps1 * ( (pow(I3,eps2)*eps2)/I3 - eps2/( pow(I3,eps2)*I3));
  S.Scale(scalar1);
  double scalar2 = c1*pow(I3,-1.0/3.0);
  Epetra_SerialDenseMatrix Iscale(I);
  Iscale.Scale(scalar2);
  S += Iscale;
  
  // ***** fiber part *******
  if (K3 >= 2.0){
    Epetra_SerialDenseMatrix Mscale(M);
    Mscale.Scale(fiberscalar*I1);
    S += Mscale;
    // compute CM + MC
    Epetra_SerialDenseMatrix CMMCscale(CMMC);
    CMMCscale.Scale(-fiberscalar);
    S += CMMCscale;
    
    Iscale = I;
    Iscale.Scale(fiberscalar*J4);
    S += Iscale;
  }
  
  // the wellknown factor 2!
  S.Scale(2.0);
  
  (*stress)(0) = S(0,0);
  (*stress)(1) = S(1,1);
  (*stress)(2) = S(2,2);
  (*stress)(3) = S(0,1);
  (*stress)(4) = S(1,2);
  (*stress)(5) = S(0,2);
  // end of ******* evaluate 2nd PK stress ********************
  
  // ********** evaluate C-Matrix *****************************
  // isotropic part
  double d2W_dI3dI3 = 4.0/9.0 * c1 * I1 * pow(I3,-7.0/3.0)
                    + eps1 * ( pow(I3,eps2) * eps2*eps2 / (I3*I3)
                             - pow(I3,eps2) * eps2 / (I3*I3)
                             + eps2*eps2 / (pow(I3,eps2) * I3*I3)
                             + eps2 / (pow(I3,eps2) * I3*I3) );
  double d2W_dI3dI1 = -1.0/3.0 * c1 * pow(I3,-4.0/3.0);
  double dW_dI3     = -1.0/3.0 * c1 * I1 * pow(I3,-4.0/3.0)
                      + eps1 * (pow(I3,eps2)*eps2/I3 - eps2/(pow(I3,eps2)*I3));
  
  ElastSymTensorMultiply((*cmat),d2W_dI3dI3,CofC,CofC,1.0);     // CofC x CofC
  ElastSymTensorMultiplyAddSym((*cmat),d2W_dI3dI1,I,Cinv,1.0);  // I x Cinv + Cinv x I
  ElastSymTensorMultiply((*cmat),I3 * dW_dI3,Cinv,Cinv,1.0);    // Cinv x Cinv
  ElastSymTensor_o_Multiply((*cmat),-I3 * dW_dI3,Cinv,Cinv,1.0);// - Cinv o Cinv
  
  // fiber part
  if (K3 >= 2.0){
    // compute derivatives of W_ti w.r.t. Invariants
    double K3m2sq = (K3 - 2.0) * (K3 - 2.0);
    double K3m2p = pow((K3 - 2.0),alpha2);
    double d2W_dJ4dJ4 = ( alpha1 * K3m2p * alpha2*alpha2 * I1*I1 
                         -alpha1 * K3m2p * alpha2 * I1*I1 ) / K3m2sq;
    double d2W_dJ5dJ5 = ( alpha1 * K3m2p * alpha2*alpha2  
                         -alpha1 * K3m2p * alpha2 ) / K3m2sq;
    double d2W_dI1dJ4 = ( alpha1 * K3m2p * alpha2*alpha2 * J4 * I1  
                         -alpha1 * K3m2p * alpha2 * J4 * I1 ) / K3m2sq
                         +alpha1 * K3m2p * alpha2 / (K3 - 2.0);
    double d2W_dI1dJ5 = (-alpha1 * K3m2p * alpha2*alpha2 * J4
                         +alpha1 * K3m2p * alpha2 * J4 ) / K3m2sq;
    double d2W_dJ4dJ5 = (-alpha1 * K3m2p * alpha2*alpha2 * I1
                         +alpha1 * K3m2p * alpha2 * I1 ) / K3m2sq;
    double dW_dJ5     = - alpha1 * K3m2p * alpha2 / (K3 - 2.0);
    // multiply these with corresponding tensor products
    ElastSymTensorMultiply((*cmat),d2W_dJ4dJ4,M,M,1.0);         // M x M
    ElastSymTensorMultiply((*cmat),d2W_dJ5dJ5,CMMC,CMMC,1.0);   // CM+MC x CM+MC
    ElastSymTensorMultiplyAddSym((*cmat),d2W_dI1dJ4,I,M,1.0);   // I x M + M x I
    ElastSymTensorMultiplyAddSym((*cmat),d2W_dI1dJ5,CMMC,I,1.0);// (CM+MC) x I + I x (CM+MC)
    ElastSymTensorMultiplyAddSym((*cmat),d2W_dJ4dJ5,CMMC,M,1.0);// (CM+MC) x M + M x (CM+MC)
    ElastSymTensor_o_Multiply((*cmat),dW_dJ5,I,M,1.0);          // I o M
    ElastSymTensor_o_Multiply((*cmat),dW_dJ5,M,I,1.0);          // M o I
  }
  
  // the factor 4!
  (*cmat).Scale(4.0);
  // end of ********** evaluate C-Matrix *****************************
  
  return;
}


#endif
