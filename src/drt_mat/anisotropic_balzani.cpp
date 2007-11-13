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
  return matdata_->m.hyper_polyconvex->density;  // density, returned to evaluate mass matrix
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
  double c = matdata_->m.hyper_polyconvex->c;             //parameter for ground substance
  double k1 = matdata_->m.hyper_polyconvex->k1;           //parameter for fiber potential
  double k2 = matdata_->m.hyper_polyconvex->k2;           //parameter for fiber potential
  double gamma = matdata_->m.hyper_polyconvex->gamma;     //penalty parameter
  double epsilon = matdata_->m.hyper_polyconvex->epsilon; //penalty parameter
  
  double kappa = 1.0/3.0;   //Dispersions Parameter
  //double kappa = 1.0;   //Dispersions Parameter
  //double phi = 0.0;         //Angle for Anisotropic Fiber Orientation
  //double theta = 0.0;       //Angle for Anisotropic Fiber Orientation
  
  // Vector of Preferred Direction
  Epetra_SerialDenseVector ad(3);
  ad(0) = 1.0;
  
  // Orientation Tensor
  Epetra_SerialDenseMatrix M(3,3);
  M.Multiply('N','T',1.0,ad,ad,0.0);
  
  // Identity Matrix
  Epetra_SerialDenseMatrix I(3,3);
  for (int i = 0; i < 3; ++i) I(i,i) = 1.0;
  
//  // Structural Tensor
//  Epetra_SerialDenseMatrix H(I);
//  H.Scale(kappa);
  
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
  Epetra_SerialDenseVector Inv(3);
  Inv(0) = lambda(0) + lambda(1) + lambda(2);
  Inv(1) = lambda(0) * lambda(1) + lambda(0) * lambda(2) + lambda(1) * lambda(2);
  Inv(2) = lambda(0) * lambda(1) * lambda(2);
  
  //if ((gp==0) && (ele_ID == 0)) cout;// << "I3: " << Inv(2) << endl;
  
  // compute C^-1
  Epetra_SerialDenseMatrix Cinv(C);
  Epetra_SerialDenseSolver solve_for_inverseC;
  solve_for_inverseC.SetMatrix(Cinv);
  int err2 = solve_for_inverseC.Factor();        
  int err = solve_for_inverseC.Invert();
  if ((err != 0) && (err2!=0)) dserror("Inversion of Cauchy-Green failed");
  Cinv.SetUseTranspose(true);
  
  // Structural Tensor H defined implicitly as H=kappa*I
  Epetra_SerialDenseMatrix HxC(3,3);
  HxC.Multiply('N','N',kappa,I,C,0.0);
  
  // Anisotropic Invariant K
  double K = HxC(0,0) + HxC(1,1) + HxC(2,2);
  
  // ******* evaluate 2nd PK stress ********************
  Epetra_SerialDenseMatrix S(Cinv);   // S = C^{-T}
  
  // Ground Substance
  double scalar = -1.0/3.0 * Inv(0);
  S.Scale(scalar);                    // S = -1/3 I_1 * C^{-T}
  S += I;                             // S = -1/3 I_1 * C^{-T} + I
  scalar = 2.0 * c * pow(Inv(2),-1.0/3.0);
  S.Scale(scalar);                    // S = 2cI_3^{-1/3} * (-1/3 I_1 * C^{-T} + I) = S_GS
  
  // Penalty
  scalar = 2.0 * epsilon * gamma * (pow(Inv(2),gamma) - pow(Inv(2),-gamma));
  Epetra_SerialDenseMatrix S_pen(Cinv);// S_pen = C^{-T}
  S_pen.Scale(scalar);                 // S_pen = 2*eps*gam*(I_3^gam - I_3^{-gam}) * C^{-T}
  S += S_pen;                          // S = S_GS + S_pen
  
  // Fiber
  double deltafib = 0.0;
  if (K >= 1.0){
    scalar = 2.0 * k1 * exp(k2 * pow( K-1.0 ,2.0)) * (K-1.0);
    Epetra_SerialDenseMatrix S_fiber(I);  // scalar_fib = 2k1 e^{k2(K-1)^2} (K-1)
    S_fiber.Scale(scalar * kappa);        // S_fiber = scalar_fib * (kappa*I)
    S += S_fiber;
    deltafib = 4.0 * k1 * exp(k2 * pow( K-1.0 ,2.0)) * ((2.0 * k2 * pow( K-1.0 ,2)) + 1.0);
  }
  
  (*stress)(0) = S(0,0);
  (*stress)(1) = S(1,1);
  (*stress)(2) = S(2,2);
  (*stress)(3) = S(0,1);
  (*stress)(4) = S(1,2);
  (*stress)(5) = S(0,2);
  // end of ******* evaluate 2nd PK stress ********************
  
  Epetra_SerialDenseVector delta(7);          // deltas
  // ground substance
  delta(2) += -c * 4.0 / 3.0 * pow(Inv(2),-1.0/3.0);
  delta(5) += 4.0 / 9.0 * c * Inv(0) * pow(Inv(2),-1.0/3.0);
  delta(6) += 4.0 / 3.0 * c * Inv(0) * pow(Inv(2),-1.0/3.0);
  // penalty
  delta(5) += 4.0 * epsilon * pow(gamma,2.0) * (pow(Inv(2),gamma) + pow(Inv(2),-gamma));
  delta(6) += -4.0 * epsilon * gamma * (pow(Inv(2),gamma) - pow(Inv(2),-gamma));
  
  // *** new "faster" evaluate of C-Matrix
  ElastSymTensorMultiply((*cmat),delta(0),I,I,1.0);           // I x I
  ElastSymTensorMultiplyAddSym((*cmat),delta(1),I,C,1.0);     // I x C + C x I
  ElastSymTensorMultiplyAddSym((*cmat),delta(2),I,Cinv,1.0);  // I x Cinv + Cinv x I
  ElastSymTensorMultiply((*cmat),delta(3),C,C,1.0);           // C x C
  ElastSymTensorMultiplyAddSym((*cmat),delta(4),C,Cinv,1.0);  // C x Cinv + Cinv x C
  ElastSymTensorMultiply((*cmat),delta(5),Cinv,Cinv,1.0);     // Cinv x Cinv
  ElastSymTensor_o_Multiply((*cmat),delta(6),Cinv,Cinv,1.0);  // Cinv o Cinv
  // fiber part
  ElastSymTensorMultiply((*cmat),(deltafib*kappa*kappa),I,I,1.0);
  
  
  
  return;
}


#endif
