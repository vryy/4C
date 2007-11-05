/*!----------------------------------------------------------------------
\file hyperpolyconvex.cpp
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
#include "hyperpolyconvex.H"

extern struct _MATERIAL *mat;

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)     maf 07/07|
 *----------------------------------------------------------------------*/
MAT::HyperPolyconvex::HyperPolyconvex()
  : matdata_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)      maf 07/07|
 *----------------------------------------------------------------------*/
MAT::HyperPolyconvex::HyperPolyconvex(MATERIAL* matdata)
  : matdata_(matdata)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)     maf 07/07|
 *----------------------------------------------------------------------*/
void MAT::HyperPolyconvex::Pack(vector<char>& data) const
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
void MAT::HyperPolyconvex::Unpack(const vector<char>& data)
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
double MAT::HyperPolyconvex::Density()
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

void MAT::HyperPolyconvex::Evaluate(const Epetra_SerialDenseVector* glstrain,
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
  LINALG::SymmetricEigen(Ccopy,lambda,3,'N');
  
  // evaluate principle Invariants of C
  Epetra_SerialDenseVector Inv(3);
  Inv(0) = lambda(0) + lambda(1) + lambda(2);
  Inv(1) = lambda(0) * lambda(1) + lambda(0) * lambda(2) + lambda(1) * lambda(2);
  Inv(2) = lambda(0) * lambda(1) * lambda(2);
  
  if ((gp==0) && (ele_ID == 0)) cout << "I3: " << Inv(2) << endl;
  
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
  
  
  // ************* evaluate C-matrix ***************************
  Epetra_SerialDenseMatrix I9     = tensorproduct(I,I,1.0,1.0);
  Epetra_SerialDenseMatrix IC     = tensorproduct(I,C,1.0,1.0);
  Epetra_SerialDenseMatrix CI     = tensorproduct(C,I,1.0,1.0);
  Epetra_SerialDenseMatrix ICinv  = tensorproduct(I,Cinv,1.0,1.0);
  Epetra_SerialDenseMatrix CinvI  = tensorproduct(Cinv,I,1.0,1.0);
  Epetra_SerialDenseMatrix CC     = tensorproduct(C,C,1.0,1.0);
  Epetra_SerialDenseMatrix CCinv  = tensorproduct(C,Cinv,1.0,1.0);
  Epetra_SerialDenseMatrix CinvC  = tensorproduct(Cinv,C,1.0,1.0);
  Epetra_SerialDenseMatrix CiCi   = tensorproduct(Cinv,Cinv,1.0,1.0);
  Epetra_SerialDenseMatrix HH     = tensorproduct(I,I,kappa,kappa);
  Epetra_SerialDenseMatrix HC     = tensorproduct(I,C,kappa,1.0);
  Epetra_SerialDenseMatrix CH     = tensorproduct(C,I,1.0,kappa);
  Epetra_SerialDenseMatrix HCinv  = tensorproduct(I,Cinv,kappa,1.0);
  Epetra_SerialDenseMatrix CinvH  = tensorproduct(Cinv,I,1.0,kappa);
  
  Epetra_SerialDenseMatrix CinvoCinv(9,9);
  for (int k = 0; k < 3; ++k) {
    for (int l = 0; l < 3; ++l) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          CinvoCinv(i+3*k,j+3*l) = 0.5 * (Cinv(k,i)*Cinv(l,j) + Cinv(k,j)*Cinv(l,i));
        }
      }
    }
  }
  
  Epetra_SerialDenseVector delta(7);          // deltas
  // ground substance
  delta(2) += -c * 4.0 / 3.0 * pow(Inv(2),-1.0/3.0);
  delta(5) += 4.0 / 9.0 * c * Inv(0) * pow(Inv(2),-1.0/3.0);
  delta(6) += 4.0 / 3.0 * c * Inv(0) * pow(Inv(2),-1.0/3.0);
  // penalty
  delta(5) += 4.0 * epsilon * pow(gamma,2.0) * (pow(Inv(2),gamma) + pow(Inv(2),-gamma));
  delta(6) += -4.0 * epsilon * gamma * (pow(Inv(2),gamma) - pow(Inv(2),-gamma));
  
  Epetra_SerialDenseMatrix Celast(9,9);
  for (int k = 0; k < 3; ++k) {
    for (int l = 0; l < 3; ++l) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          Celast(i+3*k,j+3*l) += delta(0) * I9(i+3*k,j+3*l) + deltafib * HH(i+3*k,j+3*l);
          Celast(i+3*k,j+3*l) += delta(1) * (IC(i+3*k,j+3*l) + CI(i+3*k,j+3*l));
          Celast(i+3*k,j+3*l) += delta(2) * (ICinv(i+3*k,j+3*l) + CinvI(i+3*k,j+3*l));
          Celast(i+3*k,j+3*l) += delta(3) * CC(i+3*k,j+3*l);
          Celast(i+3*k,j+3*l) += delta(4) * (CCinv(i+3*k,j+3*l) + CinvC(i+3*k,j+3*l));
          Celast(i+3*k,j+3*l) += delta(5) * CiCi(i+3*k,j+3*l);
          Celast(i+3*k,j+3*l) += delta(6) * CinvoCinv(i+3*k,j+3*l);
          Celast(i+3*k,j+3*l) += delta(7);
        }
      }
    }
  }
  (*cmat)(0,0)=Celast(0,0);
  (*cmat)(0,1)=Celast(1,1);
  (*cmat)(0,2)=Celast(2,2);
  (*cmat)(0,3)=Celast(1,0);
  (*cmat)(0,4)=Celast(2,1);
  (*cmat)(0,5)=Celast(2,0);
  
  (*cmat)(1,0)=Celast(3,3);
  (*cmat)(1,1)=Celast(4,4);
  (*cmat)(1,2)=Celast(5,5);
  (*cmat)(1,3)=Celast(4,3);
  (*cmat)(1,4)=Celast(5,4);
  (*cmat)(1,5)=Celast(5,3);
  
  (*cmat)(2,0)=Celast(6,6);
  (*cmat)(2,1)=Celast(7,7);
  (*cmat)(2,2)=Celast(8,8);
  (*cmat)(2,3)=Celast(7,6);
  (*cmat)(2,4)=Celast(8,7);
  (*cmat)(2,5)=Celast(8,6);

  (*cmat)(3,0)=Celast(3,0);
  (*cmat)(3,1)=Celast(4,1);
  (*cmat)(3,2)=Celast(5,2);
  (*cmat)(3,3)=Celast(4,0);
  (*cmat)(3,4)=Celast(5,1);
  (*cmat)(3,5)=Celast(5,0);
  
  (*cmat)(4,0)=Celast(6,3);
  (*cmat)(4,1)=Celast(7,4);
  (*cmat)(4,2)=Celast(8,5);
  (*cmat)(4,3)=Celast(7,3);
  (*cmat)(4,4)=Celast(8,4);
  (*cmat)(4,5)=Celast(8,3);
  
  (*cmat)(5,0)=Celast(6,0);
  (*cmat)(5,1)=Celast(7,1);
  (*cmat)(5,2)=Celast(8,2);
  (*cmat)(5,3)=Celast(7,0);
  (*cmat)(5,4)=Celast(8,1);
  (*cmat)(5,5)=Celast(8,0);
  
  return;
}

/*----------------------------------------------------------------------*
 |  Calculate tensor product of two tensors of 2nd order       maf 04/07|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix MAT::HyperPolyconvex::tensorproduct(
                                             const Epetra_SerialDenseMatrix A,
                                             const Epetra_SerialDenseMatrix B,
                                             const double scalarA,
                                             const double scalarB)
{
  Epetra_SerialDenseMatrix AB(9,9);
  for (int k = 0; k < 3; ++k) {
    for (int l = 0; l < 3; ++l) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          AB(i+3*k,j+3*l) = scalarA * A(k,l) * scalarB * B(i,j); 
        }
      }
    }
  }
  return AB;
}


#endif
