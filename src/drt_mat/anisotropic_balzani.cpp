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


#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "Epetra_SerialDenseSolver.h"
#include "anisotropic_balzani.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::AnisotropicBalzani::AnisotropicBalzani(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c1_(matdata->GetDouble("C1")),
  eps1_(matdata->GetDouble("EPS1")),
  eps2_(matdata->GetDouble("EPS2")),
  alpha1_(matdata->GetDouble("ALPHA1")),
  alpha2_(matdata->GetDouble("ALPHA2")),
  dens_(matdata->GetDouble("DENS")),
  aloc_(matdata->GetInt("ALOC")),
  a1_(),
  alpha1_2_(matdata->GetDouble("ALPHA1_2")),
  alpha2_2_(matdata->GetDouble("ALPHA2_2")),
  a2_()
{
  a1_[0] = matdata->GetDouble("A1X");
  a1_[1] = matdata->GetDouble("A1Y");
  a1_[2] = matdata->GetDouble("A1Z");

  a2_[0] = matdata->GetDouble("A2X");
  a2_[1] = matdata->GetDouble("A2Y");
  a2_[2] = matdata->GetDouble("A2Z");
}

Teuchos::RCP<MAT::Material> MAT::PAR::AnisotropicBalzani::CreateMaterial()
{
  return Teuchos::rcp(new MAT::AnisotropicBalzani(this));
}

MAT::AnisotropicBalzaniType MAT::AnisotropicBalzaniType::instance_;


DRT::ParObject* MAT::AnisotropicBalzaniType::Create( const std::vector<char> & data )
{
  MAT::AnisotropicBalzani* anba = new MAT::AnisotropicBalzani();
  anba->Unpack(data);
  return anba;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)     maf 07/07|
 *----------------------------------------------------------------------*/
MAT::AnisotropicBalzani::AnisotropicBalzani()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                  (public)      maf 07/07|
 *----------------------------------------------------------------------*/
MAT::AnisotropicBalzani::AnisotropicBalzani(MAT::PAR::AnisotropicBalzani* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)     maf 07/07|
 *----------------------------------------------------------------------*/
void MAT::AnisotropicBalzani::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
  // fibers
  AddtoPack(data,a1_);  // fiber vector 1
  AddtoPack(data,a2_);  // fiber vector 2
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)     maf 07/07|
 *----------------------------------------------------------------------*/
void MAT::AnisotropicBalzani::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::AnisotropicBalzani*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // fibres
  ExtractfromPack(position,data,a1_);  // fiber vector 1
  ExtractfromPack(position,data,a2_);  // fiber vector 2

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

void MAT::AnisotropicBalzani::Setup(DRT::INPUT::LineDefinition* linedef)
{
  // check whether fiber are based on local cosy
  if (params_->aloc_ == 1){
    // fibers aligned in local element cosy with gamma_i around circumferential direction
    vector<double> rad;
    vector<double> axi;
    vector<double> cir;
    // read local (cylindrical) cosy-directions at current element
    linedef->ExtractDoubleVector("RAD",rad);
    linedef->ExtractDoubleVector("AXI",axi);
    linedef->ExtractDoubleVector("CIR",cir);
    Epetra_SerialDenseMatrix locsys(3,3);
    // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
    for (int i=0; i<3; ++i){
      locsys(i,0) = rad[i];
      locsys(i,1) = axi[i];
      locsys(i,2) = cir[i];
    }

    // alignment angles gamma_i are read from first entry of then unnecessary vectors a1 and a2
    double gamma1 = params_->a1_[0];
    gamma1 = (gamma1 * M_PI)/180.0;  // convert to radians
    double gamma2 = params_->a2_[0];
    gamma2 = (gamma2 * M_PI)/180.0;  // convert to radians

    a1_.resize(3);
    a2_.resize(3);
    vector<double> a(3);
    for (int i = 0; i < 3; ++i) {
      // a = cos gamma e1 +- sin gamma e2
      a1_.at(i) = cos(gamma1)*locsys(i,2) + sin(gamma1)*locsys(i,1);
      a2_.at(i) = cos(gamma2)*locsys(i,2) - sin(gamma2)*locsys(i,1);
    }
  } else if (params_->aloc_ == 0){
    a1_.resize(3);
    a2_.resize(3);
    for (int i = 0; i < 3; ++i) {
      a1_.at(i) = params_->a1_[i];
      a2_.at(i) = params_->a2_[i];
    }

  }
}



/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)     maf 07/07|
 *----------------------------------------------------------------------*

This routine establishes a local material law, stress-strain relationship
for hyperelastic, anisotropic material for a 3D-hex-element.
Used for biological, soft, collagenous tissues. 2 Fiber families account for
anisotropy where fiber angles need to be specified

Based on Dissertation of Daniel Balzani
*/

void MAT::AnisotropicBalzani::Evaluate(
        const LINALG::Matrix<NUM_STRESS_3D,1> * glstrain,
        const int gp, const int ele_ID, const double time,
        LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
        LINALG::Matrix<NUM_STRESS_3D,1> * stress)
{
  // wrapper for FixedSizeMatrix
  Epetra_SerialDenseMatrix cmat_e(View,cmat->A(),cmat->Rows(),cmat->Rows(),cmat->Columns());
  // stress and glstrain are copied value by value and are thus not necessary

  // get material parameters
  double c1 = params_->c1_;          //parameter for ground substance
  double eps1 = params_->eps1_;      //parameter for incomp. penalty
  double eps2 = params_->eps2_;      //parameter for incomp. penalty
  double alpha1 = params_->alpha1_;  //parameter for 1st fiber potential
  double alpha2 = params_->alpha2_;  //parameter for 1st fiber potential
  double alpha1_2 = params_->alpha1_2_;//parameter for 2nd fiber potential
  double alpha2_2 = params_->alpha2_2_;//parameter for 2nd fiber potential

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
//  //C_lock.Multiply('T','N',1.0,*defgrd,*defgrd,0.0);
//  cout << "Disp-based Cauchy-Green: " << C_lock;

//  // compute eigenvalues of C
//  Epetra_SerialDenseMatrix Ccopy(C);
//  Epetra_SerialDenseVector lambda(3);
//  SymmetricEigenValues(Ccopy,lambda);
//
//  // evaluate principle Invariants of C
//  //Epetra_SerialDenseVector Inv(3);
//  double I1 = lambda(0) + lambda(1) + lambda(2);
//  //double I2 = lambda(0) * lambda(1) + lambda(0) * lambda(2) + lambda(1) * lambda(2);
//  double I3 = lambda(0) * lambda(1) * lambda(2);  // = det(C)


  double I1 = C(0,0) + C(1,1) + C(2,2);  // I1 = trace(C)
  // compute determinant of C by Sarrus' rule
  double I3  = C(0,0) * C(1,1) * C(2,2)    // I3 = det(C)
             + C(0,1) * C(1,2) * C(2,0)
             + C(0,2) * C(1,0) * C(2,1)
             - C(0,0) * C(1,2) * C(2,1)
             - C(0,1) * C(1,0) * C(2,2)
             - C(0,2) * C(1,1) * C(2,0);

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
  if (params_->aloc_ != 1){
    a(0) = params_->a1_[0];  // first fiber vector from input
    a(1) = params_->a1_[1];
    a(2) = params_->a1_[2];
  } else {
    a(0) = a1_.at(0); a(1) = a1_.at(1); a(2) = a1_.at(2);
  }
  // normalize a
  double norma = a.Norm2();
  a.Scale(1.0/norma);
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

  /* Second fiber part ******************************************************/
  // Structural Tensor M_2, defined by a_2 x a_2
  Epetra_SerialDenseVector a_2(3);
  if (params_->aloc_ != 1){
    a_2(0) = params_->a2_[0];  // 2nd fiber vector from input
    a_2(1) = params_->a2_[1];
    a_2(2) = params_->a2_[2];
  } else {
    a_2(0) = a2_.at(0); a_2(1) = a2_.at(1); a_2(2) = a2_.at(2);
  }
  // normalize a_2
  double norma_2 = a_2.Norm2();
  a_2.Scale(1.0/norma_2);
  Epetra_SerialDenseMatrix M_2(3,3);
  M_2.Multiply('N','T',1.0,a_2,a_2,0.0);

  Epetra_SerialDenseMatrix CM_2(3,3);
  CM_2.Multiply('N','N',1.0,C,M_2,0.0);
  // compute CM_2 + M_2C
  Epetra_SerialDenseMatrix CMMC_2(CM_2);  // CMMC_2 = CM_2
  CMMC_2.Multiply('N','N',1.0,M_2,C,1.0); // CMMC_2 = CM_2 + MC_2
  // anisotropic Invariant J4_2 = tr[CM_2]
  double J4_2 = CM_2(0,0) + CM_2(1,1) + CM_2(2,2);

  Epetra_SerialDenseMatrix CsqM_2(3,3);
  CsqM_2.Multiply('N','N',1.0,Csq,M_2,0.0);
  // anisotropic Invariant J5_2 = tr[C^2M_2]
  double J5_2 = CsqM_2(0,0) + CsqM_2(1,1) + CsqM_2(2,2);
  // anisotropic associated reducible Invariant K3_2 = I1*J4_2 - J5_2
  double K3_2 = I1 * J4_2 - J5_2;
  /* end of Second fiber part ***********************************************/

  // cofC = det(C)*C^-T
  Epetra_SerialDenseMatrix CofC(Cinv);
  CofC.Scale(I3);

  /* Underlying strain-energy function
   * Wiso = c1*(\frac{I_1}{I_3^{\frac{1}{3}} -3)                  // ground substance SEF
   *      + eps1 * (I_3^eps2 + \frac{1}{I_3^eps2} - 2)            // penalty function
   *      + alpha1*(I_1*J4 - J5 - 2)^alpha2 if ... >= 2; else 0   // fiber SEF
   */

  // ******* evaluate 2nd PK stress ********************
  Epetra_SerialDenseMatrix S(CofC);   // S = CofC
  double scalar1 = -1.0/3.0 * c1*I1 * pow(I3,-4.0/3.0)
                 + eps1 * eps2 * ( pow(I3,eps2-1.0) - pow(I3,(-eps2-1.0)));
  S.Scale(scalar1);
  double scalar2 = c1*pow(I3,-1.0/3.0);
  if (isnan(scalar1) || isnan(scalar2)){
    cout << "anisotropic_balzani Material evaluated NaN" << endl;
    cout << "I1 = " << I1 << "; I3 = " << I3 << endl;
    cout << "scalar1: " << scalar1 << "; scalar2: " << scalar2 << endl;
    cout << "ElementID: " << ele_ID << "; gp: " << gp << "; time: " << time << endl;
    cout << "GL-Strain: " << *glstrain;
    cout << "Cauchy-Green: " << C;
    dserror("Material computed NaN");
  }
  Epetra_SerialDenseMatrix Iscale(I);
  Iscale.Scale(scalar2);
  S += Iscale;

  // ***** fiber part *******
  if ( (K3 - 2.0) > 1.0E-15){
    // 1st fiber active
    double fiberscalar = alpha1*pow( (K3-2.0),(alpha2-1.0) )*alpha2;
    Epetra_SerialDenseMatrix Mscale(M);
    Mscale.Scale(fiberscalar*I1);       // dW/dJ4 M
    S += Mscale;

    Epetra_SerialDenseMatrix CMMCscale(CMMC);
    CMMCscale.Scale(-fiberscalar);      // dW/dJ5 (CM + MC)
    S += CMMCscale;

    Iscale = I;
    Iscale.Scale(fiberscalar*J4);       // dW/dI1 I just fiber part
    S += Iscale;
  }
  if ( (K3_2 - 2.0) > 1.0E-15){
    // 2nd fiber active
    double fiberscalar = alpha1_2*pow( (K3_2-2.0),(alpha2_2-1.0) )*alpha2_2;
    Epetra_SerialDenseMatrix Mscale(M_2);
    Mscale.Scale(fiberscalar*I1);       // dW/dJ4_2 M_2
    S += Mscale;

    Epetra_SerialDenseMatrix CMMCscale(CMMC_2);
    CMMCscale.Scale(-fiberscalar);      // dW/dJ5_2 (CM_2 + M_2C)
    S += CMMCscale;

    Iscale = I;
    Iscale.Scale(fiberscalar*J4_2);     // dW/dI1 I just fiber part
    S += Iscale;
  }

  // the wellknown factor 2!
  S.Scale(2.0);

  //cout << std::setprecision(10) << S;

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
                    + eps1 * ( pow(I3,eps2-2.0) * eps2*eps2
                             - pow(I3,eps2-2.0) * eps2
                             + eps2*eps2 * pow(I3,-eps2-2.0)
                             + eps2 * pow(I3,-eps2-2.0) );
  double d2W_dI3dI1 = -1.0/3.0 * c1 * pow(I3,-4.0/3.0);
  double dW_dI3     = -1.0/3.0 * c1 * I1 * pow(I3,-4.0/3.0)
                      + eps1 * (pow(I3,eps2-1.0)*eps2 - eps2 * pow(I3,-eps2-1.0));

  ElastSymTensorMultiply(cmat_e,d2W_dI3dI3,CofC,CofC,1.0);     // CofC x CofC
  ElastSymTensorMultiplyAddSym(cmat_e,d2W_dI3dI1,I,CofC,1.0);  // I x CofC + CofC x I
  ElastSymTensorMultiply(cmat_e,I3 * dW_dI3,Cinv,Cinv,1.0);    // Cinv x Cinv
  ElastSymTensor_o_Multiply(cmat_e,-I3 * dW_dI3,Cinv,Cinv,1.0);// - Cinv o Cinv

  // fiber part
  // 1st fiber active
  if ( (K3 - 2.0) > 1.0E-15){
    // compute derivatives of W_ti w.r.t. Invariants
    double K3fac = std::pow( K3-2.0 , alpha2-2.0);
    double d2W_dJ4dJ4 = alpha1 * alpha2 * I1*I1 * (alpha2-1.0) * K3fac;
    double d2W_dJ5dJ5 = alpha1 * alpha2 * (alpha2-1.0) * K3fac;
    double d2W_dI1dJ4 = alpha1 * alpha2 * J4 * I1 * (alpha2-1.0) * K3fac
                      + alpha1 * alpha2 * (K3-2.0) * K3fac;
    double d2W_dI1dJ5 = alpha1 * alpha2 * J4 * (1.0-alpha2) * K3fac;
    double d2W_dJ4dJ5 = alpha1 * alpha2 * I1 * (1.0-alpha2) * K3fac;
    double dW_dJ5     = - alpha1 * alpha2 * (K3-2.0) * K3fac;
    double d2W_dI1dI1 = alpha1 * alpha2 * (alpha2-1.0) * J4*J4 * K3fac;

//    double K3m2sq = 1.0 / (K3 - 2.0) * (K3 - 2.0);
//    double K3m2p = std::pow((K3 - 2.0),alpha2);
//    double d2W_dJ4dJ4 = ( alpha1 * K3m2p * alpha2*alpha2 * I1*I1
//                         -alpha1 * K3m2p * alpha2 * I1*I1 ) * K3m2sq;
//    double d2W_dJ5dJ5 = ( alpha1 * K3m2p * alpha2*alpha2
//                         -alpha1 * K3m2p * alpha2 ) * K3m2sq;
//    double d2W_dI1dJ4 = ( alpha1 * K3m2p * alpha2*alpha2 * J4 * I1
//                         -alpha1 * K3m2p * alpha2 * J4 * I1 ) * K3m2sq
//                         +alpha1 * K3m2p * alpha2 / (K3 - 2.0);
//    double d2W_dI1dJ5 = (-alpha1 * K3m2p * alpha2*alpha2 * J4
//                         +alpha1 * K3m2p * alpha2 * J4 ) * K3m2sq;
//    double d2W_dJ4dJ5 = (-alpha1 * K3m2p * alpha2*alpha2 * I1
//                         +alpha1 * K3m2p * alpha2 * I1 ) * K3m2sq;
//    double dW_dJ5     = - alpha1 * K3m2p * alpha2 / (K3 - 2.0);
//
//    double dW_dI1dI1 = alpha1*alpha2*(alpha2-1) * pow(K3-2.0,alpha2-2.0) * J4*J4;

    // multiply these with corresponding tensor products
    ElastSymTensorMultiply(cmat_e,d2W_dI1dI1,I,I,1.0);         // I x I
    ElastSymTensorMultiply(cmat_e,d2W_dJ4dJ4,M,M,1.0);         // M x M
    ElastSymTensorMultiply(cmat_e,d2W_dJ5dJ5,CMMC,CMMC,1.0);   // (CM+MC) x (CM+MC)
    ElastSymTensorMultiplyAddSym(cmat_e,d2W_dI1dJ4,I,M,1.0);   // I x M + M x I
    ElastSymTensorMultiplyAddSym(cmat_e,d2W_dI1dJ5,CMMC,I,1.0);// (CM+MC) x I + I x (CM+MC)
    ElastSymTensorMultiplyAddSym(cmat_e,d2W_dJ4dJ5,CMMC,M,1.0);// (CM+MC) x M + M x (CM+MC)
    ElastSymTensor_o_Multiply(cmat_e,dW_dJ5,I,M,1.0);          // I o M
    ElastSymTensor_o_Multiply(cmat_e,dW_dJ5,M,I,1.0);          // M o I
  }
  //2nd fiber active
  if ( (K3_2 - 2.0) > 1.0E-15) {
    // temporay factors are NOT renamed to *_2 but refer to 2nd fiber
    // compute derivatives of W_ti w.r.t. Invariants
    double K3fac = std::pow( K3_2-2.0 , alpha2_2-2.0);
    double d2W_dJ4dJ4 = alpha1_2 * alpha2_2 * I1*I1 * (alpha2_2-1.0) * K3fac;
    double d2W_dJ5dJ5 = alpha1_2 * alpha2_2 * (alpha2_2-1.0) * K3fac;
    double d2W_dI1dJ4 = alpha1_2 * alpha2_2 * J4_2 * I1 * (alpha2_2-1.0) * K3fac
                      + alpha1_2 * alpha2_2 * (K3_2-2.0) * K3fac;
    double d2W_dI1dJ5 = alpha1_2 * alpha2_2 * J4_2 * (1.0-alpha2_2) * K3fac;
    double d2W_dJ4dJ5 = alpha1_2 * alpha2_2 * I1 * (1.0-alpha2_2) * K3fac;
    double dW_dJ5     = - alpha1_2 * alpha2_2 * (K3_2-2.0) * K3fac;
    double d2W_dI1dI1 = alpha1_2 * alpha2_2 * (alpha2_2-1.0) * J4_2*J4_2 * K3fac;

    // multiply these with corresponding tensor products
    ElastSymTensorMultiply(cmat_e,d2W_dI1dI1,I,I,1.0);         // I x I
    ElastSymTensorMultiply(cmat_e,d2W_dJ4dJ4,M_2,M_2,1.0);         // M x M
    ElastSymTensorMultiply(cmat_e,d2W_dJ5dJ5,CMMC_2,CMMC_2,1.0);   // (CM+MC) x (CM+MC)
    ElastSymTensorMultiplyAddSym(cmat_e,d2W_dI1dJ4,I,M_2,1.0);   // I x M + M x I
    ElastSymTensorMultiplyAddSym(cmat_e,d2W_dI1dJ5,CMMC_2,I,1.0);// (CM+MC) x I + I x (CM+MC)
    ElastSymTensorMultiplyAddSym(cmat_e,d2W_dJ4dJ5,CMMC_2,M_2,1.0);// (CM+MC) x M + M x (CM+MC)
    ElastSymTensor_o_Multiply(cmat_e,dW_dJ5,I,M_2,1.0);          // I o M
    ElastSymTensor_o_Multiply(cmat_e,dW_dJ5,M_2,I,1.0);          // M o I
  }

  cmat_e.Scale(4.0);

//  cout << (*cmat);
  // end of ********** evaluate C-Matrix *****************************

  return;
}



