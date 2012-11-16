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


#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "Epetra_SerialDenseSolver.h"
#include "mooneyrivlin.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::MooneyRivlin::MooneyRivlin(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c1_(matdata->GetDouble("C1")),
  c2_(matdata->GetDouble("C2")),
  kap_(matdata->GetDouble("KAPPA")),
  lambda_(matdata->GetDouble("LAMBDA")),
  density_(matdata->GetDouble("DENS"))

{
}

Teuchos::RCP<MAT::Material> MAT::PAR::MooneyRivlin::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MooneyRivlin(this));
}


MAT::MooneyRivlinType MAT::MooneyRivlinType::instance_;


DRT::ParObject* MAT::MooneyRivlinType::Create( const std::vector<char> & data )
{
  MAT::MooneyRivlin* moon = new MAT::MooneyRivlin();
  moon->Unpack(data);
  return moon;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)     maf 04/08|
 *----------------------------------------------------------------------*/
MAT::MooneyRivlin::MooneyRivlin()
  : params_(NULL)
{
  dserror("This material law - MOONEY-RIVLIN - is maintained only inside the Elasthyper Toolbox in its regular form (the one in this particular file has an additional summand).\n"
    "If you want to use this law, the material input line should read :\n"
    "MAT 1   MAT_ElastHyper   NUMMAT 1 MATIDS 2 DENS 0\n"
    "MAT 2   ELAST_CoupMooneyRivlin C1 1 C2 1 C3 1 \n");
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)      maf 04/08|
 *----------------------------------------------------------------------*/
MAT::MooneyRivlin::MooneyRivlin(MAT::PAR::MooneyRivlin* params)
  : params_(params)
{
  dserror("This material law - MOONEY-RIVLIN - is maintained only inside the Elasthyper Toolbox in its regular form (the one in this particular file has an additional summand).\n"
    "If you want to use this law, the material input line should read :\n"
    "MAT 1   MAT_ElastHyper   NUMMAT 1 MATIDS 2 DENS 0\n"
    "MAT 2   ELAST_CoupMooneyRivlin C1 1 C2 1 C3 1 \n");
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)     maf 04/08|
 *----------------------------------------------------------------------*/
void MAT::MooneyRivlin::Pack(DRT::PackBuffer& data) const
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
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)     maf 04/08|
 *----------------------------------------------------------------------*/
void MAT::MooneyRivlin::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
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
        params_ = static_cast<MAT::PAR::MooneyRivlin*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
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
  const double c1  = params_->c1_;
  const double c2  = params_->c2_;
  const double kappa_q1 = params_->kap_; // kappa_q1*(J-1)^2
  const double lambda = params_->lambda_;

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


