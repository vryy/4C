/*!----------------------------------------------------------------------
\file yeoh.cpp
\brief

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*----------------------------------------------------------------------*/


#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "Epetra_SerialDenseSolver.h"
#include "yeoh.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::Yeoh::Yeoh(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c1_(matdata->GetDouble("C1")),
  c2_(matdata->GetDouble("C2")),
  c3_(matdata->GetDouble("C3")),
  kap_(matdata->GetDouble("KAPPA")),
  density_(matdata->GetDouble("DENS"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::Yeoh::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Yeoh(this));
}

MAT::YeohType MAT::YeohType::instance_;


DRT::ParObject* MAT::YeohType::Create( const std::vector<char> & data )
{
  MAT::Yeoh* yeoh = new MAT::Yeoh();
  yeoh->Unpack(data);
  return yeoh;
}

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)      tk 01/09|
 *----------------------------------------------------------------------*/
MAT::Yeoh::Yeoh()
  : params_(NULL)
{
  dserror("This material law - YEOH - is maintained only inside the Elasthyper Toolbox.\n"
    "If you want to use this law, the material input line should read :\n"
    "MAT 1   MAT_ElastHyper   NUMMAT 2 MATIDS 2 3 DENS 0\n"
    "MAT 2   ELAST_IsoYeoh C1 1 C2 1 C3 1 \n"
    "MAT 3   ELAST_Vol... the volumetric summand of your choice \n");
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)       tk 01/09|
 *----------------------------------------------------------------------*/
MAT::Yeoh::Yeoh(MAT::PAR::Yeoh* params)
  : params_(params)
{
  dserror("This material law - YEOH - is maintained only inside the Elasthyper Toolbox.\n"
    "If you want to use this law, the material input line should read :\n"
    "MAT 1   MAT_ElastHyper   NUMMAT 2 MATIDS 2 3 DENS 0\n"
    "MAT 2   ELAST_IsoYeoh C1 1 C2 1 C3 1 \n"
    "MAT 3   ELAST_Vol... the volumetric summand of your choice \n");
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)      tk 01/09|
 *----------------------------------------------------------------------*/
void MAT::Yeoh::Pack(DRT::PackBuffer& data) const
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
 |  Unpack                                        (public)      tk 01/09|
 *----------------------------------------------------------------------*/
void MAT::Yeoh::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::Yeoh*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}


/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)      tk 01/09|
 *----------------------------------------------------------------------*/
void MAT::Yeoh::Evaluate(
        const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
        LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
        LINALG::Matrix<NUM_STRESS_3D,1> * stress)
{
  // get material parameters
  const double c1  = params_->c1_;
  const double c2  = params_->c2_;
  const double c3  = params_->c3_;
  const double kappa = params_->kap_;

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
        + 2.0*C(3)*C(4)*C(5)
        - C(1)*C(5)*C(5)
        - C(2)*C(3)*C(3)
        - C(0)*C(4)*C(4);    // 3rd invariant, determinant

  const double J = sqrt(I3);
  const double incJ = std::pow(I3,-1.0/3.0);  // J^{-2/3}

  const double I1bar= incJ*I1; //first invariant of modified right Cauchy-Green


  // invert C
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv(false);

  Cinv(0) = C(1)*C(2) - C(4)*C(4);
  Cinv(1) = C(0)*C(2) - C(5)*C(5);
  Cinv(2) = C(0)*C(1) - C(3)*C(3);
  Cinv(3) = C(5)*C(4) - C(3)*C(2);
  Cinv(4) = C(3)*C(5) - C(0)*C(4);
  Cinv(5) = C(3)*C(4) - C(5)*C(1);
  Cinv.Scale(1.0/I3);

  // *** strain energy function  ************************************************
  // Yeoh (Holzapfel, p.248)
  // W = c1 (I1bar-3) + c2 (I1bar-3)^2 + c3 (I1bar-3)^3 +1/2 kappa (J-1)^2

  // *** 2nd PK stresses ********************************************************
  // S = Svol + Siso
  // Svol = J*kappa*(J-1)
  // Isochoric (deviatoric) part via projection PP:Sbar, see Holzapfel p. 230
  // Siso = J^{-2/3}  Dev[Sbar] = J^{-2/3} [Sbar - 1/3 trace(Sbar C) Cinv]
  // here: Sbar = (2* c1 + 4 * c2* (I1bar-3) + 6* c3 * (I1bar-3)^2) * Id  (Holzapfel, p. 249)

  const double third = 1./3.;
  const double p = kappa*(J-1);
  const double gamma1 = 2.*c1 + 4.*c2*(I1bar-3.) + 6.*c3*(I1bar-3.)*(I1bar-3.);
  for (int i = 0; i < 6; ++i) {
    (*stress)(i) = J*p * Cinv(i);  // volumetric part
    (*stress)(i) += incJ* (gamma1*Id(i) - third*gamma1*I1*Cinv(i));  //isochoric part
  }

  // *** Elasticity =  CCvol + CCiso *******************************************
  // CCvol = J(p + J dp/dJ) Cinv x Cinv  -  2 J p Cinv o Cinv
  // CCiso = PP:CCbar::PP^T + 2/3 J^{-2/3} Sbar:C Psl - 2/3 (Cinv x Siso + Siso x Cinv)
  // with CCbar = (2*J^{-4/3}* dSbar/dCbar) (see Holzapfel p. 255)

  //first part of volumetric CC
  AddtoCmatHolzapfelProduct((*cmat),Cinv,(-2*J*p));  // -2 J p Cinv o Cinv

  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>  Psl(true);        // Psl = Cinv o Cinv - 1/3 Cinv x Cinv
  AddtoCmatHolzapfelProduct(Psl,Cinv,1.0);  // first part Psl = Cinv o Cinv

  const double delta1 = 8.*c2+24.*c3*(I1bar-3.); //dSbar/dCbar
  const double alpha = incJ*incJ*delta1;
  const double fac = 2.*third*incJ*gamma1*I1;  // 2/3 J^{-2/3} Sbar:C

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      //second part of volumetric CC
      (*cmat)(i,j) += J*(p+J*kappa) * Cinv(i) * Cinv(j);  // J(p + J dp/dJ) Cinv x Cinv
      // on the fly complete Psl needed later
      Psl(i,j) += (-third) * Cinv(i) * Cinv(j);

      //isochoric part: PP:CCbar::PP^T
      (*cmat)(i,j) += alpha * ( Id(i)*Id(j) //alpha* Id x Id
                      - third * I1 * Id(i)*Cinv(j) // alpha*I1 / 3 * ID x Cinv
                      - third * I1 * Id(j)*Cinv(i) // alpha*I1 / 3 * Cinv x Id
                      + third * third * I1 * I1 * Cinv(i) * Cinv(j)); // alpha*I1^2/9*Cinv x Cinv

      //isochoric parts: 2/3 J^{-2/3} Sbar:C Psl - 2/3 (Cinv x Siso + Siso x Cinv)
      (*cmat)(i,j) += fac * Psl(i,j)                            // fac Psl
                      - 2*third * Cinv(i) * incJ * (gamma1*Id(j) - third*gamma1*I1*Cinv(j)) // -2/3 Cinv x Siso
                      - 2*third * Cinv(j) * incJ * (gamma1*Id(i) - third*gamma1*I1*Cinv(i));// -2/3 Siso x Cinv
    }
  }




  // do the dirty scaling!!!

//  (*stress).Scale(c2);
//  (*cmat).Scale(c2);

  return;
}


