/*!----------------------------------------------------------------------
\brief
This file contains the routines required for aneurysmatic artery wall following
Raghavan and Vorp [2000]

The input line should read
 MAT 1 MAT_Struct_AAA_MixedEffects AGE 67 REFDIA 22.5 NUE 0.49 DENS 0.0001

\level 3

\maintainer Christoph Schmidt

*----------------------------------------------------------------------*/

#include "aaa_mixedeffects.H"
#include "matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::AAA_mixedeffects::AAA_mixedeffects(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      nue_(matdata->GetDouble("NUE")),
      age_(matdata->GetDouble("AGE")),
      refdia_(matdata->GetDouble("REFDIA")),
      density_(matdata->GetDouble("DENS"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::AAA_mixedeffects::CreateMaterial()
{
  return Teuchos::rcp(new MAT::AAA_mixedeffects(this));
}

MAT::AAA_mixedeffectsType MAT::AAA_mixedeffectsType::instance_;


DRT::ParObject* MAT::AAA_mixedeffectsType::Create(const std::vector<char>& data)
{
  MAT::AAA_mixedeffects* aaa = new MAT::AAA_mixedeffects();
  aaa->Unpack(data);
  return aaa;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  chfoe 03/08 |
 *----------------------------------------------------------------------*/
MAT::AAA_mixedeffects::AAA_mixedeffects() : params_(NULL) {}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   chfoe 03/08 |
 *----------------------------------------------------------------------*/
MAT::AAA_mixedeffects::AAA_mixedeffects(MAT::PAR::AAA_mixedeffects* params) : params_(params) {}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)  chfoe 03/08 |
 *----------------------------------------------------------------------*/
void MAT::AAA_mixedeffects::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
 |  Unpack                                        (public)  chfoe 03/08 |
 *----------------------------------------------------------------------*/
void MAT::AAA_mixedeffects::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;

  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::AAA_mixedeffects*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }


  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*
 |  Evaluate Material                   (public)  chfoe 03/08 gee 10/08 |
 *----------------------------------------------------------------------*

 plain strain energy function

 W    = alpha (Ic*IIIc^(-1/3) -3) + beta (Ic*IIIc^(-1/3)-3)²

 taken from
 M.L. Raghavan, D.A. Vorp: Toward a biomechanical tool to evaluate rupture potential
 of abdominal aortic aneurysm: identification of a finite strain constitutive model
 and evaluation of its applicability, J. of Biomechanics 33 (2000) 475-482.

 and modified to slight compressibility

 here

 Ic   .. first invariant of right Cauchy-Green tensor C
 IIIc .. third invariant of right Cauchy-Green tensor C

 The volumetric part is done by a volumetric strain energy function taken from
 Holzapfel

 W_vol = K beta2^(-2) ( beta2 ln (J) + J^(-beta2) -1 )

 where

 K    .. bulk modulus
 beta2 = 9.0 a parameter according to Holzapfel
 J    .. det(F) determinante of the Jacobian matrix


 Note: Young's modulus is in the input just for convenience. Actually we need the
       parameter alpha (see W above) which is related to E by

     E = 6.0 * alpha.

       Correspondingly the bulk modulus is given by

     K = E / (3-6*nu) = 2*alpha / (1-2*nu)

     with nu = 0.495 we have K = 200 alpha
     with nu = 0.45  we have K =  20 alpha

 */
void MAT::AAA_mixedeffects::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int eleGID)
{
  double elelocalrad = params.get("localrad meanvalue", -999.0);
  if (elelocalrad == -999.0) dserror("Aneurysm local radii not found");

  // material parameters for isochoric part
  const double alpha = 1E6 * (0.09631 + 0.03329 * (elelocalrad * 2 / params_->refdia_ - 2.55));
  const double beta =
      1E6 * (-0.9553 * (elelocalrad * 2 / (params_->refdia_) - 2.55) + (0.06721 * (params_->age_)));
  // std::cout << elelocalrad*2 << "  " << beta << std::endl;
  const double nue = params_->nue_;  // Poisson's ratio
  // material parameters for volumetric part
  const double beta2 = -2.0;  // parameter from Holzapfel
  const double komp = (nue != 0.5) ? 2.0 * alpha / (1.0 - 2.0 * nue) : 0.0;  // bulk modulus

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<6, 1> identity(true);
  for (int i = 0; i < 3; i++) identity(i) = 1.0;

  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<6, 1> rcg(*glstrain);
  rcg.Scale(2.0);
  rcg += identity;

  // invariants
  double inv = rcg(0) + rcg(1) + rcg(2);  // 1st invariant, trace
  double iiinv = rcg(0) * rcg(1) * rcg(2) + 0.25 * rcg(3) * rcg(4) * rcg(5) -
                 0.25 * rcg(1) * rcg(5) * rcg(5) - 0.25 * rcg(2) * rcg(3) * rcg(3) -
                 0.25 * rcg(0) * rcg(4) * rcg(4);  // 3rd invariant, determinante

  double detf = 0.0;
  if (iiinv < 0.0)
    dserror("fatal failure in aneurysmatic artery wall material");
  else
    detf = sqrt(iiinv);  // determinate of deformation gradient

  //--------------------------------------------------------------------------------------
  // invert C
  LINALG::Matrix<6, 1> invc(false);

  double invdet = 1. / iiinv;

  invc(0) = rcg(1) * rcg(2) - 0.25 * rcg(4) * rcg(4);
  invc(1) = rcg(0) * rcg(2) - 0.25 * rcg(5) * rcg(5);
  invc(2) = rcg(0) * rcg(1) - 0.25 * rcg(3) * rcg(3);
  invc(3) = 0.25 * rcg(5) * rcg(4) - 0.5 * rcg(3) * rcg(2);
  invc(4) = 0.25 * rcg(3) * rcg(5) - 0.5 * rcg(0) * rcg(4);
  invc(5) = 0.25 * rcg(3) * rcg(4) - 0.5 * rcg(5) * rcg(1);

  invc.Scale(invdet);

  //--- prepare some constants -----------------------------------------------------------
  const double third = 1.0 / 3.0;
  const double twthi = 2.0 / 3.0;


  //--- determine 2nd Piola Kirchhoff stresses pktwo -------------------------------------
  // 1st step: isochoric part
  //=========================
  double isochor1 =
      2.0 * (alpha * pow(iiinv, third) + 2.0 * beta * inv - 6.0 * beta * pow(iiinv, third)) *
      pow(iiinv, -twthi);
  double isochor2 =
      -twthi * inv *
      (alpha * pow(iiinv, third) + 2.0 * beta * inv - 6.0 * beta * pow(iiinv, third)) *
      pow(iiinv, -twthi);

  // contribution: Cinv
  LINALG::Matrix<6, 1> pktwoiso(invc);
  pktwoiso.Scale(isochor2);

  // contribution: I
  for (int i = 0; i < 3; i++) pktwoiso(i) += isochor1;


  // 2nd step: volumetric part
  //==========================
  double scalar = komp / beta2 * (1.0 - pow(detf, -beta2));

  // initialise PKtwo with volumetric part
  LINALG::Matrix<6, 1> pktwovol(invc);
  pktwovol.Scale(scalar);

  // 3rd step: add everything up
  //============================
  (*stress) = pktwoiso;
  (*stress) += pktwovol;


  //--- do elasticity matrix -------------------------------------------------------------
  // ensure that cmat is zero when it enters the computation
  // It is an implicit law that cmat is zero upon input
  // cmat.PutScalar(0.0);

  // 1st step: isochoric part
  //=========================

  // deltas (see also Holzapfel p.261)
  // note that these deltas serve for the isochoric part only
  double delta1 = 8.0 * beta * pow(iiinv, -twthi);
  double delta3 = -4. / 3 *
                  (alpha * pow(iiinv, third) + 4. * beta * inv - 6 * beta * pow(iiinv, third)) *
                  pow(iiinv, -twthi);
  double delta6 = 4. / 9 * inv *
                  (alpha * pow(iiinv, third) + 4. * beta * inv - 6 * beta * pow(iiinv, third)) *
                  pow(iiinv, -twthi);
  double delta7 = 4. / 3 * inv *
                  (alpha * pow(iiinv, third) + 2. * beta * inv - 6 * beta * pow(iiinv, third)) *
                  pow(iiinv, -twthi);

  // contribution: I \obtimes I
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) (*cmat)(i, j) = delta1;

  // contribution: Cinv \otimes Cinv
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
    {
      // contribution: Cinv \otimes I + I \otimes Cinv
      (*cmat)(i, j) += delta3 * (identity(i) * invc(j) + invc(i) * identity(j));
      // contribution: Cinv \otimes Cinv
      (*cmat)(i, j) += delta6 * invc(i) * invc(j);
    }

  // contribution: boeppel-product
  AddtoCmatHolzapfelProduct(*cmat, invc, delta7);

  // 2nd step: volumetric part
  //==========================
  delta6 = komp * pow(detf, -beta2);
  delta7 = -2.0 * scalar;

  // contribution: Cinv \otimes Cinv
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++) (*cmat)(i, j) += delta6 * invc(i) * invc(j);

  // contribution: boeppel-product
  AddtoCmatHolzapfelProduct(*cmat, invc, delta7);

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)  chfoe 03/08 |
 *----------------------------------------------------------------------*

 plain strain energy function

 W    = alpha (Ic*IIIc^(-1/3) -3) + beta (Ic*IIIc^(-1/3)-3)²

 taken from
 M.L. Raghavan, D.A. Vorp: Toward a biomechanical tool to evaluate rupture potential
 of abdominal aortic aneurysm: identification of a finite strain constitutive model
 and evaluation of its applicability, J. of Biomechanics 33 (2000) 475-482.

 and modified to slight compressibility

 here

 Ic   .. first invariant of right Cauchy-Green tensor C
 IIIc .. third invariant of right Cauchy-Green tensor C

 The volumetric part is done by a volumetric strain engergy function taken from
 Holzapfel

 W_vol = K beta2^(-2) ( beta2 ln (J) + J^(-beta2) -1 )

 where

 K    .. bulk modulus
 beta2 = 9.0 a parameter according to Holzapfel
 J    .. det(F) determinante of the Jacobian matrix


 Note: Young's modulus is in the input just for convenience. Actually we need the
       parameter alpha (see W above) which is related to E by

     E = 6.0 * alpha.

       Correspondingly the bulk modulus is given by

     K = E / (3-6*nu) = 2*alpha / (1-2*nu)

     with nu = 0.495 we have K = 200 alpha
     with nu = 0.45  we have K =  20 alpha

 */
void MAT::AAA_mixedeffects::Evaluate(const Epetra_SerialDenseVector* glstrain_e,
    Epetra_SerialDenseMatrix* cmat_e, Epetra_SerialDenseVector* stress_e, double elelocalrad)
{
  // this is temporary as long as the material does not have a
  // Matrix-type interface
  const LINALG::Matrix<6, 1> glstrain(glstrain_e->A(), true);
  LINALG::Matrix<6, 6> cmat(cmat_e->A(), true);
  LINALG::Matrix<6, 1> stress(stress_e->A(), true);

  // material parameters for isochoric part
  const double alpha = 1E6 * (0.09631 + 0.03329 * (elelocalrad * 2 / params_->refdia_ - 2.55));
  const double beta =
      1E6 * (-0.9553 * (elelocalrad / (params_->refdia_) - 2.55) + (0.06721 * (params_->age_)));
  const double nue = params_->nue_;  // Poisson's ratio
  // material parameters for volumetric part
  double beta2 = 9.0;                             // parameter from Holzapfel
  double komp = 2.0 * alpha / (1.0 - 2.0 * nue);  // bulk modulus

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<6, 1> identity(true);
  for (int i = 0; i < 3; i++) identity(i) = 1.0;

  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<6, 1> rcg(glstrain);
  rcg.Scale(2.0);
  rcg += identity;

  // invariants
  double inv = rcg(0) + rcg(1) + rcg(2);  // 1st invariant, trace
  double iiinv = rcg(0) * rcg(1) * rcg(2) + 0.25 * rcg(3) * rcg(4) * rcg(5) -
                 0.25 * rcg(1) * rcg(5) * rcg(5) - 0.25 * rcg(2) * rcg(3) * rcg(3) -
                 0.25 * rcg(0) * rcg(4) * rcg(4);  // 3rd invariant, determinante

  double detf = 0.0;
  if (iiinv < 0.0)
    dserror("fatal failure in aneurysmatic artery wall material");
  else
    detf = sqrt(iiinv);  // determinate of deformation gradient

  //--------------------------------------------------------------------------------------
  // invert C
  LINALG::Matrix<6, 1> invc(false);

  double invdet = 1. / iiinv;

  invc(0) = rcg(1) * rcg(2) - 0.25 * rcg(4) * rcg(4);
  invc(1) = rcg(0) * rcg(2) - 0.25 * rcg(5) * rcg(5);
  invc(2) = rcg(0) * rcg(1) - 0.25 * rcg(3) * rcg(3);
  invc(3) = 0.25 * rcg(5) * rcg(4) - 0.5 * rcg(3) * rcg(2);
  invc(4) = 0.25 * rcg(3) * rcg(5) - 0.5 * rcg(0) * rcg(4);
  invc(5) = 0.25 * rcg(3) * rcg(4) - 0.5 * rcg(5) * rcg(1);

  invc.Scale(invdet);

  //--- prepare some constants -----------------------------------------------------------
  const double third = 1.0 / 3.0;
  const double twthi = 2.0 / 3.0;


  //--- determine 2nd Piola Kirchhoff stresses pktwo -------------------------------------
  // 1st step: isochoric part
  //=========================
  double isochor1 =
      2.0 * (alpha * pow(iiinv, third) + 2.0 * beta * inv - 6.0 * beta * pow(iiinv, third)) *
      pow(iiinv, -twthi);
  double isochor2 =
      -twthi * inv *
      (alpha * pow(iiinv, third) + 2.0 * beta * inv - 6.0 * beta * pow(iiinv, third)) *
      pow(iiinv, -twthi);

  // contribution: Cinv
  LINALG::Matrix<6, 1> pktwoiso(invc);
  pktwoiso.Scale(isochor2);

  // contribution: I
  for (int i = 0; i < 3; i++) pktwoiso(i) += isochor1;


  // 2nd step: volumetric part
  //==========================
  double scalar = komp / beta2 * (1.0 - pow(detf, -beta2));

  // initialise PKtwo with volumetric part
  LINALG::Matrix<6, 1> pktwovol(invc);
  pktwovol.Scale(scalar);

  // 3rd step: add everything up
  //============================
  stress = pktwoiso;
  stress += pktwovol;


  //--- do elasticity matrix -------------------------------------------------------------
  // ensure that cmat is zero when it enters the computation
  // It is an implicit law that cmat is zero upon input
  // cmat.PutScalar(0.0);

  // 1st step: isochoric part
  //=========================

  // deltas (see also Holzapfel p.261)
  // note that these deltas serve for the isochoric part only
  double delta1 = 8.0 * beta * pow(iiinv, -twthi);
  double delta3 = -4. / 3 *
                  (alpha * pow(iiinv, third) + 4. * beta * inv - 6 * beta * pow(iiinv, third)) *
                  pow(iiinv, -twthi);
  double delta6 = 4. / 9 * inv *
                  (alpha * pow(iiinv, third) + 4. * beta * inv - 6 * beta * pow(iiinv, third)) *
                  pow(iiinv, -twthi);
  double delta7 = 4. / 3 * inv *
                  (alpha * pow(iiinv, third) + 2. * beta * inv - 6 * beta * pow(iiinv, third)) *
                  pow(iiinv, -twthi);

  // contribution: I \obtimes I
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) cmat(i, j) = delta1;

  // contribution: Cinv \otimes Cinv
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
    {
      // contribution: Cinv \otimes I + I \otimes Cinv
      cmat(i, j) += delta3 * (identity(i) * invc(j) + invc(i) * identity(j));
      // contribution: Cinv \otimes Cinv
      cmat(i, j) += delta6 * invc(i) * invc(j);
    }

  // contribution: boeppel-product
  AddtoCmatHolzapfelProduct(cmat, invc, delta7);

  // 2nd step: volumetric part
  //==========================
  delta6 = komp * pow(detf, -beta2);
  delta7 = -2.0 * scalar;

  // contribution: Cinv \otimes Cinv
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++) cmat(i, j) += delta6 * invc(i) * invc(j);

  // contribution: boeppel-product
  AddtoCmatHolzapfelProduct(cmat, invc, delta7);

  return;
}
