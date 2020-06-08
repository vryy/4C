/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required for OGDEN-like material model
according to GASSER: "Failure properties of intraluminal thrombus in abdominal
aortic aneurysm under static and pulsating mechanical load", Journal of Vascular
Surgery Volume 48, Number 1, July 2008.

the input line should read:
  MAT 1 MAT_Struct_AAAGasser DENS 0.0001 VOL OgSiMi NUE 0.49 BETA -2.0 CLUM 2.62e-3 CMED 2.13e-3
CABLUM 1.98e-3

\level 3

\maintainer Christoph Schmidt


*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include "aaagasser.H"
#include "matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::AAAgasser::AAAgasser(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),

      density_(matdata->GetDouble("DENS")),
      vol_(matdata->Get<std::string>("VOL")),
      nue_(matdata->GetDouble("NUE")),
      beta_(matdata->GetDouble("BETA")),
      Clum_(matdata->GetDouble("CLUM")),
      Cmed_(matdata->GetDouble("CMED")),
      Cablum_(matdata->GetDouble("CABLUM"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::AAAgasser::CreateMaterial()
{
  return Teuchos::rcp(new MAT::AAAgasser(this));
}

MAT::AAAgasserType MAT::AAAgasserType::instance_;


DRT::ParObject* MAT::AAAgasserType::Create(const std::vector<char>& data)
{
  MAT::AAAgasser* aaa = new MAT::AAAgasser();
  aaa->Unpack(data);
  return aaa;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)              |
 *----------------------------------------------------------------------*/
MAT::AAAgasser::AAAgasser() : params_(NULL) {}


/*----------------------------------------------------------------------*
 |  Constructor                                    (public)             |
 *----------------------------------------------------------------------*/
MAT::AAAgasser::AAAgasser(MAT::PAR::AAAgasser* params) : params_(params) {}

/*----------------------------------------------------------------------*
 |  Pack                                           (public)             |
 *----------------------------------------------------------------------*/
void MAT::AAAgasser::Pack(DRT::PackBuffer& data) const
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
 |  Unpack                                         (public)             |
 *----------------------------------------------------------------------*/
void MAT::AAAgasser::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::AAAgasser*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                         (public)                  |
 *----------------------------------------------------------------------*/
void MAT::AAAgasser::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int gp, const int eleGID)
{
  double normdist = params.get("iltthick meanvalue", -999.0);
  if (normdist == -999.0) dserror("Aneurysm mean ilt distance not found");

  // material parameters for isochoric part:
  // calculate element stiffness parameter, in dependence of 'normdist':
  double Cele = 0.0;
  if (0.0 <= normdist and normdist <= 0.5)
    Cele = ((normdist - 0.5) / (-0.5)) * params_->Clum_ + (normdist / 0.5) * params_->Cmed_;
  else if (0.5 < normdist and normdist <= 1.0)
    Cele =
        ((normdist - 1.0) / (-0.5)) * params_->Cmed_ + ((normdist - 0.5) / 0.5) * params_->Cablum_;
  else
    dserror("Unable to calculate valid stiffness parameter in material AAAGasser");

  // double Cele        = 1./3.*(params_->Clum_ + params_->Cmed_ + params_->) * normdist/normdist;

  // material parameters for volumetric part:
  const std::string vol = *params_->vol_;
  // kappa = youngs / (3-6*nue); youngs = 24*Cele (see GASSER p.184);
  const double kappa = 24 * Cele / (3 - 6 * params_->nue_);  // dilatational modulus
  const double beta = params_->beta_;                        // parameter from Holzapfel

  //--------------------------------------------------------------------------------------
  // build second order identity tensor:
  LINALG::Matrix<6, 1> id2(true);
  for (int i = 0; i < 3; i++) id2(i) = 1.0;

  // build fourth order identity tensor S (see HOLZAPFEL p. 261):
  LINALG::Matrix<6, 6> id4(true);
  for (int i = 0; i < 3; i++) id4(i, i) = 1.0;
  for (int i = 3; i < 6; i++) id4(i, i) = 0.5;

  // GREEN-LAGRANGE strain tensor:
  //  [    glstrain(0)    0.5*glstrain(3)   0.5*glstrain(5)]
  //  [0.5*glstrain(3)        glstrain(1)   0.5*glstrain(4)]
  //  [0.5*glstrain(5)    0.5*glstrain(4)       glstrain(2)]

  // right Cauchy-Green Tensor:
  LINALG::Matrix<6, 1> rcg(*glstrain);
  rcg.Scale(2.0);
  rcg.Update(1.0, id2, 1.0);

  // 'contra-variant' right Cauchy-Green Tensor (originally from bborn 08/09)
  LINALG::Matrix<6, 1> scg(rcg);
  for (int i = 3; i < 6; i++) scg(i) *= 0.5;

  //--------------------------------------------------------------------------------------
  // principal invariants
  // 1st invariant, trace:
  double inv = rcg(0) + rcg(1) + rcg(2);
  // 2nd invariant:
  double iinv = 0.5 * (inv * inv - rcg(0) * rcg(0) - rcg(1) * rcg(1) - rcg(2) * rcg(2) -
                          0.5 * rcg(3) * rcg(3) - 0.5 * rcg(4) * rcg(4) - 0.5 * rcg(5) * rcg(5));
  // 3rd invariant, determinante (SARRUS):
  double iiinv = rcg(0) * rcg(1) * rcg(2) + 0.25 * rcg(3) * rcg(4) * rcg(5) -
                 0.25 * rcg(1) * rcg(5) * rcg(5) - 0.25 * rcg(2) * rcg(3) * rcg(3) -
                 0.25 * rcg(0) * rcg(4) * rcg(4);

  double detf = 0.0;
  if (iiinv < 0.0)
    dserror("fatal failure in GASSER thrombus material (detf<0)");
  else
    detf = sqrt(iiinv);  // determinate of deformation gradient -> J

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


  //--- determine 2nd Piola Kirchhoff stresses pktwo -------------------------------------
  stress->Clear();

  // 1st step: isochoric part (HOLZAPFEL S.248)
  //=========================
  double gamma1 = 0.0;
  double gamma2 = 4.0 * Cele * pow(iiinv, -2. / 3.);
  double gamma3 = 4. / 3. * Cele * pow(iiinv, -2. / 3.) * (-inv * inv + 2 * iinv);

  // contribution: I
  stress->Update(gamma1, id2, 1.0);
  // contribution: C
  stress->Update(gamma2, scg, 1.0);
  // contribution: Cinv
  stress->Update(gamma3, invc, 1.0);


  // 2nd step: volumetric part (HOLZAPFEL S.230)
  //==========================
  double pres = -999.0;
  double prestild = -999.0;

  // choose volumetric strain energy form
  if (vol.compare("OSM") == 0)
  {
    pres = kappa * (1 - pow(detf, -beta)) / (detf * beta);
    prestild = kappa * pow(detf, -beta - 1);
  }
  else if (vol.compare("SuBa") == 0)
  {
    pres = kappa * (detf - 1);
    prestild = kappa * (2 * detf - 1);
  }
  else if (vol.compare("SiTa") == 0)
  {
    pres = kappa * (detf + (1 / detf) * log(detf) - 1) / 2;
    prestild = kappa * (2 * detf + (1 / detf) - 1) / 2;
  }
  else
    dserror("Choose OSM, SuBa or SiTa for the volumetric part! See reference...!");

  stress->Update(pres * detf, invc, 1.0);



  //--- do elasticity matrix -------------------------------------------------------------
  cmat->Clear();

  // 1st step: isochoric part (HOLZAPFEL S.261)
  //=========================
  double delta5 = -16. / 3. * Cele * pow(iiinv, -2. / 3.);
  double delta6 = 16. / 9. * Cele * pow(iiinv, -2. / 3.) * (inv * inv - 2 * iinv);
  double delta7 = 8. / 3. * Cele * pow(iiinv, -2. / 3.) * (inv * inv - 2 * iinv);
  double delta8 = 8.0 * Cele * pow(iiinv, -2. / 3.);

  // contribution: (C \otimes Cinv + Cinv \otimes C)
  cmat->MultiplyNT(delta5, scg, invc, 1.0);
  cmat->MultiplyNT(delta5, invc, scg, 1.0);
  // contribution: Cinv \otimes Cinv
  cmat->MultiplyNT(delta6, invc, invc, 1.0);
  // contribution: Cinv \odot Cinv
  AddtoCmatHolzapfelProduct(*cmat, invc, delta7);
  // contribution: S
  cmat->Update(delta8, id4, 1.0);


  // 2nd step: volumetric part (HOLZAPFEL S.254f)
  //==========================
  // contribution: J prestild Cinv \otimes Cinv
  cmat->MultiplyNT(detf * prestild, invc, invc, 1.0);
  // contribution: -2 J*p Cinv \odot Cinv
  AddtoCmatHolzapfelProduct(*cmat, invc, -2 * detf * pres);


  return;
}
