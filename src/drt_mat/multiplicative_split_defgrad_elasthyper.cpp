/*----------------------------------------------------------------------*/
/*! \file
\brief evaluation of a generic material whose deformation gradient is modeled to be split
multiplicatively into elastic and inelastic parts

\level 3

*/
/*----------------------------------------------------------------------*/

#include "multiplicative_split_defgrad_elasthyper.H"

#include "anisotropy.H"
#include "elasthyper_service.H"
#include "inelastic_defgrad_factors.H"
#include "material_service.H"
#include "matpar_bundle.H"
#include "multiplicative_split_defgrad_elasthyper_service.H"

#include "../drt_inpar/inpar_ssi.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_structure_new/str_enum_lists.H"

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::MultiplicativeSplitDefgrad_ElastHyper::MultiplicativeSplitDefgrad_ElastHyper(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      nummat_elast_(matdata->GetInt("NUMMATEL")),
      matids_elast_(matdata->Get<std::vector<int>>("MATIDSEL")),
      numfac_inel_(matdata->GetInt("NUMFACINEL")),
      inel_defgradfacids_(matdata->Get<std::vector<int>>("INELDEFGRADFACIDS")),
      density_(matdata->GetDouble("DENS"))
{
  // check if sizes fit
  if (nummat_elast_ != static_cast<int>(matids_elast_->size()))
    dserror("number of elastic materials %d does not fit to size of elastic material ID vector %d",
        nummat_elast_, matids_elast_->size());

  if (numfac_inel_ != static_cast<int>(inel_defgradfacids_->size()))
  {
    dserror(
        "number of inelastic deformation gradient factors %d does not fit to size of inelastic "
        "deformation gradient ID vector %d",
        numfac_inel_, inel_defgradfacids_->size());
  }
}

Teuchos::RCP<MAT::Material> MAT::PAR::MultiplicativeSplitDefgrad_ElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MultiplicativeSplitDefgrad_ElastHyper(this));
}

MAT::MultiplicativeSplitDefgrad_ElastHyperType
    MAT::MultiplicativeSplitDefgrad_ElastHyperType::instance_;

DRT::ParObject* MAT::MultiplicativeSplitDefgrad_ElastHyperType::Create(
    const std::vector<char>& data)
{
  auto* splitdefgrad_elhy = new MAT::MultiplicativeSplitDefgrad_ElastHyper();
  splitdefgrad_elhy->Unpack(data);

  return splitdefgrad_elhy;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::MultiplicativeSplitDefgrad_ElastHyper::MultiplicativeSplitDefgrad_ElastHyper()
    : anisotropy_(Teuchos::rcp(new MAT::Anisotropy())),
      inelastic_(Teuchos::rcp(new MAT::InelasticFactorsHandler())),
      params_(nullptr),
      potsumel_(0)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::MultiplicativeSplitDefgrad_ElastHyper::MultiplicativeSplitDefgrad_ElastHyper(
    MAT::PAR::MultiplicativeSplitDefgrad_ElastHyper* params)
    : anisotropy_(Teuchos::rcp(new MAT::Anisotropy())),
      inelastic_(Teuchos::rcp(new MAT::InelasticFactorsHandler())),
      params_(params),
      potsumel_(0)
{
  // elastic materials
  for (int matid_elastic : *(params_->matids_elast_))
  {
    auto elastic_summand = MAT::ELASTIC::Summand::Factory(matid_elastic);
    if (elastic_summand == Teuchos::null) dserror("Failed to allocate");
    potsumel_.push_back(elastic_summand);
    elastic_summand->RegisterAnisotropyExtensions(*anisotropy_);
  }

  inelastic_->Setup(params);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::MultiplicativeSplitDefgrad_ElastHyper::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  anisotropy_->PackAnisotropy(data);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (const auto& p : potsumel_) p->PackSummand(data);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::MultiplicativeSplitDefgrad_ElastHyper::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  params_ = nullptr;
  potsumel_.clear();

  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      auto* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<MAT::PAR::MultiplicativeSplitDefgrad_ElastHyper*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  anisotropy_->UnpackAnisotropy(data, position);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // elastic materials
    for (const auto& matid_elastic : *(params_->matids_elast_))
    {
      auto elastic_summand = MAT::ELASTIC::Summand::Factory(matid_elastic);
      if (elastic_summand == Teuchos::null) dserror("Failed to allocate");
      potsumel_.push_back(elastic_summand);
    }
    // loop map of associated potential summands
    for (const auto& elastic_summand : potsumel_)
    {
      elastic_summand->UnpackSummand(data, position);
      elastic_summand->RegisterAnisotropyExtensions(*anisotropy_);
    }

    // inelastic deformation gradient factors
    inelastic_->Setup(params_);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::MultiplicativeSplitDefgrad_ElastHyper::Evaluate(const LINALG::Matrix<3, 3>* const defgrad,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int gp, const int eleGID)
{
  // do all stuff that only has to be done once per Evaluate() call
  PreEvaluate(params, gp);

  // static variables
  static LINALG::Matrix<6, 6> cmatiso(true);
  static LINALG::Matrix<6, 9> dSdiFin(true);
  static LINALG::Matrix<6, 6> cmatadd(true);

  // build inverse inelastic deformation gradient
  static LINALG::Matrix<3, 3> iFinM(true);
  inelastic_->EvaluateInverseInelasticDefGrad(defgrad, iFinM);

  // determinante of inelastic deformation gradient
  const double detFin = 1.0 / iFinM.Determinant();

  // static variables of kinetic quantities
  static LINALG::Matrix<6, 1> iCV(true);
  static LINALG::Matrix<6, 1> iCinV(true);
  static LINALG::Matrix<6, 1> iCinCiCinV(true);
  static LINALG::Matrix<3, 3> iCinCM(true);
  static LINALG::Matrix<3, 3> iFinCeM(true);
  static LINALG::Matrix<9, 1> CiFin9x1(true);
  static LINALG::Matrix<9, 1> CiFinCe9x1(true);
  static LINALG::Matrix<9, 1> CiFiniCe9x1(true);
  static LINALG::Matrix<3, 1> prinv(true);
  EvaluateKinQuantElast(defgrad, iFinM, iCinV, iCinCiCinV, iCV, iCinCM, iFinCeM, CiFin9x1,
      CiFinCe9x1, CiFiniCe9x1, prinv);

  // derivatives of principle invariants
  static LINALG::Matrix<3, 1> dPIe(true);
  static LINALG::Matrix<6, 1> ddPIIe(true);
  EvaluateInvariantDerivatives(prinv, gp, eleGID, dPIe, ddPIIe);

  // 2nd Piola Kirchhoff stresses factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
  static LINALG::Matrix<3, 1> gamma(true);
  // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
  static LINALG::Matrix<8, 1> delta(true);
  // compose coefficients
  CalculateGammaDelta(gamma, delta, prinv, dPIe, ddPIIe);

  // evaluate dSdiFin
  EvaluatedSdiFin(gamma, delta, iFinM, iCinCM, iCinV, CiFin9x1, CiFinCe9x1, iCinCiCinV, CiFiniCe9x1,
      iCV, iFinCeM, detFin, dSdiFin);

  // if cmat != nullptr, we are evaluating the structural residual and linearizations, so we need to
  // calculate the stresses and the cmat if you like to evaluate the off-diagonal block of your
  // monolithic system (structural residual w.r.t. dofs of another field), you need to pass NULL as
  // the cmat when you call Evaluate() in the element
  if (cmat != nullptr)
  {
    // cmat = 2 dS/dC = 2 \frac{\partial S}{\partial C} + 2 \frac{\partial S}{\partial F_{in}^{-1}}
    // : \frac{\partial F_{in}^{-1}}{\partial C} = cmatiso + cmatadd
    EvaluateStressCmatIso(iCV, iCinV, iCinCiCinV, gamma, delta, detFin, *stress, cmatiso);
    cmat->Update(1.0, cmatiso, 0.0);

    // evaluate additional terms for the elasticity tensor
    // cmatadd = 2 \frac{\partial S}{\partial F_{in}^{-1}} : \frac{\partial F_{in}^{-1}}{\partial
    // C}, where F_{in}^{-1} can be multiplicatively composed of several inelastic contributions
    EvaluateAdditionalCmat(defgrad, iCV, dSdiFin, cmatadd);
    cmat->Update(1.0, cmatadd, 1.0);
  }
  // evaluate OD Block
  else
  {
    // get source of deformation for this OD block depending on the differentiation type
    auto source(PAR::InelasticSource::none);
    const int differentiationtype =
        params.get<int>("differentiationtype", static_cast<int>(STR::DifferentiationType::none));
    if (differentiationtype == static_cast<int>(STR::DifferentiationType::elch))
      source = PAR::InelasticSource::concentration;
    else if (differentiationtype == static_cast<int>(STR::DifferentiationType::temp))
      source = PAR::InelasticSource::temperature;
    else
      dserror("unknown scalaratype");

    EvaluateODStiffMat(source, defgrad, dSdiFin, *stress);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::MultiplicativeSplitDefgrad_ElastHyper::EvaluateStressCmatIso(
    const LINALG::Matrix<6, 1>& iCV, const LINALG::Matrix<6, 1>& iCinV,
    const LINALG::Matrix<6, 1>& iCinCiCinV, const LINALG::Matrix<3, 1>& gamma,
    const LINALG::Matrix<8, 1>& delta, const double detFin, LINALG::Matrix<6, 1>& stress,
    LINALG::Matrix<6, 6>& cmatiso) const
{
  // clear variables
  stress.Clear();
  cmatiso.Clear();

  // 2nd Piola Kirchhoff stresses
  stress.Update(gamma(0), iCinV, 1.0);
  stress.Update(gamma(1), iCinCiCinV, 1.0);
  stress.Update(gamma(2), iCV, 1.0);
  stress.Scale(detFin);

  // constitutive tensor
  cmatiso.MultiplyNT(delta(0), iCinV, iCinV, 1.);
  cmatiso.MultiplyNT(delta(1), iCinCiCinV, iCinV, 1.);
  cmatiso.MultiplyNT(delta(1), iCinV, iCinCiCinV, 1.);
  cmatiso.MultiplyNT(delta(2), iCinV, iCV, 1.);
  cmatiso.MultiplyNT(delta(2), iCV, iCinV, 1.);
  cmatiso.MultiplyNT(delta(3), iCinCiCinV, iCinCiCinV, 1.);
  cmatiso.MultiplyNT(delta(4), iCinCiCinV, iCV, 1.);
  cmatiso.MultiplyNT(delta(4), iCV, iCinCiCinV, 1.);
  cmatiso.MultiplyNT(delta(5), iCV, iCV, 1.);
  AddtoCmatHolzapfelProduct(cmatiso, iCV, delta(6));
  AddtoCmatHolzapfelProduct(cmatiso, iCinV, delta(7));
  cmatiso.Scale(detFin);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::MultiplicativeSplitDefgrad_ElastHyper::EvaluateKinQuantElast(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinM,
    LINALG::Matrix<6, 1>& iCinV, LINALG::Matrix<6, 1>& iCinCiCinV, LINALG::Matrix<6, 1>& iCV,
    LINALG::Matrix<3, 3>& iCinCM, LINALG::Matrix<3, 3>& iFinCeM, LINALG::Matrix<9, 1>& CiFin9x1,
    LINALG::Matrix<9, 1>& CiFinCe9x1, LINALG::Matrix<9, 1>& CiFiniCe9x1,
    LINALG::Matrix<3, 1>& prinv) const
{
  // inverse inelastic right Cauchy-Green
  static LINALG::Matrix<3, 3> iCinM(true);
  iCinM.MultiplyNT(1.0, iFinM, iFinM, 0.0);
  UTILS::VOIGT::Stresses::MatrixToVector(iCinM, iCinV);

  // inverse right Cauchy-Green
  static LINALG::Matrix<3, 3> iCM(true);
  static LINALG::Matrix<3, 3> CM(true);
  CM.MultiplyTN(1.0, *defgrad, *defgrad, 0.0);
  iCM.Invert(CM);
  UTILS::VOIGT::Stresses::MatrixToVector(iCM, iCV);

  // C_{in}^{-1} * C * C_{in}^{-1}
  static LINALG::Matrix<3, 3> tmp(true);
  static LINALG::Matrix<3, 3> iCinCiCinM;
  MAT::EvaluateiCinCiCin(CM, iCinM, iCinCiCinM);
  UTILS::VOIGT::Stresses::MatrixToVector(iCinCiCinM, iCinCiCinV);

  // elastic right Cauchy-Green in strain-like Voigt notation.
  static LINALG::Matrix<3, 3> CeM(true);
  MAT::EvaluateCe(*defgrad, iFinM, CeM);
  static LINALG::Matrix<6, 1> CeV_strain(true);
  UTILS::VOIGT::Strains::MatrixToVector(CeM, CeV_strain);

  // principal invariants of elastic right Cauchy-Green strain
  UTILS::VOIGT::Strains::InvariantsPrincipal(prinv, CeV_strain);

  // C_{in}^{-1} * C
  iCinCM.MultiplyNN(1.0, iCinM, CM, 0.0);

  // F_{in}^{-1} * C_e
  iFinCeM.MultiplyNN(1.0, iFinM, CeM, 0.0);

  // C * F_{in}^{-1}
  static LINALG::Matrix<3, 3> CiFinM(true);
  CiFinM.MultiplyNN(1.0, CM, iFinM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(CiFinM, CiFin9x1);

  // C * F_{in}^{-1} * C_e
  static LINALG::Matrix<3, 3> CiFinCeM(true);
  tmp.MultiplyNN(1.0, CM, iFinM, 0.0);
  CiFinCeM.MultiplyNN(1.0, tmp, CeM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(CiFinCeM, CiFinCe9x1);

  // C * F_{in}^{-1} * C_e^{-1}
  static LINALG::Matrix<3, 3> CiFiniCeM(true);
  static LINALG::Matrix<3, 3> iCeM(true);
  iCeM.Invert(CeM);
  tmp.MultiplyNN(1.0, CM, iFinM, 0.0);
  CiFiniCeM.MultiplyNN(1.0, tmp, iCeM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(CiFiniCeM, CiFiniCe9x1);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::MultiplicativeSplitDefgrad_ElastHyper::EvaluateInvariantDerivatives(
    const LINALG::Matrix<3, 1>& prinv, const int gp, const int eleGID, LINALG::Matrix<3, 1>& dPI,
    LINALG::Matrix<6, 1>& ddPII) const
{
  // clear variables
  dPI.Clear();
  ddPII.Clear();

  // loop over map of associated potential summands
  // derivatives of strain energy function w.r.t. principal invariants
  for (const auto& p : potsumel_)
  {
    p->AddDerivativesPrincipal(dPI, ddPII, prinv, gp, eleGID);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::MultiplicativeSplitDefgrad_ElastHyper::EvaluatedSdiFin(const LINALG::Matrix<3, 1>& gamma,
    const LINALG::Matrix<8, 1>& delta, const LINALG::Matrix<3, 3>& iFinM,
    const LINALG::Matrix<3, 3>& iCinCM, const LINALG::Matrix<6, 1>& iCinV,
    const LINALG::Matrix<9, 1>& CiFin9x1, const LINALG::Matrix<9, 1>& CiFinCe9x1,
    const LINALG::Matrix<6, 1>& iCinCiCinV, const LINALG::Matrix<9, 1>& CiFiniCe9x1,
    const LINALG::Matrix<6, 1>& iCV, const LINALG::Matrix<3, 3>& iFinCeM, const double detFin,
    LINALG::Matrix<6, 9>& dSdiFin) const
{
  // clear variable
  dSdiFin.Clear();

  // calculate identity tensor
  static LINALG::Matrix<3, 3> id(true);
  for (int i = 0; i < 3; ++i) id(i, i) = 1.0;

  // derivative of second Piola Kirchhoff stresses w.r.t. inverse growth deformation gradient
  // (contribution from iFin)
  MAT::AddRightNonSymmetricHolzapfelProduct(dSdiFin, id, iFinM, gamma(0));
  MAT::AddRightNonSymmetricHolzapfelProduct(dSdiFin, iCinCM, iFinM, gamma(1));
  dSdiFin.MultiplyNT(delta(0), iCinV, CiFin9x1, 1.);
  dSdiFin.MultiplyNT(delta(1), iCinV, CiFinCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(1), iCinCiCinV, CiFin9x1, 1.);
  dSdiFin.MultiplyNT(delta(2), iCinV, CiFiniCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(2), iCV, CiFin9x1, 1.);
  dSdiFin.MultiplyNT(delta(3), iCinCiCinV, CiFinCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(4), iCinCiCinV, CiFiniCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(4), iCV, CiFinCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(5), iCV, CiFiniCe9x1, 1.);
  MAT::AddRightNonSymmetricHolzapfelProduct(dSdiFin, id, iFinCeM, gamma(1));
  dSdiFin.Scale(detFin);

  // derivative of second Piola Kirchhoff stresses w.r.t. inverse growth deformation gradient
  // (contribution from det(Fin))

  // dS/d(det(Fin))
  LINALG::Matrix<6, 1> dSddetFin(true);
  dSddetFin.Update(gamma(0), iCinV, 0.0);
  dSddetFin.Update(gamma(1), iCinCiCinV, 1.0);
  dSddetFin.Update(gamma(2), iCV, 1.0);

  // d(det(Fin))/diFin
  LINALG::Matrix<9, 1> ddetFindiFinV(true);
  LINALG::Matrix<3, 3> ddetFindiFinM(true);
  LINALG::Matrix<3, 3> FinM(true);
  FinM.Invert(iFinM);
  ddetFindiFinM.UpdateT((-1.0) * detFin, FinM);
  UTILS::VOIGT::Matrix3x3to9x1(ddetFindiFinM, ddetFindiFinV);

  // chain rule to get dS/d(det(Fin)) * d(det(Fin))/diFin
  dSdiFin.MultiplyNT(1.0, dSddetFin, ddetFindiFinV, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::MultiplicativeSplitDefgrad_ElastHyper::EvaluateAdditionalCmat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<6, 1>& iCV,
    const LINALG::Matrix<6, 9>& dSdiFin, LINALG::Matrix<6, 6>& cmatadd)
{
  // clear variable
  cmatadd.Clear();

  const auto& facdefgradin = inelastic_->FacDefGradIn();
  const auto& iFinjM = inelastic_->GetiFinj();
  const int num_contributions = inelastic_->NumInelasticDefGrad();

  // check amount of factors the inelastic deformation gradient consists of and choose
  // implementation accordingly
  if (num_contributions == 1)
  {
    facdefgradin[0].second->EvaluateAdditionalCmat(
        defgrad, iFinjM[0].second, iCV, dSdiFin, cmatadd);
  }
  else if (num_contributions > 1)
  {
    // static variables
    // dSdiFinj = dSdiFin : diFindiFinj
    // diFindiFinj = \Pi_(k=num_contributions_part)^(j+1) iFin_k : dFinjdFinj : \Pi_(k=j-1)^(0). The
    // double contraction and derivative in index notation for three contributions: diFindiFinj_abcd
    // = iFin_(j-1)_ae \delta_ec \delta_df iFin_(j+1)_fb = iFin_(j-1)_ac iFin_(j+1)_db. This is
    // performed by nonsymmetric product \Pi_(k=num_contributions_part)^(j+1) iFin_k (x)
    // \Pi_(k=j-1)^(0) iFin_k, where (x) denots the nonsymmetric product and \Pi is the
    // multiplication operator.
    static LINALG::Matrix<6, 9> dSdiFinj(true);
    static LINALG::Matrix<9, 9> diFindiFinj(true);
    static LINALG::Matrix<3, 3> id(true);
    for (int i = 0; i < 3; ++i) id(i, i) = 1.0;

    // product of all iFinj, except for the one, that is currently evaluated. In case of inner
    // (neither first or last) we have two products
    static LINALG::Matrix<3, 3> producta(true);
    static LINALG::Matrix<3, 3> productb(true);
    static LINALG::Matrix<3, 3> producta_temp(true);
    static LINALG::Matrix<3, 3> productb_temp(true);

    for (int i = 0; i < num_contributions; ++i)
    {
      // clear static variable
      diFindiFinj.Clear();

      // multiply all inelastic deformation gradients except for range between first and current
      producta = id;
      for (int j = num_contributions - 1; j > i; --j)
      {
        producta_temp.Multiply(1.0, producta, iFinjM[j].second, 0.0);
        producta.Update(1.0, producta_temp, 0.0);
      }

      // multiply all inelastic deformation gradients except for range between last and current
      productb = id;
      if (i > 0)
      {
        for (int j = i - 1; j >= 0; --j)
        {
          productb_temp.Multiply(1.0, productb, iFinjM[j].second, 0.0);
          productb.Update(1.0, productb_temp, 0.0);
        }
      }

      // evaluate additional contribution to C by applying chain rule
      AddNonSymmetricProduct(1.0, producta, productb, diFindiFinj);
      dSdiFinj.Multiply(1.0, dSdiFin, diFindiFinj, 0.0);
      facdefgradin[i].second->EvaluateAdditionalCmat(
          defgrad, iFinjM[i].second, iCV, dSdiFinj, cmatadd);
    }
  }
  else
    dserror("You should not be here");
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::MultiplicativeSplitDefgrad_ElastHyper::Setup(
    const int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Read anisotropy
  anisotropy_->SetNumberOfGaussPoints(numgp);
  anisotropy_->ReadAnisotropyFromElement(linedef);

  // elastic materials
  for (const auto& summand : potsumel_) summand->Setup(numgp, linedef);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::MultiplicativeSplitDefgrad_ElastHyper::Update()
{
  // loop map of associated potential summands
  for (const auto& summand : potsumel_) summand->Update();
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::MultiplicativeSplitDefgrad_ElastHyper::EvaluateODStiffMat(PAR::InelasticSource source,
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<6, 9>& dSdiFin,
    LINALG::Matrix<6, 1>& dstressdx)
{
  // clear variable
  dstressdx.Clear();

  // References to vector of inelastic contributions and inelastic deformation gradients
  const auto& facdefgradin = inelastic_->FacDefGradIn();
  const auto& iFinjM = inelastic_->GetiFinj();

  // number of contributions for this source
  const int num_contributions = inelastic_->NumInelasticDefGrad();

  // check number of factors the inelastic deformation gradient consists of and choose
  // implementation accordingly
  if (num_contributions == 1)
  {
    facdefgradin[0].second->EvaluateODStiffMat(defgrad, iFinjM[0].second, dSdiFin, dstressdx);
  }
  else if (num_contributions > 1)
  {
    // static variables
    // dSdiFinj = dSdiFin : diFindiFinj
    // diFindiFinj = \Pi_(k=num_contributions_part)^(j+1) iFin_k : dFinjdFinj : \Pi_(k=j-1)^(0). The
    // double contraction and derivative in index notation for three contributions: diFindiFinj_abcd
    // = iFin_(j-1)_ae \delta_ec \delta_df iFin_(j+1)_fb = iFin_(j-1)_ac iFin_(j+1)_db. This is
    // performed by nonsymmetric product \Pi_(k=num_contributions_part)^(j+1) iFin_k (x)
    // \Pi_(k=j-1)^(0) iFin_k, where (x) denots the nonsymmetric product and \Pi is the
    // multiplication operator.
    static LINALG::Matrix<6, 9> dSdiFinj(true);
    static LINALG::Matrix<9, 9> diFindiFinj(true);
    static LINALG::Matrix<3, 3> id(true);
    for (int i = 0; i < 3; ++i) id(i, i) = 1.0;

    // product of all iFinj, except for the one, that is currently evaluated. In case of inner
    // (neither first or last) we have two products
    static LINALG::Matrix<3, 3> producta(true);
    static LINALG::Matrix<3, 3> productb(true);
    static LINALG::Matrix<3, 3> producta_temp(true);
    static LINALG::Matrix<3, 3> productb_temp(true);

    for (int i = 0; i < num_contributions; ++i)
    {
      // only if the contribution is from this source, the derivative is non-zero
      if (facdefgradin[i].first == source)
      {
        // clear static variable
        diFindiFinj.Clear();

        // multiply all inelastic deformation gradients except for range between first and current
        // in reverse order
        producta = id;
        for (int j = num_contributions - 1; j > i; --j)
        {
          producta_temp.Multiply(1.0, producta, iFinjM[j].second, 0.0);
          producta.Update(1.0, producta_temp, 0.0);
        }

        // multiply all inelastic deformation gradients except for range between last and current
        productb = id;
        if (i > 0)
        {
          for (int j = i - 1; j >= 0; --j)
          {
            productb_temp.Multiply(1.0, productb, iFinjM[j].second, 0.0);
            productb.Update(1.0, productb_temp, 0.0);
          }
        }

        // evaluate additional contribution to OD block by applying chain rule
        AddNonSymmetricProduct(1.0, producta, productb, diFindiFinj);
        dSdiFinj.Multiply(1.0, dSdiFin, diFindiFinj, 0.0);
        facdefgradin[i].second->EvaluateODStiffMat(defgrad, iFinjM[i].second, dSdiFinj, dstressdx);
      }
    }
  }
  else
    dserror("You should not be here");
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::MultiplicativeSplitDefgrad_ElastHyper::PreEvaluate(
    Teuchos::ParameterList& params, const int gp) const
{
  // loop over all inelastic contributions
  for (int p = 0; p < inelastic_->NumInelasticDefGrad(); ++p)
    inelastic_->FacDefGradIn()[p].second->PreEvaluate(params, gp);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticFactorsHandler::Setup(MAT::PAR::MultiplicativeSplitDefgrad_ElastHyper* params)
{
  facdefgradin_.clear();
  iFinj_.clear();

  // get inelastic deformation gradient factors and assign them to their source
  for (int inelastic_matnum : *(params->inel_defgradfacids_))
  {
    auto inelastic_factor = MAT::InelasticDefgradFactors::Factory(inelastic_matnum);
    if (inelastic_factor == Teuchos::null) dserror("Failed to allocate!");
    std::pair<PAR::InelasticSource, Teuchos::RCP<MAT::InelasticDefgradFactors>> temppair(
        inelastic_factor->GetInelasticSource(), inelastic_factor);
    facdefgradin_.push_back(temppair);
  }

  iFinj_.resize(facdefgradin_.size());

  // safety checks
  // get the scatra structure control parameter list
  const auto& ssicontrol = DRT::Problem::Instance()->SSIControlParams();
  if (Teuchos::getIntegralValue<INPAR::SSI::SolutionSchemeOverFields>(ssicontrol, "COUPALGO") ==
      INPAR::SSI::SolutionSchemeOverFields::ssi_Monolithic)
  {
    for (const auto& inelasitc_factor : facdefgradin_)
    {
      const auto materialtype = inelasitc_factor.second->MaterialType();
      if ((materialtype != INPAR::MAT::mfi_lin_scalar_aniso) and
          (materialtype != INPAR::MAT::mfi_lin_scalar_iso) and
          (materialtype != INPAR::MAT::mfi_poly_intercal_frac_aniso) and
          (materialtype != INPAR::MAT::mfi_poly_intercal_frac_iso) and
          (materialtype != INPAR::MAT::mfi_lin_temp_iso))
      {
        dserror(
            "When you use the 'COUPALGO' 'ssi_Monolithic' from the 'SSI CONTROL' section, you need "
            "to use the material 'MAT_InelasticDefgradLinScalarAniso' "
            "'MAT_InelasticDefgradLinScalarIso', 'MAT_InelasticDefgradPolyScalarIso' or "
            "'MAT_InelasticDefgradPolyScalarAniso'!"
            " If you want to use a different material, feel free to implement it! ;-)");
      }
    }
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticFactorsHandler::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* const defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  // temporary variables
  static LINALG::Matrix<3, 3> iFinp(true);
  static LINALG::Matrix<3, 3> iFin_init_store(true);

  // clear variables
  iFinM.Clear();
  iFin_init_store.Clear();

  for (int i = 0; i < 3; ++i) iFin_init_store(i, i) = 1.0;

  for (int i = 0; i < NumInelasticDefGrad(); ++i)
  {
    // clear tmp variable
    iFinp.Clear();

    // calculate inelastic deformation gradient and its inverse
    facdefgradin_[i].second->EvaluateInverseInelasticDefGrad(defgrad, iFinp);

    // store inelastic deformation gradient of p-th inelastic contribution
    iFinj_[i].second = iFinp;

    // update inverse inelastic deformation gradient
    iFinM.Multiply(iFin_init_store, iFinp);

    // store result for next evaluation
    iFin_init_store.Update(1.0, iFinM, 0.0);
  }
}
