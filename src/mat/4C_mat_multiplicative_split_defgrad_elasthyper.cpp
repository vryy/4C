/*----------------------------------------------------------------------*/
/*! \file
\brief evaluation of a generic material whose deformation gradient is modeled to be split
multiplicatively into elastic and inelastic parts

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_multiplicative_split_defgrad_elasthyper.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_ssi.hpp"
#include "4C_mat_anisotropy.hpp"
#include "4C_mat_elasthyper_service.hpp"
#include "4C_mat_inelastic_defgrad_factors.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper_service.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_structure_new_enum_lists.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::MultiplicativeSplitDefgradElastHyper::MultiplicativeSplitDefgradElastHyper(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      nummat_elast_(matdata.parameters.Get<int>("NUMMATEL")),
      matids_elast_(matdata.parameters.Get<std::vector<int>>("MATIDSEL")),
      numfac_inel_(matdata.parameters.Get<int>("NUMFACINEL")),
      inel_defgradfacids_(matdata.parameters.Get<std::vector<int>>("INELDEFGRADFACIDS")),
      density_(matdata.parameters.Get<double>("DENS"))
{
  // check if sizes fit
  if (nummat_elast_ != static_cast<int>(matids_elast_.size()))
    FOUR_C_THROW(
        "number of elastic materials %d does not fit to size of elastic material ID vector %d",
        nummat_elast_, matids_elast_.size());

  if (numfac_inel_ != static_cast<int>(inel_defgradfacids_.size()))
  {
    FOUR_C_THROW(
        "number of inelastic deformation gradient factors %d does not fit to size of inelastic "
        "deformation gradient ID vector %d",
        numfac_inel_, inel_defgradfacids_.size());
  }
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::MultiplicativeSplitDefgradElastHyper::create_material()
{
  return Teuchos::rcp(new Mat::MultiplicativeSplitDefgradElastHyper(this));
}

Mat::MultiplicativeSplitDefgradElastHyperType
    Mat::MultiplicativeSplitDefgradElastHyperType::instance_;

Core::Communication::ParObject* Mat::MultiplicativeSplitDefgradElastHyperType::Create(
    const std::vector<char>& data)
{
  auto* splitdefgrad_elhy = new Mat::MultiplicativeSplitDefgradElastHyper();
  splitdefgrad_elhy->Unpack(data);

  return splitdefgrad_elhy;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::MultiplicativeSplitDefgradElastHyper::MultiplicativeSplitDefgradElastHyper()
    : anisotropy_(Teuchos::rcp(new Mat::Anisotropy())),
      inelastic_(Teuchos::rcp(new Mat::InelasticFactorsHandler())),
      params_(nullptr),
      potsumel_(0)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::MultiplicativeSplitDefgradElastHyper::MultiplicativeSplitDefgradElastHyper(
    Mat::PAR::MultiplicativeSplitDefgradElastHyper* params)
    : anisotropy_(Teuchos::rcp(new Mat::Anisotropy())),
      inelastic_(Teuchos::rcp(new Mat::InelasticFactorsHandler())),
      params_(params),
      potsumel_(0)
{
  // elastic materials
  for (int matid_elastic : params_->matids_elast_)
  {
    auto elastic_summand = Mat::Elastic::Summand::Factory(matid_elastic);
    if (elastic_summand == Teuchos::null) FOUR_C_THROW("Failed to allocate");
    potsumel_.push_back(elastic_summand);
    elastic_summand->register_anisotropy_extensions(*anisotropy_);
  }

  inelastic_->Setup(params);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  anisotropy_->pack_anisotropy(data);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (const auto& p : potsumel_) p->PackSummand(data);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  params_ = nullptr;
  potsumel_.clear();

  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      auto* mat = Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<Mat::PAR::MultiplicativeSplitDefgradElastHyper*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  anisotropy_->unpack_anisotropy(data, position);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // elastic materials
    for (const auto& matid_elastic : params_->matids_elast_)
    {
      auto elastic_summand = Mat::Elastic::Summand::Factory(matid_elastic);
      if (elastic_summand == Teuchos::null) FOUR_C_THROW("Failed to allocate");
      potsumel_.push_back(elastic_summand);
    }
    // loop map of associated potential summands
    for (const auto& elastic_summand : potsumel_)
    {
      elastic_summand->UnpackSummand(data, position);
      elastic_summand->register_anisotropy_extensions(*anisotropy_);
    }

    // inelastic deformation gradient factors
    inelastic_->Setup(params_);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::Evaluate(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<6, 1>* glstrain,
    Teuchos::ParameterList& params, Core::LinAlg::Matrix<6, 1>* stress,
    Core::LinAlg::Matrix<6, 6>* cmat, const int gp, const int eleGID)
{
  // do all stuff that only has to be done once per Evaluate() call
  pre_evaluate(params, gp);

  // static variables
  static Core::LinAlg::Matrix<6, 6> cmatiso(true);
  static Core::LinAlg::Matrix<6, 9> dSdiFin(true);
  static Core::LinAlg::Matrix<6, 6> cmatadd(true);

  // build inverse inelastic deformation gradient
  static Core::LinAlg::Matrix<3, 3> iFinM(true);
  inelastic_->evaluate_inverse_inelastic_def_grad(defgrad, iFinM);

  // determinante of inelastic deformation gradient
  const double detFin = 1.0 / iFinM.Determinant();

  // static variables of kinetic quantities
  static Core::LinAlg::Matrix<6, 1> iCV(true);
  static Core::LinAlg::Matrix<6, 1> iCinV(true);
  static Core::LinAlg::Matrix<6, 1> iCinCiCinV(true);
  static Core::LinAlg::Matrix<3, 3> iCinCM(true);
  static Core::LinAlg::Matrix<3, 3> iFinCeM(true);
  static Core::LinAlg::Matrix<9, 1> CiFin9x1(true);
  static Core::LinAlg::Matrix<9, 1> CiFinCe9x1(true);
  static Core::LinAlg::Matrix<9, 1> CiFiniCe9x1(true);
  static Core::LinAlg::Matrix<3, 1> prinv(true);
  evaluate_kin_quant_elast(defgrad, iFinM, iCinV, iCinCiCinV, iCV, iCinCM, iFinCeM, CiFin9x1,
      CiFinCe9x1, CiFiniCe9x1, prinv);

  // derivatives of principle invariants
  static Core::LinAlg::Matrix<3, 1> dPIe(true);
  static Core::LinAlg::Matrix<6, 1> ddPIIe(true);
  evaluate_invariant_derivatives(prinv, gp, eleGID, dPIe, ddPIIe);

  // 2nd Piola Kirchhoff stresses factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
  static Core::LinAlg::Matrix<3, 1> gamma(true);
  // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
  static Core::LinAlg::Matrix<8, 1> delta(true);
  // compose coefficients
  CalculateGammaDelta(gamma, delta, prinv, dPIe, ddPIIe);

  // evaluate dSdiFin
  EvaluatedSdiFin(gamma, delta, iFinM, iCinCM, iCinV, CiFin9x1, CiFinCe9x1, iCinCiCinV, CiFiniCe9x1,
      iCV, iFinCeM, detFin, dSdiFin);

  // if cmat != nullptr, we are evaluating the structural residual and linearizations, so we need to
  // calculate the stresses and the cmat if you like to evaluate the off-diagonal block of your
  // monolithic system (structural residual w.r.t. dofs of another field), you need to pass nullptr
  // as the cmat when you call Evaluate() in the element
  if (cmat != nullptr)
  {
    // cmat = 2 dS/dC = 2 \frac{\partial S}{\partial C} + 2 \frac{\partial S}{\partial F_{in}^{-1}}
    // : \frac{\partial F_{in}^{-1}}{\partial C} = cmatiso + cmatadd
    evaluate_stress_cmat_iso(iCV, iCinV, iCinCiCinV, gamma, delta, detFin, *stress, cmatiso);
    cmat->Update(1.0, cmatiso, 0.0);

    // evaluate additional terms for the elasticity tensor
    // cmatadd = 2 \frac{\partial S}{\partial F_{in}^{-1}} : \frac{\partial F_{in}^{-1}}{\partial
    // C}, where F_{in}^{-1} can be multiplicatively composed of several inelastic contributions
    evaluate_additional_cmat(defgrad, iCV, dSdiFin, cmatadd);
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
      FOUR_C_THROW("unknown scalaratype");

    EvaluateODStiffMat(source, defgrad, dSdiFin, *stress);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::evaluate_cauchy_n_dir_and_derivatives(
    const Core::LinAlg::Matrix<3, 3>& defgrd, const Core::LinAlg::Matrix<3, 1>& n,
    const Core::LinAlg::Matrix<3, 1>& dir, double& cauchy_n_dir,
    Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn, Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir,
    Core::LinAlg::Matrix<9, 1>* d_cauchyndir_dF, Core::LinAlg::Matrix<9, 9>* d2_cauchyndir_dF2,
    Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_dn,
    Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_ddir, int gp, int eleGID,
    const double* concentration, const double* temp, double* d_cauchyndir_dT,
    Core::LinAlg::Matrix<9, 1>* d2_cauchyndir_dF_dT)
{
  if (concentration != nullptr) SetConcentrationGP(*concentration);

  // reset sigma contracted with n and dir
  cauchy_n_dir = 0.0;

  static Core::LinAlg::Matrix<6, 1> idV(true);
  for (int i = 0; i < 3; ++i) idV(i) = 1.0;
  static Core::LinAlg::Matrix<3, 3> idM(true);
  for (int i = 0; i < 3; ++i) idM(i, i) = 1.0;
  static Core::LinAlg::Matrix<3, 3> iFinM(true);
  inelastic_->evaluate_inverse_inelastic_def_grad(&defgrd, iFinM);
  static Core::LinAlg::Matrix<3, 3> FeM(true);
  FeM.MultiplyNN(1.0, defgrd, iFinM, 0.0);

  // get elastic left cauchy-green tensor and corresponding principal invariants
  static Core::LinAlg::Matrix<3, 3> beM(true);
  beM.MultiplyNT(1.0, FeM, FeM, 0.0);
  static Core::LinAlg::Matrix<6, 1> beV_strain(true);
  Core::LinAlg::Voigt::Strains::matrix_to_vector(beM, beV_strain);
  static Core::LinAlg::Matrix<3, 1> prinv(true);
  Core::LinAlg::Voigt::Strains::invariants_principal(prinv, beV_strain);
  static Core::LinAlg::Matrix<6, 1> beV_stress(true);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(beM, beV_stress);

  static Core::LinAlg::Matrix<3, 1> beMdn(true);
  beMdn.Multiply(1.0, beM, n, 0.0);
  const double beMdnddir = beMdn.Dot(dir);
  static Core::LinAlg::Matrix<3, 1> beMddir(true);
  beMddir.Multiply(1.0, beM, dir, 0.0);

  static Core::LinAlg::Matrix<3, 3> ibeM(true);
  ibeM.Invert(beM);
  static Core::LinAlg::Matrix<6, 1> ibeV_stress(true);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(ibeM, ibeV_stress);
  static Core::LinAlg::Matrix<3, 1> ibeMdn(true);
  ibeMdn.Multiply(1.0, ibeM, n, 0.0);
  const double ibeMdnddir = ibeMdn.Dot(dir);
  static Core::LinAlg::Matrix<3, 1> ibeMddir(true);
  ibeMddir.Multiply(1.0, ibeM, dir, 0.0);

  // derivatives of principle invariants of elastic left cauchy-green tensor
  static Core::LinAlg::Matrix<3, 1> dPI(true);
  static Core::LinAlg::Matrix<6, 1> ddPII(true);
  evaluate_invariant_derivatives(prinv, gp, eleGID, dPI, ddPII);

  const double detFe = FeM.Determinant();
  const double nddir = n.Dot(dir);
  const double prefac = 2.0 / detFe;

  // calculate \mat{\sigma} \cdot \vec{n} \cdot \vec{v}
  cauchy_n_dir = prefac * (prinv(1) * dPI(1) * nddir + prinv(2) * dPI(2) * nddir +
                              dPI(0) * beMdnddir - prinv(2) * dPI(1) * ibeMdnddir);

  if (d_cauchyndir_dn)
  {
    d_cauchyndir_dn->Update(prinv(1) * dPI(1) + prinv(2) * dPI(2), dir, 0.0);
    d_cauchyndir_dn->Update(dPI(0), beMddir, 1.0);
    d_cauchyndir_dn->Update(-prinv(2) * dPI(1), ibeMddir, 1.0);
    d_cauchyndir_dn->Scale(prefac);
  }

  if (d_cauchyndir_ddir)
  {
    d_cauchyndir_ddir->Update(prinv(1) * dPI(1) + prinv(2) * dPI(2), n, 0.0);
    d_cauchyndir_ddir->Update(dPI(0), beMdn, 1.0);
    d_cauchyndir_ddir->Update(-prinv(2) * dPI(1), ibeMdn, 1.0);
    d_cauchyndir_ddir->Scale(prefac);
  }

  if (d_cauchyndir_dF)
  {
    static Core::LinAlg::Matrix<6, 1> d_I1_dbe(true);
    d_I1_dbe = idV;
    static Core::LinAlg::Matrix<6, 1> d_I2_dbe(true);
    d_I2_dbe.Update(prinv(0), idV, -1.0, beV_stress);
    static Core::LinAlg::Matrix<6, 1> d_I3_dbe(true);
    d_I3_dbe.Update(prinv(2), ibeV_stress, 0.0);

    // calculation of \partial b_{el} / \partial F (elastic left cauchy-green w.r.t. deformation
    // gradient)
    static Core::LinAlg::Matrix<6, 9> d_be_dFe(true);
    d_be_dFe.Clear();
    AddRightNonSymmetricHolzapfelProductStrainLike(d_be_dFe, idM, FeM, 1.0);
    static Core::LinAlg::Matrix<9, 9> d_Fe_dF(true);
    d_Fe_dF.Clear();
    AddNonSymmetricProduct(1.0, idM, iFinM, d_Fe_dF);
    static Core::LinAlg::Matrix<6, 9> d_be_dF(true);
    d_be_dF.Multiply(1.0, d_be_dFe, d_Fe_dF, 0.0);

    // calculation of \partial I_i / \partial F (Invariants of b_{el} w.r.t. deformation gradient)
    static Core::LinAlg::Matrix<9, 1> d_I1_dF(true);
    static Core::LinAlg::Matrix<9, 1> d_I2_dF(true);
    static Core::LinAlg::Matrix<9, 1> d_I3_dF(true);
    d_I1_dF.MultiplyTN(1.0, d_be_dF, d_I1_dbe, 0.0);
    d_I2_dF.MultiplyTN(1.0, d_be_dF, d_I2_dbe, 0.0);
    d_I3_dF.MultiplyTN(1.0, d_be_dF, d_I3_dbe, 0.0);

    // add d_cauchyndir_dI1 \odot d_I1_dF and clear static matrix
    d_cauchyndir_dF->Update(prefac * (prinv(1) * ddPII(5) * nddir + prinv(2) * ddPII(4) * nddir +
                                         ddPII(0) * beMdnddir - prinv(2) * ddPII(5) * ibeMdnddir),
        d_I1_dF, 0.0);
    // add d_cauchyndir_dI2 \odot d_I2_dF
    d_cauchyndir_dF->Update(
        prefac * (dPI(1) * nddir + prinv(1) * ddPII(1) * nddir + prinv(2) * ddPII(3) * nddir +
                     ddPII(5) * beMdnddir - prinv(2) * ddPII(1) * ibeMdnddir),
        d_I2_dF, 1.0);
    // add d_cauchyndir_dI3 \odot d_I3_dF
    d_cauchyndir_dF->Update(
        prefac * (prinv(1) * ddPII(3) * nddir + dPI(2) * nddir + prinv(2) * ddPII(2) * nddir +
                     ddPII(4) * beMdnddir - dPI(1) * ibeMdnddir - prinv(2) * ddPII(3) * ibeMdnddir),
        d_I3_dF, 1.0);

    // next three updates add partial derivative of snt w.r.t. the deformation gradient F for
    // constant invariants first part is term arising from \partial Je^{-1} / \partial F
    static Core::LinAlg::Matrix<3, 3> iFeM(true);
    static Core::LinAlg::Matrix<3, 3> iFeTM(true);
    iFeM.Invert(FeM);
    iFeTM.UpdateT(1.0, iFeM, 0.0);
    static Core::LinAlg::Matrix<9, 1> iFeTV(true);
    Core::LinAlg::Voigt::matrix_3x3_to_9x1(iFeTM, iFeTV);
    static Core::LinAlg::Matrix<1, 9> d_iJe_dFV(true);
    d_iJe_dFV.MultiplyTN(1.0, iFeTV, d_Fe_dF, 0.0);
    d_cauchyndir_dF->UpdateT(-cauchy_n_dir, d_iJe_dFV, 1.0);

    // second part is term arising from \partial b_el * n * v / \partial F
    static Core::LinAlg::Matrix<3, 3> FeMiFinTM(true);
    FeMiFinTM.MultiplyNT(1.0, FeM, iFinM, 0.0);
    static Core::LinAlg::Matrix<3, 1> tempvec(true);
    tempvec.MultiplyTN(1.0, FeMiFinTM, n, 0.0);
    static Core::LinAlg::Matrix<3, 3> d_bednddir_dF(true);
    d_bednddir_dF.MultiplyNT(1.0, dir, tempvec, 0.0);
    // now reuse tempvec
    tempvec.MultiplyTN(1.0, FeMiFinTM, dir, 0.0);
    d_bednddir_dF.MultiplyNT(1.0, n, tempvec, 1.0);
    static Core::LinAlg::Matrix<9, 1> d_bednddir_dFV(true);
    Core::LinAlg::Voigt::matrix_3x3_to_9x1(d_bednddir_dF, d_bednddir_dFV);
    d_cauchyndir_dF->Update(prefac * dPI(0), d_bednddir_dFV, 1.0);

    // third part is term arising from \partial b_el^{-1} * n * v / \partial F
    static Core::LinAlg::Matrix<3, 3> iFM(true);
    iFM.Invert(defgrd);
    static Core::LinAlg::Matrix<3, 1> tempvec2(true);
    tempvec.Multiply(1.0, ibeM, dir, 0.0);
    tempvec2.Multiply(1.0, iFM, n, 0.0);
    static Core::LinAlg::Matrix<3, 3> d_ibednddir_dFM(true);
    d_ibednddir_dFM.MultiplyNT(1.0, tempvec, tempvec2, 0.0);
    // now reuse tempvecs
    tempvec.Multiply(1.0, ibeM, n, 0.0);
    tempvec2.Multiply(1.0, iFM, dir, 0.0);
    d_ibednddir_dFM.MultiplyNT(1.0, tempvec, tempvec2, 1.0);
    d_ibednddir_dFM.Scale(-1.0);
    static Core::LinAlg::Matrix<9, 1> d_ibednddir_dFV(true);
    Core::LinAlg::Voigt::matrix_3x3_to_9x1(d_ibednddir_dFM, d_ibednddir_dFV);
    d_cauchyndir_dF->Update(-prefac * prinv(2) * dPI(1), d_ibednddir_dFV, 1.0);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::evaluate_linearization_od(
    const Core::LinAlg::Matrix<3, 3>& defgrd, const double concentration,
    Core::LinAlg::Matrix<9, 1>* d_F_dx)
{
  SetConcentrationGP(concentration);

  // References to vector of inelastic contributions and inelastic deformation gradients
  auto facdefgradin = inelastic_->FacDefGradIn();

  // number of contributions for this source
  const int num_contributions = inelastic_->NumInelasticDefGrad();

  // build inverse inelastic deformation gradient
  static Core::LinAlg::Matrix<3, 3> iFinM(true);
  inelastic_->evaluate_inverse_inelastic_def_grad(&defgrd, iFinM);

  static Core::LinAlg::Matrix<3, 3> idM(true);
  for (int i = 0; i < 3; ++i) idM(i, i) = 1.0;
  static Core::LinAlg::Matrix<3, 3> FeM(true);
  FeM.MultiplyNN(1.0, defgrd, iFinM, 0.0);

  // calculate the derivative of the deformation gradient w.r.t. the inelastic deformation gradient
  static Core::LinAlg::Matrix<9, 9> d_F_dFin(true);
  d_F_dFin.Clear();
  AddNonSymmetricProduct(1.0, FeM, idM, d_F_dFin);

  static Core::LinAlg::Matrix<9, 1> d_Fin_dx(true);

  // check number of factors the inelastic deformation gradient consists of and choose
  // implementation accordingly
  if (num_contributions == 1)
  {
    facdefgradin[0].second->evaluate_inelastic_def_grad_derivative(defgrd.Determinant(), d_Fin_dx);
  }
  else
    FOUR_C_THROW("NOT YET IMPLEMENTED");

  d_F_dx->MultiplyNN(1.0, d_F_dFin, d_Fin_dx, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::evaluate_stress_cmat_iso(
    const Core::LinAlg::Matrix<6, 1>& iCV, const Core::LinAlg::Matrix<6, 1>& iCinV,
    const Core::LinAlg::Matrix<6, 1>& iCinCiCinV, const Core::LinAlg::Matrix<3, 1>& gamma,
    const Core::LinAlg::Matrix<8, 1>& delta, const double detFin,
    Core::LinAlg::Matrix<6, 1>& stress, Core::LinAlg::Matrix<6, 6>& cmatiso) const
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
void Mat::MultiplicativeSplitDefgradElastHyper::evaluate_kin_quant_elast(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<3, 3>& iFinM,
    Core::LinAlg::Matrix<6, 1>& iCinV, Core::LinAlg::Matrix<6, 1>& iCinCiCinV,
    Core::LinAlg::Matrix<6, 1>& iCV, Core::LinAlg::Matrix<3, 3>& iCinCM,
    Core::LinAlg::Matrix<3, 3>& iFinCeM, Core::LinAlg::Matrix<9, 1>& CiFin9x1,
    Core::LinAlg::Matrix<9, 1>& CiFinCe9x1, Core::LinAlg::Matrix<9, 1>& CiFiniCe9x1,
    Core::LinAlg::Matrix<3, 1>& prinv) const
{
  // inverse inelastic right Cauchy-Green
  static Core::LinAlg::Matrix<3, 3> iCinM(true);
  iCinM.MultiplyNT(1.0, iFinM, iFinM, 0.0);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(iCinM, iCinV);

  // inverse right Cauchy-Green
  static Core::LinAlg::Matrix<3, 3> iCM(true);
  static Core::LinAlg::Matrix<3, 3> CM(true);
  CM.MultiplyTN(1.0, *defgrad, *defgrad, 0.0);
  iCM.Invert(CM);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(iCM, iCV);

  // C_{in}^{-1} * C * C_{in}^{-1}
  static Core::LinAlg::Matrix<3, 3> tmp(true);
  static Core::LinAlg::Matrix<3, 3> iCinCiCinM;
  Mat::EvaluateiCinCiCin(CM, iCinM, iCinCiCinM);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(iCinCiCinM, iCinCiCinV);

  // elastic right Cauchy-Green in strain-like Voigt notation.
  static Core::LinAlg::Matrix<3, 3> CeM(true);
  Mat::EvaluateCe(*defgrad, iFinM, CeM);
  static Core::LinAlg::Matrix<6, 1> CeV_strain(true);
  Core::LinAlg::Voigt::Strains::matrix_to_vector(CeM, CeV_strain);

  // principal invariants of elastic right Cauchy-Green strain
  Core::LinAlg::Voigt::Strains::invariants_principal(prinv, CeV_strain);

  // C_{in}^{-1} * C
  iCinCM.MultiplyNN(1.0, iCinM, CM, 0.0);

  // F_{in}^{-1} * C_e
  iFinCeM.MultiplyNN(1.0, iFinM, CeM, 0.0);

  // C * F_{in}^{-1}
  static Core::LinAlg::Matrix<3, 3> CiFinM(true);
  CiFinM.MultiplyNN(1.0, CM, iFinM, 0.0);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(CiFinM, CiFin9x1);

  // C * F_{in}^{-1} * C_e
  static Core::LinAlg::Matrix<3, 3> CiFinCeM(true);
  tmp.MultiplyNN(1.0, CM, iFinM, 0.0);
  CiFinCeM.MultiplyNN(1.0, tmp, CeM, 0.0);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(CiFinCeM, CiFinCe9x1);

  // C * F_{in}^{-1} * C_e^{-1}
  static Core::LinAlg::Matrix<3, 3> CiFiniCeM(true);
  static Core::LinAlg::Matrix<3, 3> iCeM(true);
  iCeM.Invert(CeM);
  tmp.MultiplyNN(1.0, CM, iFinM, 0.0);
  CiFiniCeM.MultiplyNN(1.0, tmp, iCeM, 0.0);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(CiFiniCeM, CiFiniCe9x1);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::evaluate_invariant_derivatives(
    const Core::LinAlg::Matrix<3, 1>& prinv, const int gp, const int eleGID,
    Core::LinAlg::Matrix<3, 1>& dPI, Core::LinAlg::Matrix<6, 1>& ddPII) const
{
  // clear variables
  dPI.Clear();
  ddPII.Clear();

  // loop over map of associated potential summands
  // derivatives of strain energy function w.r.t. principal invariants
  for (const auto& p : potsumel_)
  {
    p->add_derivatives_principal(dPI, ddPII, prinv, gp, eleGID);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::EvaluatedSdiFin(
    const Core::LinAlg::Matrix<3, 1>& gamma, const Core::LinAlg::Matrix<8, 1>& delta,
    const Core::LinAlg::Matrix<3, 3>& iFinM, const Core::LinAlg::Matrix<3, 3>& iCinCM,
    const Core::LinAlg::Matrix<6, 1>& iCinV, const Core::LinAlg::Matrix<9, 1>& CiFin9x1,
    const Core::LinAlg::Matrix<9, 1>& CiFinCe9x1, const Core::LinAlg::Matrix<6, 1>& iCinCiCinV,
    const Core::LinAlg::Matrix<9, 1>& CiFiniCe9x1, const Core::LinAlg::Matrix<6, 1>& iCV,
    const Core::LinAlg::Matrix<3, 3>& iFinCeM, const double detFin,
    Core::LinAlg::Matrix<6, 9>& dSdiFin) const
{
  // clear variable
  dSdiFin.Clear();

  // calculate identity tensor
  static Core::LinAlg::Matrix<3, 3> id(true);
  for (int i = 0; i < 3; ++i) id(i, i) = 1.0;

  // derivative of second Piola Kirchhoff stresses w.r.t. inverse growth deformation gradient
  // (contribution from iFin)
  Mat::AddRightNonSymmetricHolzapfelProduct(dSdiFin, id, iFinM, gamma(0));
  Mat::AddRightNonSymmetricHolzapfelProduct(dSdiFin, iCinCM, iFinM, gamma(1));
  dSdiFin.MultiplyNT(delta(0), iCinV, CiFin9x1, 1.);
  dSdiFin.MultiplyNT(delta(1), iCinV, CiFinCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(1), iCinCiCinV, CiFin9x1, 1.);
  dSdiFin.MultiplyNT(delta(2), iCinV, CiFiniCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(2), iCV, CiFin9x1, 1.);
  dSdiFin.MultiplyNT(delta(3), iCinCiCinV, CiFinCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(4), iCinCiCinV, CiFiniCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(4), iCV, CiFinCe9x1, 1.);
  dSdiFin.MultiplyNT(delta(5), iCV, CiFiniCe9x1, 1.);
  Mat::AddRightNonSymmetricHolzapfelProduct(dSdiFin, id, iFinCeM, gamma(1));
  dSdiFin.Scale(detFin);

  // derivative of second Piola Kirchhoff stresses w.r.t. inverse growth deformation gradient
  // (contribution from det(Fin))

  // dS/d(det(Fin))
  Core::LinAlg::Matrix<6, 1> dSddetFin(true);
  dSddetFin.Update(gamma(0), iCinV, 0.0);
  dSddetFin.Update(gamma(1), iCinCiCinV, 1.0);
  dSddetFin.Update(gamma(2), iCV, 1.0);

  // d(det(Fin))/diFin
  Core::LinAlg::Matrix<9, 1> ddetFindiFinV(true);
  Core::LinAlg::Matrix<3, 3> ddetFindiFinM(true);
  Core::LinAlg::Matrix<3, 3> FinM(true);
  FinM.Invert(iFinM);
  ddetFindiFinM.UpdateT((-1.0) * detFin, FinM);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(ddetFindiFinM, ddetFindiFinV);

  // chain rule to get dS/d(det(Fin)) * d(det(Fin))/diFin
  dSdiFin.MultiplyNT(1.0, dSddetFin, ddetFindiFinV, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::evaluate_additional_cmat(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<6, 1>& iCV,
    const Core::LinAlg::Matrix<6, 9>& dSdiFin, Core::LinAlg::Matrix<6, 6>& cmatadd)
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
    facdefgradin[0].second->evaluate_additional_cmat(
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
    static Core::LinAlg::Matrix<6, 9> dSdiFinj(true);
    static Core::LinAlg::Matrix<9, 9> diFindiFinj(true);
    static Core::LinAlg::Matrix<3, 3> id(true);
    for (int i = 0; i < 3; ++i) id(i, i) = 1.0;

    // product of all iFinj, except for the one, that is currently evaluated. In case of inner
    // (neither first or last) we have two products
    static Core::LinAlg::Matrix<3, 3> producta(true);
    static Core::LinAlg::Matrix<3, 3> productb(true);
    static Core::LinAlg::Matrix<3, 3> producta_temp(true);
    static Core::LinAlg::Matrix<3, 3> productb_temp(true);

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
      facdefgradin[i].second->evaluate_additional_cmat(
          defgrad, iFinjM[i].second, iCV, dSdiFinj, cmatadd);
    }
  }
  else
    FOUR_C_THROW("You should not be here");
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::Setup(
    const int numgp, Input::LineDefinition* linedef)
{
  // Read anisotropy
  anisotropy_->set_number_of_gauss_points(numgp);
  anisotropy_->read_anisotropy_from_element(linedef);

  // elastic materials
  for (const auto& summand : potsumel_) summand->Setup(numgp, linedef);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::Update()
{
  // loop map of associated potential summands
  for (const auto& summand : potsumel_) summand->Update();
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::EvaluateODStiffMat(PAR::InelasticSource source,
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<6, 9>& dSdiFin,
    Core::LinAlg::Matrix<6, 1>& dstressdx)
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
    static Core::LinAlg::Matrix<6, 9> dSdiFinj(true);
    static Core::LinAlg::Matrix<9, 9> diFindiFinj(true);
    static Core::LinAlg::Matrix<3, 3> id(true);
    for (int i = 0; i < 3; ++i) id(i, i) = 1.0;

    // product of all iFinj, except for the one, that is currently evaluated. In case of inner
    // (neither first or last) we have two products
    static Core::LinAlg::Matrix<3, 3> producta(true);
    static Core::LinAlg::Matrix<3, 3> productb(true);
    static Core::LinAlg::Matrix<3, 3> producta_temp(true);
    static Core::LinAlg::Matrix<3, 3> productb_temp(true);

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
    FOUR_C_THROW("You should not be here");
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::pre_evaluate(
    Teuchos::ParameterList& params, const int gp) const
{
  // loop over all inelastic contributions
  for (int p = 0; p < inelastic_->NumInelasticDefGrad(); ++p)
    inelastic_->FacDefGradIn()[p].second->pre_evaluate(params, gp);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::SetConcentrationGP(const double concentration)
{
  for (int p = 0; p < inelastic_->NumInelasticDefGrad(); ++p)
    inelastic_->FacDefGradIn()[p].second->SetConcentrationGP(concentration);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticFactorsHandler::Setup(Mat::PAR::MultiplicativeSplitDefgradElastHyper* params)
{
  facdefgradin_.clear();
  i_finj_.clear();

  // get inelastic deformation gradient factors and assign them to their source
  for (int inelastic_matnum : params->inel_defgradfacids_)
  {
    auto inelastic_factor = Mat::InelasticDefgradFactors::Factory(inelastic_matnum);
    if (inelastic_factor == Teuchos::null) FOUR_C_THROW("Failed to allocate!");
    std::pair<PAR::InelasticSource, Teuchos::RCP<Mat::InelasticDefgradFactors>> temppair(
        inelastic_factor->GetInelasticSource(), inelastic_factor);
    facdefgradin_.push_back(temppair);
  }

  i_finj_.resize(facdefgradin_.size());

  // safety checks
  // get the scatra structure control parameter list
  const auto& ssicontrol = Global::Problem::Instance()->SSIControlParams();
  if (Teuchos::getIntegralValue<Inpar::SSI::SolutionSchemeOverFields>(ssicontrol, "COUPALGO") ==
      Inpar::SSI::SolutionSchemeOverFields::ssi_Monolithic)
  {
    for (const auto& inelasitc_factor : facdefgradin_)
    {
      const auto materialtype = inelasitc_factor.second->MaterialType();
      if ((materialtype != Core::Materials::mfi_lin_scalar_aniso) and
          (materialtype != Core::Materials::mfi_lin_scalar_iso) and
          (materialtype != Core::Materials::mfi_lin_temp_iso) and
          (materialtype != Core::Materials::mfi_no_growth) and
          (materialtype != Core::Materials::mfi_time_funct) and
          (materialtype != Core::Materials::mfi_poly_intercal_frac_aniso) and
          (materialtype != Core::Materials::mfi_poly_intercal_frac_iso))
      {
        FOUR_C_THROW(
            "When you use the 'COUPALGO' 'ssi_Monolithic' from the 'SSI CONTROL' section, you need "
            "to use one of the materials derived from 'Mat::InelasticDefgradFactors'!"
            " If you want to use a different material, feel free to implement it! ;-)");
      }
    }
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticFactorsHandler::evaluate_inverse_inelastic_def_grad(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, Core::LinAlg::Matrix<3, 3>& iFinM)
{
  // temporary variables
  static Core::LinAlg::Matrix<3, 3> iFinp(true);
  static Core::LinAlg::Matrix<3, 3> iFin_init_store(true);

  // clear variables
  iFinM.Clear();
  iFin_init_store.Clear();

  for (int i = 0; i < 3; ++i) iFin_init_store(i, i) = 1.0;

  for (int i = 0; i < NumInelasticDefGrad(); ++i)
  {
    // clear tmp variable
    iFinp.Clear();

    // calculate inelastic deformation gradient and its inverse
    facdefgradin_[i].second->evaluate_inverse_inelastic_def_grad(defgrad, iFinp);

    // store inelastic deformation gradient of p-th inelastic contribution
    i_finj_[i].second = iFinp;

    // update inverse inelastic deformation gradient
    iFinM.Multiply(iFin_init_store, iFinp);

    // store result for next evaluation
    iFin_init_store.Update(1.0, iFinM, 0.0);
  }
}

FOUR_C_NAMESPACE_CLOSE
