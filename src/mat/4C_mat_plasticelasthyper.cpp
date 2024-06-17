/*----------------------------------------------------------------------*/
/*! \file
\brief This file contains the hyperelastic toolbox with application to finite
strain plasticity using a semi-smooth Newton method. It allows summing up
several summands of isotropic non-splitted type to build
a hyperelastic strain energy function.

The input line should read
MAT 1 MAT_PlasticElastHyper NUMMAT 1 MATIDS 2 DENS 1.0 INITYIELD 0.45 ISOHARD 0.12924
EXPISOHARD 16.93 INFYIELD 0.715 KINHARD 0.0 CTE 1.0e-5 INITTEMP 293 YIELDSOFT 0.002 HARDSOFT 0.002
VISC 1e-4 VISC_TEMP 0.003 PL_SPIN_CHI -50 rY_11 1.0 rY_22 0.9 rY_33 0.9 rY_12 0.7 rY_23 0.57385
rY_13 0.7

\level 2
*/

/*----------------------------------------------------------------------*/

#include "4C_mat_plasticelasthyper.hpp"

#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::PlasticElastHyper::PlasticElastHyper(const Core::Mat::PAR::Parameter::Data& matdata)
    : Mat::PAR::ElastHyper(matdata),
      inityield_(matdata.parameters.Get<double>("INITYIELD")),
      isohard_(matdata.parameters.Get<double>("ISOHARD")),
      expisohard_(matdata.parameters.Get<double>("EXPISOHARD")),
      infyield_(matdata.parameters.Get<double>("INFYIELD")),
      kinhard_(matdata.parameters.Get<double>("KINHARD")),
      visc_(matdata.parameters.Get<double>("VISC")),
      rate_dependency_(matdata.parameters.Get<double>("RATE_DEPENDENCY")),
      visc_soft_(matdata.parameters.Get<double>("VISC_SOFT")),
      cte_(matdata.parameters.Get<double>("CTE")),
      inittemp_(matdata.parameters.Get<double>("INITTEMP")),
      yieldsoft_(matdata.parameters.Get<double>("YIELDSOFT")),
      hardsoft_(matdata.parameters.Get<double>("HARDSOFT")),
      taylor_quinney_(matdata.parameters.Get<double>("TAYLOR_QUINNEY")),
      plspin_chi_(-1. * matdata.parameters.Get<double>("PL_SPIN_CHI")),
      rY_11_(matdata.parameters.Get<double>("rY_11")),
      rY_22_(matdata.parameters.Get<double>("rY_22")),
      rY_33_(matdata.parameters.Get<double>("rY_33")),
      rY_12_(matdata.parameters.Get<double>("rY_12")),
      rY_23_(matdata.parameters.Get<double>("rY_23")),
      rY_13_(matdata.parameters.Get<double>("rY_13")),
      cpl_(0.),
      stab_s_(0.),
      dis_mode_(Inpar::TSI::pl_multiplier)
{
  // check if sizes fit
  if (nummat_ != (int)matids_.size())
    FOUR_C_THROW("number of materials %d does not fit to size of material vector %d", nummat_,
        matids_.size());

  // check plastic parameter validity
  if (inityield_ <= 0.)
    std::cout << "***************************************************\n"
                 "Warning: negative initial yield stress detected!!!\n"
                 "this will be converted to a purely elastic response"
              << "\n***************************************************\n"
              << std::endl;

  // no infyield provided 0. is default
  if (infyield_ == 0.)
    if (expisohard_ != 0.) FOUR_C_THROW("hardening exponent provided without inf yield stress");
  if (expisohard_ < 0.) FOUR_C_THROW("Nonlinear hardening exponent must be non-negative");

  // polyconvexity check is just implemented for isotropic hyperlastic materials
  if (polyconvex_)
    FOUR_C_THROW(
        "This polyconvexity-check is just implemented for isotropic "
        "hyperelastic-materials (do not use for plastic materials).");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::PlasticElastHyper::create_material()
{
  return Teuchos::rcp(new Mat::PlasticElastHyper(this));
}


Mat::PlasticElastHyperType Mat::PlasticElastHyperType::instance_;


Core::Communication::ParObject* Mat::PlasticElastHyperType::Create(const std::vector<char>& data)
{
  Mat::PlasticElastHyper* elhy = new Mat::PlasticElastHyper();
  elhy->Unpack(data);

  return elhy;
}


/*----------------------------------------------------------------------*
 |  initialise static arrays                                 seitz 05/14|
 *----------------------------------------------------------------------*/
Core::LinAlg::Matrix<3, 1> Mat::PlasticElastHyper::prinv_;
Core::LinAlg::Matrix<9, 1> Mat::PlasticElastHyper::CFpiCei_;
Core::LinAlg::Matrix<9, 1> Mat::PlasticElastHyper::CFpi_;
Core::LinAlg::Matrix<9, 1> Mat::PlasticElastHyper::CFpiCe_;
Core::LinAlg::Matrix<6, 1> Mat::PlasticElastHyper::Cpi_;
Core::LinAlg::Matrix<6, 1> Mat::PlasticElastHyper::CpiCCpi_;
Core::LinAlg::Matrix<6, 1> Mat::PlasticElastHyper::ircg_;
Core::LinAlg::Matrix<6, 1> Mat::PlasticElastHyper::Ce_;
Core::LinAlg::Matrix<6, 1> Mat::PlasticElastHyper::Ce2_;
Core::LinAlg::Matrix<6, 1> Mat::PlasticElastHyper::id2V_;
Core::LinAlg::Matrix<6, 1> Mat::PlasticElastHyper::bev_;
Core::LinAlg::Matrix<6, 1> Mat::PlasticElastHyper::be2v_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::invpldefgrd_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::CeM_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::id2_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::CpiC_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::FpiCe_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::FpiTC_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::CeFpiTC_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::be_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::be2_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::Fe_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::beF_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::beFe_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::FCpi_;
Core::LinAlg::Matrix<3, 3> Mat::PlasticElastHyper::beFCpi_;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PlasticElastHyper::PlasticElastHyper()
    : params_(nullptr),
      last_plastic_defgrd_inverse_(Teuchos::null),
      last_alpha_isotropic_(Teuchos::null),
      last_alpha_kinematic_(Teuchos::null),
      activity_state_(Teuchos::null)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PlasticElastHyper::PlasticElastHyper(Mat::PAR::PlasticElastHyper* params)
    : params_(params),
      HepDiss_(Teuchos::null),
      dHepDissdd_(Teuchos::null),
      dHepDissdT_(Teuchos::null),
      dHepDissdTeas_(Teuchos::null)
{
  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;
  for (m = MatParams()->matids_.begin(); m != MatParams()->matids_.end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<Mat::Elastic::Summand> sum = Mat::Elastic::Summand::Factory(matid);
    if (sum == Teuchos::null) FOUR_C_THROW("Failed to allocate");
    potsum_.push_back(sum);
    sum->register_anisotropy_extensions(anisotropy_);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // matid
  int matid = -1;
  if (MatParams() != nullptr) matid = MatParams()->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);
  summandProperties_.Pack(data);

  // plastic anisotropy
  add_to_pack(data, PlAniso_full_);
  add_to_pack(data, InvPlAniso_full_);

  if (MatParams() != nullptr)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->PackSummand(data);
    }
  }

  // plastic history data
  add_to_pack<3, 3>(data, last_plastic_defgrd_inverse_);
  add_to_pack(data, last_alpha_isotropic_);
  add_to_pack<3, 3>(data, last_alpha_kinematic_);

  add_to_pack(data, (int)activity_state_.size());
  for (int i = 0; i < (int)activity_state_.size(); ++i)
    add_to_pack(data, (int)activity_state_.at(i));

  add_to_pack(data, (int)delta_alpha_i_.size());
  for (int i = 0; i < (int)delta_alpha_i_.size(); ++i) add_to_pack(data, delta_alpha_i_.at(i));

  // tsi data
  bool tsi = HepDiss_ != Teuchos::null;
  add_to_pack(data, (int)tsi);
  bool tsi_eas = dHepDissdTeas_ != Teuchos::null;
  add_to_pack(data, (int)tsi_eas);
  if (tsi) add_to_pack(data, (int)dHepDissdd_->at(0).numRows());

  // dissipation mode
  add_to_pack(data, (int)DisMode());

  add_to_pack(data, cpl());
  add_to_pack(data, s());

  anisotropy_.pack_anisotropy(data);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  params_ = nullptr;
  potsum_.clear();

  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover MatParams()
  int matid;
  extract_from_pack(position, data, matid);
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const unsigned int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::PlasticElastHyper*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  summandProperties_.Unpack(position, data);

  // plastic anisotropy
  extract_from_pack(position, data, PlAniso_full_);
  extract_from_pack(position, data, InvPlAniso_full_);

  if (MatParams() != nullptr)  // summands are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = MatParams()->matids_.begin(); m != MatParams()->matids_.end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<Mat::Elastic::Summand> sum = Mat::Elastic::Summand::Factory(matid);
      if (sum == Teuchos::null) FOUR_C_THROW("Failed to allocate");
      potsum_.push_back(sum);
    }

    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->UnpackSummand(data, position);
      potsum_[p]->register_anisotropy_extensions(anisotropy_);
    }
  }

  // plastic history data
  extract_from_pack<3, 3>(position, data, last_plastic_defgrd_inverse_);
  extract_from_pack(position, data, last_alpha_isotropic_);
  extract_from_pack<3, 3>(position, data, last_alpha_kinematic_);

  activity_state_.resize(extract_int(position, data));
  for (int i = 0; i < (int)activity_state_.size(); ++i)
    activity_state_.at(i) = (bool)extract_int(position, data);

  delta_alpha_i_.resize(extract_int(position, data));
  for (int i = 0; i < (int)delta_alpha_i_.size(); ++i)
    delta_alpha_i_.at(i) = extract_double(position, data);

  bool tsi = (bool)extract_int(position, data);
  bool tsi_eas = (bool)extract_int(position, data);
  if (!tsi)
  {
    HepDiss_ = Teuchos::null;
    dHepDissdd_ = Teuchos::null;
    dHepDissdT_ = Teuchos::null;
    dHepDissdTeas_ = Teuchos::null;
  }
  else
  {
    int ngp = last_alpha_isotropic_.size();
    HepDiss_ = Teuchos::rcp(new std::vector<double>(ngp, 0.0));
    int numdofperelement = extract_int(position, data);
    dHepDissdd_ = Teuchos::rcp(new std::vector<Core::LinAlg::SerialDenseVector>(
        ngp, Core::LinAlg::SerialDenseVector(numdofperelement)));
    dHepDissdT_ = Teuchos::rcp(new std::vector<double>(ngp, 0.0));
    if (tsi_eas)
      dHepDissdTeas_ = Teuchos::rcp(new std::vector<Core::LinAlg::SerialDenseVector>(
          ngp, Core::LinAlg::SerialDenseVector(numdofperelement / 3)));
  }

  // dissipation mode
  Inpar::TSI::DissipationMode mode = (Inpar::TSI::DissipationMode)extract_int(position, data);
  SetDissipationMode(mode);

  double cpl = extract_double(position, data);
  double s = extract_double(position, data);
  GetParams(s, cpl);

  anisotropy_.unpack_anisotropy(data, position);

  // in the postprocessing mode, we do not unpack everything we have packed
  // -> position check cannot be done in this case
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::Setup(int numgp, Input::LineDefinition* linedef)
{
  // Read anisotropy
  anisotropy_.set_number_of_gauss_points(numgp);
  anisotropy_.read_anisotropy_from_element(linedef);

  // Setup summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->Setup(numgp, linedef);
  }

  // find out which formulations are used

  summandProperties_.Clear();

  ElastHyperProperties(potsum_, summandProperties_);

  // in this case the mandel stress become non-symmetric and the
  // calculated derivatives have to be extended.
  if (summandProperties_.anisomod == true || summandProperties_.anisoprinc == true)
    FOUR_C_THROW("PlasticElastHyper only for isotropic elastic material!");

  // no visco-elasticity (yet?)
  if (summandProperties_.viscoGeneral)
    FOUR_C_THROW("no visco-elasticity in PlasticElastHyper...yet(?)");

  // check if either zero or three fiber directions are given
  if (linedef->has_named("FIBER1") || linedef->has_named("FIBER2") || linedef->has_named("FIBER3"))
    if (!linedef->has_named("FIBER1") || !linedef->has_named("FIBER2") ||
        !linedef->has_named("FIBER3"))
      FOUR_C_THROW("so3 expects no fibers or 3 fiber directions");

  // plastic anisotropy
  SetupHillPlasticity(linedef);

  // setup plastic history variables
  Core::LinAlg::Matrix<3, 3> tmp(true);
  last_alpha_isotropic_.resize(numgp, 0.);
  last_alpha_kinematic_.resize(numgp, tmp);
  for (int i = 0; i < 3; i++) tmp(i, i) = 1.;
  last_plastic_defgrd_inverse_.resize(numgp, tmp);
  activity_state_.resize(numgp, false);
  delta_alpha_i_.resize(numgp, 0.);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::SetupTSI(const int numgp, const int numdofperelement, const bool eas,
    const Inpar::TSI::DissipationMode mode)
{
  // dissipation mode
  if (mode == Inpar::TSI::pl_multiplier)
    if (MatParams()->rY_11_ != 0. || MatParams()->rY_22_ != 0. || MatParams()->rY_33_ != 0. ||
        MatParams()->rY_12_ != 0. || MatParams()->rY_23_ != 0. || MatParams()->rY_13_ != 0.)
      FOUR_C_THROW("TSI with Hill plasticity not available with DISSIPATION_MODE pl_multiplier");
  SetDissipationMode(mode);

  // allocate memory
  HepDiss_ = Teuchos::rcp(new std::vector<double>(numgp, 0.0));
  dHepDissdd_ = Teuchos::rcp(new std::vector<Core::LinAlg::SerialDenseVector>(
      numgp, Core::LinAlg::SerialDenseVector(numdofperelement)));
  dHepDissdT_ = Teuchos::rcp(new std::vector<double>(numgp, 0.0));
  if (eas)
    dHepDissdTeas_ = Teuchos::rcp(new std::vector<Core::LinAlg::SerialDenseVector>(
        numgp, Core::LinAlg::SerialDenseVector(numdofperelement / 3)));

  // no TSI with kinematic hardening yet
  // be aware that in that case, another dependency of the NCP function on the
  // temperature arises, namely via H^k=H^k(T) in the computation of the effective
  // stress eta. Without this term, the only dependency of the NCP function is
  // via the effective yield stress Y^pl=Y^pl(T)
  if (Kinhard() != 0.) FOUR_C_THROW("no kinematic hardening for TSI (yet)");
  // no TSI with plastic spin yet
  // besides the kinematic hardening term (see above) there is not much to do
  // just add the derivatives of theating and NCP function in the
  // evaluate_nc_pand_spin(...) function
  if (PlSpinChi() != 0.) FOUR_C_THROW("no thermo-plasticitiy with plastic spin");

  /// Hill TSI only with pl_flow dissipation
  if (MatParams()->rY_11_ != 0. && mode == Inpar::TSI::pl_multiplier)
    FOUR_C_THROW("hill thermo plasticity not with dissipation mode pl_multiplier");

  /// viscoplastic TSI only with  pl_flow dissipation
  if (Visc() != 0. && mode == Inpar::TSI::pl_multiplier)
    FOUR_C_THROW("thermo-visco-plasticity not with dissipation mode pl_multiplier");
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::SetupHillPlasticity(Input::LineDefinition* linedef)
{
  // check if parameters are valid
  if (MatParams()->rY_11_ != 0. || MatParams()->rY_22_ != 0. || MatParams()->rY_33_ != 0. ||
      MatParams()->rY_12_ != 0. || MatParams()->rY_23_ != 0. || MatParams()->rY_13_ != 0.)
    if (MatParams()->rY_11_ <= 0. || MatParams()->rY_22_ <= 0. || MatParams()->rY_33_ <= 0. ||
        MatParams()->rY_12_ <= 0. || MatParams()->rY_23_ <= 0. || MatParams()->rY_13_ <= 0.)
      FOUR_C_THROW("Hill parameters all must be positive (incomplete set?)");

  // all (optional) Hill parameters are zero (default value)
  // --> we want to do von Mises plasticity
  if (MatParams()->rY_11_ == 0. && MatParams()->rY_22_ == 0. && MatParams()->rY_33_ == 0. &&
      MatParams()->rY_12_ == 0. && MatParams()->rY_23_ == 0. && MatParams()->rY_13_ == 0.)
  {
    PlAniso_full_.Clear();
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        if (i == j)
          PlAniso_full_(i, j) = 2. / 3.;
        else
          PlAniso_full_(i, j) = -1. / 3.;
    for (int i = 3; i < 6; i++) PlAniso_full_(i, i) = 1.;
    InvPlAniso_full_.Update(PlAniso_full_);
  }
  // we do Hill plasticity
  else
  {
    // anisotropy directions
    std::vector<Core::LinAlg::Matrix<3, 1>> directions(3);

    // first anisotropy direction
    if (linedef->has_named("FIBER1"))
    {
      std::vector<double> fiber;
      linedef->extract_double_vector("FIBER1", fiber);
      double fnorm = 0.;
      // normalization
      for (int i = 0; i < 3; ++i) fnorm += fiber[i] * fiber[i];
      fnorm = sqrt(fnorm);
      if (fnorm == 0.) FOUR_C_THROW("Fiber vector has norm zero");

      // fill final fiber vector
      for (int i = 0; i < 3; ++i) directions.at(0)(i) = fiber[i] / fnorm;
    }
    // second anisotropy direction
    if (linedef->has_named("FIBER2"))
    {
      std::vector<double> fiber;
      linedef->extract_double_vector("FIBER2", fiber);
      double fnorm = 0.;
      // normalization
      for (int i = 0; i < 3; ++i) fnorm += fiber[i] * fiber[i];
      fnorm = sqrt(fnorm);
      if (fnorm == 0.) FOUR_C_THROW("Fiber vector has norm zero");

      // fill final fiber vector
      for (int i = 0; i < 3; ++i) directions.at(1)(i) = fiber[i] / fnorm;
    }
    // third anisotropy direction
    if (linedef->has_named("FIBER3"))
    {
      std::vector<double> fiber;
      linedef->extract_double_vector("FIBER3", fiber);
      double fnorm = 0.;
      // normalization
      for (int i = 0; i < 3; ++i) fnorm += fiber[i] * fiber[i];
      fnorm = sqrt(fnorm);
      if (fnorm == 0.) FOUR_C_THROW("Fiber vector has norm zero");

      // fill final fiber vector
      for (int i = 0; i < 3; ++i) directions.at(2)(i) = fiber[i] / fnorm;
    }

    // check orthogonality
    Core::LinAlg::Matrix<1, 1> matrix1;
    matrix1.MultiplyTN(directions.at(0), directions.at(1));
    if (std::abs(matrix1(0, 0)) > 1.e-16) FOUR_C_THROW("fiber directions not orthogonal");
    matrix1.MultiplyTN(directions.at(0), directions.at(2));
    if (std::abs(matrix1(0, 0)) > 1.e-16) FOUR_C_THROW("fiber directions not orthogonal");
    matrix1.MultiplyTN(directions.at(2), directions.at(1));
    if (std::abs(matrix1(0, 0)) > 1.e-16) FOUR_C_THROW("fiber directions not orthogonal");

    // check right-handed trihedron
    Core::LinAlg::Matrix<3, 1> A0xA1;
    A0xA1(0) =
        (directions.at(0)(1) * directions.at(1)(2) - directions.at(0)(2) * directions.at(1)(1));
    A0xA1(1) =
        (directions.at(0)(2) * directions.at(1)(0) - directions.at(0)(0) * directions.at(1)(2));
    A0xA1(2) =
        (directions.at(0)(0) * directions.at(1)(1) - directions.at(0)(1) * directions.at(1)(0));
    A0xA1.Update(-1., directions.at(2), 1.);
    if (A0xA1.Norm2() > 1.e-8) FOUR_C_THROW("fibers don't form right-handed trihedron");

    // setup structural tensor for first and second direction
    // (as the directions are orthogonal, 2 structural tensors are sufficient)
    Core::LinAlg::Matrix<3, 3> M0;
    M0.MultiplyNT(directions.at(0), directions.at(0));
    Core::LinAlg::Matrix<3, 3> M1;
    M1.MultiplyNT(directions.at(1), directions.at(1));
    Core::LinAlg::Matrix<3, 3> M2;
    M2.MultiplyNT(directions.at(2), directions.at(2));

    double alpha1 = 2. / 3. / MatParams()->rY_11_ / MatParams()->rY_11_;
    double alpha2 = 2. / 3. / MatParams()->rY_22_ / MatParams()->rY_22_;
    double alpha3 = 2. / 3. / MatParams()->rY_33_ / MatParams()->rY_33_;
    double alpha4 = 1. / 3. / MatParams()->rY_12_ / MatParams()->rY_12_;
    double alpha5 = 1. / 3. / MatParams()->rY_23_ / MatParams()->rY_23_;
    double alpha6 = 1. / 3. / MatParams()->rY_13_ / MatParams()->rY_13_;

    // calculate plastic anisotropy tensor
    PlAniso_full_.Clear();
    ElastSymTensorMultiply(PlAniso_full_, alpha1, M0, M0, 1.);
    ElastSymTensorMultiply(PlAniso_full_, alpha2, M1, M1, 1.);
    ElastSymTensorMultiply(PlAniso_full_, alpha3, M2, M2, 1.);
    ElastSymTensorMultiplyAddSym(PlAniso_full_, 0.5 * (alpha3 - alpha1 - alpha2), M0, M1, 1.);
    ElastSymTensorMultiplyAddSym(PlAniso_full_, 0.5 * (alpha1 - alpha2 - alpha3), M1, M2, 1.);
    ElastSymTensorMultiplyAddSym(PlAniso_full_, 0.5 * (alpha2 - alpha3 - alpha1), M0, M2, 1.);
    AddSymmetricHolzapfelProduct(PlAniso_full_, M0, M1, alpha4);
    AddSymmetricHolzapfelProduct(PlAniso_full_, M1, M2, alpha5);
    AddSymmetricHolzapfelProduct(PlAniso_full_, M0, M2, alpha6);

    // we need this matrix to get rid of the zero eigenvalue to be able to invert
    // the anisotropy tensor. After the inversion we expand the tensor again to 6x6
    // so that we have the correct pseudo-inverse.
    Core::LinAlg::Matrix<6, 5> red(true);
    red(0, 0) = 1.;
    red(1, 1) = 1.;
    red(2, 0) = -1.;
    red(2, 1) = -1.;
    red(3, 2) = 1.;
    red(4, 3) = 1.;
    red(5, 4) = 1.;

    // Invert plastic anisotropy tensor
    Core::LinAlg::Matrix<6, 5> tmp;
    tmp.Multiply(PlAniso_full_, red);
    Core::LinAlg::Matrix<5, 5> tmp55;
    tmp55.MultiplyTN(red, tmp);
    Core::LinAlg::FixedSizeSerialDenseSolver<5, 5, 1> solver;
    solver.SetMatrix(tmp55);
    int err2 = solver.Factor();
    int err = solver.Invert();
    if ((err != 0) || (err2 != 0)) FOUR_C_THROW("Inversion of plastic anisotropy tensor failed");
    tmp.MultiplyNT(red, tmp55);
    InvPlAniso_full_.MultiplyNT(tmp, red);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate elastic stress and stiffness                   seitz 05/14 |
 *----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::EvaluateElast(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<3, 3>* deltaLp, Core::LinAlg::Matrix<6, 1>* pk2,
    Core::LinAlg::Matrix<6, 6>* cmat, const int gp, const int eleGID)
{
  Core::LinAlg::Matrix<3, 1> dPI(true);
  Core::LinAlg::Matrix<6, 1> ddPII(true);

  evaluate_kin_quant_elast(defgrd, deltaLp, gp);
  ElastHyperEvaluateInvariantDerivatives(
      prinv_, dPI, ddPII, potsum_, summandProperties_, gp, eleGID);

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  pk2->Clear();
  cmat->Clear();

  // isotropic elasticity in coupled strain energy format
  // isotropic elasticity in decoupled ("mod") format go here as well
  // as the modified gammas and deltas have been converted
  if (summandProperties_.isoprinc || summandProperties_.isomod)
    evaluate_isotropic_princ_elast(*pk2, *cmat, dPI, ddPII);
  else
    FOUR_C_THROW("only isotropic hyperelastic materials");

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate elastic strain energy                          seitz 03/16 |
 *----------------------------------------------------------------------*/
double Mat::PlasticElastHyper::StrainEnergyTSI(
    const Core::LinAlg::Matrix<3, 3>& defgrd, const int gp, const int eleGID, const double temp)
{
  double psi = 0.;

  Core::LinAlg::Matrix<3, 3> Fe;
  Fe.Multiply(defgrd, last_plastic_defgrd_inverse_[gp]);
  Core::LinAlg::Matrix<3, 3> elRCG;
  elRCG.MultiplyTN(Fe, Fe);
  Core::LinAlg::Matrix<6, 1> elRCGv;
  for (int i = 0; i < 3; ++i) elRCGv(i) = elRCG(i, i);
  elRCGv(3) = elRCG(0, 1) + elRCG(1, 0);
  elRCGv(4) = elRCG(2, 1) + elRCG(1, 2);
  elRCGv(5) = elRCG(0, 2) + elRCG(2, 0);
  Core::LinAlg::Matrix<3, 1> prinv;
  Core::LinAlg::Voigt::Strains::invariants_principal(prinv, elRCGv);
  Core::LinAlg::Matrix<3, 1> modinv;
  invariants_modified(modinv, prinv);

  // loop map of associated potential summands
  Core::LinAlg::Matrix<6, 1> glstrain(true);
  Core::LinAlg::Matrix<6, 1> idv(true);
  for (int i = 0; i < 3; ++i) idv(i) = 1.0;
  glstrain.Update(0.5, elRCGv, 0.0);
  glstrain.Update(-0.5, idv, 1.0);
  for (unsigned int p = 0; p < potsum_.size(); ++p)
    potsum_[p]->AddStrainEnergy(psi, prinv, modinv, glstrain, gp, eleGID);

  double dPj1 = 0.;
  for (unsigned int p = 0; p < potsum_.size(); ++p)
    potsum_[p]->AddCoupDerivVol(sqrt(prinv(0)), &dPj1, nullptr, nullptr, nullptr);
  psi -= 3. * Cte() * (temp - InitTemp()) * dPj1;



  return psi;
}

/*----------------------------------------------------------------------*
 |  evaluate thermal stress and stiffness                   seitz 06/14 |
 *----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::evaluate_thermal_stress(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const double temp, Core::LinAlg::Matrix<6, 1>* pk2, Core::LinAlg::Matrix<6, 6>* cmat,
    const int gp, const int eleGID)
{
  // do TSI only for decoupled isotropic materials. By doing so, the stresses
  // due to thermal expansion can be easily calculated by the volumetric
  // part of the strain energy function
  if (summandProperties_.anisomod || summandProperties_.anisoprinc)
    FOUR_C_THROW(
        "TSI with semi-Smooth Newton type plasticity algorithm only "
        "with isotropic strain energy functions");

  // temperature difference
  double deltaT = temp - InitTemp();

  // we are only interested in the volumetric response
  // which is for decoupled strain energy functions defined by
  // modinv_3 = J only.
  Core::LinAlg::Matrix<3, 1> modinv(true);
  modinv(2) = defgrd->Determinant();
  Core::LinAlg::Matrix<3, 1> dPmodI;
  Core::LinAlg::Matrix<6, 1> ddPmodII;
  double dddPmodIII = 0.;

  // loop map of associated potential summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->add_derivatives_modified(dPmodI, ddPmodII, modinv, gp, eleGID);
    potsum_[p]->Add3rdVolDeriv(modinv, dddPmodIII);
    potsum_[p]->AddCoupDerivVol(modinv(2), &dPmodI(2), &ddPmodII(2), &dddPmodIII, nullptr);
  }

  // inverse RCG
  Core::LinAlg::Matrix<3, 3> invRCG;
  invRCG.MultiplyTN(*defgrd, *defgrd);
  invRCG.Invert();
  Core::LinAlg::Matrix<6, 1> icg;
  for (int i = 0; i < 3; ++i) icg(i) = invRCG(i, i);
  icg(3) = invRCG(0, 1);
  icg(4) = invRCG(1, 2);
  icg(5) = invRCG(0, 2);

  pk2->Update(-3. * Cte() * deltaT * modinv(2) * ddPmodII(2), icg, 1.);
  ElastSymTensorMultiply(
      *cmat, -3. * Cte() * deltaT * modinv(2) * modinv(2) * dddPmodIII, invRCG, invRCG, 1.);
  ElastSymTensorMultiply(*cmat, -3. * Cte() * deltaT * modinv(2) * ddPmodII(2), invRCG, invRCG, 1.);
  ElastSymTensor_o_Multiply(
      *cmat, +6. * Cte() * deltaT * modinv(2) * ddPmodII(2), invRCG, invRCG, 1.);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate thermal stress and stiffness                   seitz 06/14 |
 *----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::EvaluateCTvol(const Core::LinAlg::Matrix<3, 3>* defgrd,
    Core::LinAlg::Matrix<6, 1>* cTvol, Core::LinAlg::Matrix<6, 6>* dCTvoldE, const int gp,
    const int eleGID)
{
  // do TSI only for decoupled isotropic materials. By doing so, the stresses
  // due to thermal expansion can be easily calculated by the volumetric
  // part of the strain energy function
  if (summandProperties_.anisomod || summandProperties_.anisoprinc)
    FOUR_C_THROW(
        "TSI with semi-Smooth Newton type plasticity algorithm only "
        "with isotropic strain energy functions");

  // we are only interested in the volumetric response
  // which is for decoupled strain energy functions defined by
  // modinv_3 = J only.
  Core::LinAlg::Matrix<3, 1> modinv(true);
  modinv(2) = defgrd->Determinant();
  Core::LinAlg::Matrix<3, 1> dPmodI;
  Core::LinAlg::Matrix<6, 1> ddPmodII;
  double dddPmodIII = 0.;

  // loop map of associated potential summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->add_derivatives_modified(dPmodI, ddPmodII, modinv, gp, eleGID);
    potsum_[p]->Add3rdVolDeriv(modinv, dddPmodIII);
    potsum_[p]->AddCoupDerivVol(modinv(2), &dPmodI(2), &ddPmodII(2), &dddPmodIII, nullptr);
  }

  // clear
  cTvol->Clear();
  dCTvoldE->Clear();

  // inverse RCG
  Core::LinAlg::Matrix<3, 3> invRCG;
  invRCG.MultiplyTN(*defgrd, *defgrd);
  invRCG.Invert();
  Core::LinAlg::Matrix<6, 1> icg;
  for (int i = 0; i < 3; ++i) icg(i) = invRCG(i, i);
  icg(3) = invRCG(0, 1);
  icg(4) = invRCG(1, 2);
  icg(5) = invRCG(0, 2);

  cTvol->Update(-3. * Cte() * modinv(2) * ddPmodII(2), icg, 1.);
  ElastSymTensorMultiply(
      *dCTvoldE, -3. * Cte() * modinv(2) * modinv(2) * dddPmodIII, invRCG, invRCG, 1.);
  ElastSymTensorMultiply(*dCTvoldE, -3. * Cte() * modinv(2) * ddPmodII(2), invRCG, invRCG, 1.);
  ElastSymTensor_o_Multiply(*dCTvoldE, +6. * Cte() * modinv(2) * ddPmodII(2), invRCG, invRCG, 1.);

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate Gough Joule Effect                             seitz 10/15 |
 *----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::EvaluateGoughJoule(
    const double j, const int gp, const int eleGID, double& he_fac, double& he_fac_deriv)
{
  // we are only interested in the volumetric response
  // which is for decoupled strain energy functions defined by
  // modinv_3 = J only.
  Core::LinAlg::Matrix<3, 1> modinv(true);
  modinv(2) = j;
  Core::LinAlg::Matrix<3, 1> dPmodI;
  Core::LinAlg::Matrix<6, 1> ddPmodII;
  double dddPmodIII = 0.;

  // loop map of associated potential summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->add_derivatives_modified(dPmodI, ddPmodII, modinv, gp, eleGID);
    potsum_[p]->Add3rdVolDeriv(modinv, dddPmodIII);
    potsum_[p]->AddCoupDerivVol(modinv(2), &dPmodI(2), &ddPmodII(2), &dddPmodIII, nullptr);
  }

  he_fac = -3. * Cte() * ddPmodII(2);
  he_fac_deriv = -3. * Cte() * dddPmodIII;
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate plastic stress and stiffness                   seitz 05/14 |
 *----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::EvaluatePlast(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<3, 3>* deltaDp, const double* temp, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 6>* dPK2dDp, Core::LinAlg::Matrix<6, 1>* NCP,
    Core::LinAlg::Matrix<6, 6>* dNCPdC, Core::LinAlg::Matrix<6, 6>* dNCPdDp, bool* active,
    bool* elast, bool* as_converged, const int gp, Core::LinAlg::Matrix<6, 1>* dNCPdT,
    Core::LinAlg::Matrix<6, 1>* dHdC, Core::LinAlg::Matrix<6, 1>* dHdDp, const double dt,
    const int eleGID, Core::LinAlg::Matrix<6, 1>* cauchy, Core::LinAlg::Matrix<6, 6>* d_cauchy_ddp,
    Core::LinAlg::Matrix<6, 6>* d_cauchy_dC, Core::LinAlg::Matrix<6, 9>* d_cauchy_dF,
    Core::LinAlg::Matrix<6, 1>* d_cauchy_dT)
{
  int check = +(cauchy != nullptr) + (d_cauchy_dC != nullptr) + (d_cauchy_dF != nullptr) +
              (d_cauchy_ddp != nullptr);
  if (!(check == 0 || check == 4)) FOUR_C_THROW("some inconsistency with provided variables");

  Core::LinAlg::Matrix<3, 1> dPI;
  Core::LinAlg::Matrix<6, 1> ddPII;
  Core::LinAlg::Matrix<3, 1> gamma(true);
  Core::LinAlg::Matrix<8, 1> delta(true);

  if (evaluate_kin_quant_plast(defgrd, deltaDp, gp, params)) return;
  ElastHyperEvaluateInvariantDerivatives(
      prinv_, dPI, ddPII, potsum_, summandProperties_, gp, eleGID);
  if (temp && cauchy) add_thermal_expansion_derivs(prinv_, dPI, ddPII, gp, eleGID, *temp);
  CalculateGammaDelta(gamma, delta, prinv_, dPI, ddPII);

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  dPK2dDp->Clear();
  NCP->Clear();
  dNCPdC->Clear();
  dNCPdDp->Clear();
  if (dNCPdT != nullptr) dNCPdT->Clear();

  // new temporary matrices
  Core::LinAlg::Matrix<3, 3> mStr;  // Mandel stress tensor
  Core::LinAlg::Matrix<6, 6> dMdC;  // derivative of Mandel stress w.r.t. RCG
  Core::LinAlg::Matrix<6, 9>
      dMdFpinv;  // derivative of Mandel stress w.r.t. inverse plastic deformation gradient
  Core::LinAlg::Matrix<6, 9> dPK2dFpinv;
  Core::LinAlg::Matrix<6, 9> d_cauchy_dFpi;

  // isotropic elasticity in coupled strain energy format
  // isotropic elasticity in decoupled ("mod") format go here as well
  // as the modified gammas and deltas have been converted
  if (summandProperties_.isoprinc || summandProperties_.isomod)
    evaluate_isotropic_princ_plast(dPK2dFpinv, mStr, dMdC, dMdFpinv, gamma, delta);
  else
    FOUR_C_THROW("only isotropic hypereleastic materials");

  if (cauchy)
    EvaluateCauchyPlast(
        dPI, ddPII, defgrd, *cauchy, d_cauchy_dFpi, *d_cauchy_dC, *d_cauchy_dF, d_cauchy_dT);

  evaluate_ncp(&mStr, &dMdC, &dMdFpinv, &dPK2dFpinv, deltaDp, gp, temp, NCP, dNCPdC, dNCPdDp,
      dNCPdT, dPK2dDp, active, elast, as_converged, dHdC, dHdDp, params, dt, &d_cauchy_dFpi,
      d_cauchy_ddp);

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate NCP function                                   seitz 05/14 |
 *----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::evaluate_ncp(const Core::LinAlg::Matrix<3, 3>* mStr,
    const Core::LinAlg::Matrix<6, 6>* dMdC, const Core::LinAlg::Matrix<6, 9>* dMdFpinv,
    const Core::LinAlg::Matrix<6, 9>* dPK2dFpinv, const Core::LinAlg::Matrix<3, 3>* deltaDp,
    const int gp, const double* temp, Core::LinAlg::Matrix<6, 1>* NCP,
    Core::LinAlg::Matrix<6, 6>* dNCPdC, Core::LinAlg::Matrix<6, 6>* dNCPdDp,
    Core::LinAlg::Matrix<6, 1>* dNCPdT, Core::LinAlg::Matrix<6, 6>* dPK2dDp, bool* active,
    bool* elast, bool* as_converged, Core::LinAlg::Matrix<6, 1>* dHdC,
    Core::LinAlg::Matrix<6, 1>* dHdDp, Teuchos::ParameterList& params, const double dt,
    const Core::LinAlg::Matrix<6, 9>* d_cauchy_dFpi, Core::LinAlg::Matrix<6, 6>* d_cauchy_ddp)
{
  const double sq = sqrt(2. / 3.);
  Core::LinAlg::Matrix<6, 1> tmp61;
  double dT = 0.;
  if (dNCPdT)
    dT = *temp - InitTemp();
  else
    dT = 0.;

  // deviatoric projection tensor
  Core::LinAlg::Matrix<6, 6> pdev(true);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (i == j)
        pdev(i, j) = 2. / 3.;
      else
        pdev(i, j) = -1. / 3.;
  for (int i = 3; i < 6; i++) pdev(i, i) = 1.;

  // effective stress
  Core::LinAlg::Matrix<3, 3> eta(*mStr);
  for (int i = 0; i < 3; i++)
    eta(i, i) -= 1. / 3. * ((*mStr)(0, 0) + (*mStr)(1, 1) + (*mStr)(2, 2));
  eta.Update(2. / 3. * Kinhard(), last_alpha_kinematic_[gp], 1.);
  eta.Update(-2. / 3. * Kinhard(), *deltaDp, 1.);

  // in stress-like voigt notation
  Core::LinAlg::Matrix<6, 1> eta_v;      // in stress-like voigt notation
  Core::LinAlg::Matrix<6, 1> etatr_v;    // in stress-like voigt notation
  Core::LinAlg::Matrix<6, 1> eta_s_v;    // in stress-like voigt notation
  Core::LinAlg::Matrix<6, 1> deltaDp_v;  // in stress-like voigt notation
  for (int i = 0; i < 3; i++)
  {
    eta_v(i) = eta(i, i);
    deltaDp_v(i) = (*deltaDp)(i, i);
  }
  eta_v(3) = .5 * (eta(1, 0) + eta(0, 1));
  deltaDp_v(3) = .5 * ((*deltaDp)(0, 1) + (*deltaDp)(1, 0));
  eta_v(4) = .5 * (eta(1, 2) + eta(2, 1));
  deltaDp_v(4) = .5 * ((*deltaDp)(2, 1) + (*deltaDp)(1, 2));
  eta_v(5) = .5 * (eta(2, 0) + eta(0, 2));
  deltaDp_v(5) = .5 * ((*deltaDp)(0, 2) + (*deltaDp)(2, 0));

  // trial effective stress
  etatr_v.Update(eta_v);
  etatr_v.Multiply(cpl(), InvPlAniso_full_, deltaDp_v, 1.);

  // in strain-like voigt notation
  Core::LinAlg::Matrix<6, 1> eta_v_strainlike(eta_v);          // in strain-like voigt notation
  Core::LinAlg::Matrix<6, 1> etatr_v_strainlike(etatr_v);      // in strain-like voigt notation
  Core::LinAlg::Matrix<6, 1> deltaDp_v_strainlike(deltaDp_v);  // in strain-like voigt notation
  for (int i = 3; i < 6; i++)
  {
    eta_v_strainlike(i) *= 2.;
    etatr_v_strainlike(i) *= 2.;
    deltaDp_v_strainlike(i) *= 2.;
  }

  // different tensor norms
  tmp61.Multiply(PlAniso_full_, eta_v);
  double absHeta = norm_stress_like(tmp61);
  double abseta_H = tmp61.Dot(eta_v_strainlike);
  if (abseta_H < -1.e-16)
    FOUR_C_THROW("this should not happen. eta : H : eta =%f < 0", abseta_H);
  else if (abseta_H >= 0.)
    abseta_H = sqrt(abseta_H);
  else
    abseta_H = 0.;
  double dDpHeta = tmp61.Dot(deltaDp_v_strainlike);
  tmp61.Multiply(PlAniso_full_, etatr_v);
  double absetatr_H = tmp61.Dot(etatr_v_strainlike);
  if (absetatr_H < -1.e-16)
    FOUR_C_THROW("this should not happen. eta_tr : H : eta_tr =%f < 0", absetatr_H);
  else if (absetatr_H >= 0.)
    absetatr_H = sqrt(absetatr_H);
  else
    absetatr_H = 0.;
  Core::LinAlg::Matrix<6, 1> HdDp;
  HdDp.Multiply(PlAniso_full_, deltaDp_v);
  Core::LinAlg::Matrix<6, 1> HdDp_strainlike;
  HdDp_strainlike.Multiply(PlAniso_full_, deltaDp_v_strainlike);
  Core::LinAlg::Matrix<6, 1> HetaH_strainlike;
  tmp61.Multiply(PlAniso_full_, eta_v_strainlike);
  HetaH_strainlike.Multiply(PlAniso_full_, tmp61);

  // isotropic hardening increment
  delta_alpha_i_[gp] = 0.;
  if (dDpHeta > 0. && absHeta > 0.)
    delta_alpha_i_[gp] = sq * dDpHeta * abseta_H / (absHeta * absHeta);
  // new isotropic hardening value
  const double aI = last_alpha_isotropic_[gp] + delta_alpha_i_[gp];

  // current yield stress equivalent (yield stress scaled by sqrt(2/3))
  double ypl =
      sq *
      ((Infyield() * (1. - HardSoft() * dT) - Inityield() * (1. - YieldSoft() * dT)) *
              (1. - exp(-Expisohard() * aI)) +
          Isohard() * (1. - HardSoft() * dT) * aI + Inityield() * (1. - YieldSoft() * dT)) *
      pow(1. + Visc() * (1. - ViscSoft() * dT) * delta_alpha_i_[gp] / dt, ViscRate());

  double dYpldT = sq *
                  ((Infyield() * (-HardSoft()) - Inityield() * (-YieldSoft())) *
                          (1. - exp(-Expisohard() * aI)) -
                      Isohard() * HardSoft() * aI - Inityield() * YieldSoft()) *
                  pow(1. + Visc() * (1. - ViscSoft() * dT) * delta_alpha_i_[gp] / dt, ViscRate());

  dYpldT += sq *
            ((Infyield() * (1. - HardSoft() * dT) - Inityield() * (1. - YieldSoft() * dT)) *
                    (1. - exp(-Expisohard() * aI)) +
                Isohard() * (1. - HardSoft() * dT) * aI + Inityield() * (1. - YieldSoft() * dT)) *
            pow(1. + Visc() * (1. - ViscSoft() * dT) * delta_alpha_i_[gp] / dt, ViscRate() - 1.) *
            ViscRate() * delta_alpha_i_[gp] / dt * Visc() * (-ViscSoft());

  // Factor of derivative of Y^pl w.r.t. delta alpha ^i
  // we have added the factor sqrt(2/3) from delta_alpha_i=sq*... here
  double dYplDai =
      2. / 3. *
      (+Isohard() * (1. - HardSoft() * dT) +
          (Infyield() * (1. - HardSoft() * dT) - Inityield() * (1. - YieldSoft() * dT)) *
              Expisohard() * exp(-Expisohard() * aI)) *
      pow(1. + Visc() * (1. - ViscSoft() * dT) * delta_alpha_i_[gp] / dt, ViscRate());
  dYplDai += 2. / 3. *
             ((Infyield() * (1. - HardSoft() * dT) - Inityield() * (1. - YieldSoft() * dT)) *
                     (1. - exp(-Expisohard() * aI)) +
                 Isohard() * (1. - HardSoft() * dT) * aI + Inityield() * (1. - YieldSoft() * dT)) *
             pow(1. + Visc() * (1. - ViscSoft() * dT) * delta_alpha_i_[gp] / dt, ViscRate() - 1) *
             ViscRate() * Visc() * (1. - ViscSoft() * dT) / dt;

  // safety check: due to thermal softening, we might get a negative yield stress
  // in that case, no viable solution to the plasticity probleme exists, since we
  // need abs(dev(Sigma)) - Y = 0 --> abs(dev(Sigma))=Y, but abs(dev(Sigma))>0, Y<0 --> error
  if (ypl <= 0.)
    FOUR_C_THROW(
        "negative effective yield stress! Thermal softening and large temperature increase???");

  // activity state check
  if (ypl < absetatr_H)
  {
    if (activity_state_[gp] == false)  // gp switches state
    {
      if (abs(ypl - absetatr_H) > AS_CONVERGENCE_TOL * Inityield() ||
          deltaDp->NormInf() > AS_CONVERGENCE_TOL * Inityield() / cpl())
        *as_converged = false;
    }
    activity_state_[gp] = true;
    *active = true;
  }
  else
  {
    if (activity_state_[gp] == true)  // gp switches state
    {
      if (abs(ypl - absetatr_H) > AS_CONVERGENCE_TOL * Inityield() ||
          deltaDp->NormInf() > AS_CONVERGENCE_TOL * Inityield() / cpl())
        *as_converged = false;
    }
    activity_state_[gp] = false;
    *active = false;
  }

  // check, if gp is at the corner of the root of the NCP function
  if (abs(ypl - absetatr_H) > AS_CONVERGENCE_TOL * Inityield() ||
      deltaDp->NormInf() > AS_CONVERGENCE_TOL * Inityield() / cpl())
  { /* gp not at the corner point */
  }
  else
  {
    /* gp at the corner point --> elastic and plastic branch are equally valid, take the elastic
     * one*/
    *active = false;
    activity_state_[gp] = false;
  }

  // these cases have some terms in common
  if (*active || dDpHeta > 0. || deltaDp->NormInf() > 0.)
  {
    // damping parameter apl
    double apl = 1.;
    if (ypl < abseta_H) apl = ypl / abseta_H;

    // eta_s to abbreviate calculation of the derivatives
    if (absetatr_H > 0.) eta_s_v.Update((1. - s()) * ypl / absetatr_H, etatr_v, 1.);
    eta_s_v.Update(apl * s(), eta_v, 1.);

    // matrix exponential derivative
    Core::LinAlg::Matrix<6, 6> Dexp(false);
    Core::LinAlg::Matrix<3, 3> tmp(*deltaDp);
    tmp.Scale(-1.);
    matrix_exponential_derivative_sym3x3(tmp, Dexp);

    // Derivative of inverse plastic deformation gradient
    Core::LinAlg::Matrix<9, 6> dFpiDdeltaDp(true);
    for (int A = 0; A < 3; A++)
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
          for (int i = 0; i < 6; i++)
            if (i <= 2)
              dFpiDdeltaDp(
                  Core::LinAlg::Voigt::IndexMappings::non_symmetric_tensor_to_voigt9_index(A, a),
                  i) -=
                  last_plastic_defgrd_inverse_[gp](A, b) *
                  Dexp(Core::LinAlg::Voigt::IndexMappings::symmetric_tensor_to_voigt6_index(b, a),
                      i);
            else
              dFpiDdeltaDp(
                  Core::LinAlg::Voigt::IndexMappings::non_symmetric_tensor_to_voigt9_index(A, a),
                  i) -=
                  2. * last_plastic_defgrd_inverse_[gp](A, b) *
                  Dexp(Core::LinAlg::Voigt::IndexMappings::symmetric_tensor_to_voigt6_index(b, a),
                      i);

    // derivative of mandel stress
    // we spare the deviatoric projection of the mandel stress derivative to get the effective
    // stress derivative. It is enforced implicitly as the "detaddp" is always contracted with
    // deviatoric tensors.
    Core::LinAlg::Matrix<6, 6> detaddp;
    detaddp.Multiply(*dMdFpinv, dFpiDdeltaDp);
    detaddp.Update(-2. / 3. * Kinhard(), pdev, 1.);
    dPK2dDp->Multiply(*dPK2dFpinv, dFpiDdeltaDp);
    if (d_cauchy_ddp) d_cauchy_ddp->Multiply(*d_cauchy_dFpi, dFpiDdeltaDp);

    // TSI
    if (dNCPdT != nullptr)
    {
      if (DisMode() == Inpar::TSI::Taylor_Quinney)
      {
        double plHeating = TaylorQuinney() * eta_v_strainlike.Dot(deltaDp_v);
        Core::LinAlg::Matrix<6, 1> dHpDeta(true);
        dHpDeta.Update(TaylorQuinney(), deltaDp_v_strainlike, 1.);
        dHdC->MultiplyTN(*dMdC, dHpDeta);
        dHdDp->MultiplyTN(detaddp, dHpDeta);
        dHdDp->Update(TaylorQuinney(), eta_v_strainlike, 1.);
        plHeating /= dt;
        dHdC->Scale(1. / dt);
        dHdDp->Scale(1. / dt);
        HepDiss(gp) += plHeating;
      }
      else
      {
        // plastic heating
        double plHeating = (0. - Isohard() * HardSoft() * aI -
                               (Infyield() * HardSoft() - Inityield() * YieldSoft()) *
                                   (1. - exp(-Expisohard() * aI))) *
                           (*temp) * delta_alpha_i_[gp];
        switch (DisMode())
        {
          case Inpar::TSI::pl_multiplier:
            plHeating +=
                delta_alpha_i_[gp] *
                (0. + Inityield() * (1. - YieldSoft() * dT) +
                    Isohard() * (1. - HardSoft() * dT) * aI +
                    (Infyield() * (1. - HardSoft() * dT) - Inityield() * (1. - YieldSoft() * dT)) *
                        (1. - exp(-Expisohard() * aI)));
            break;
          case Inpar::TSI::pl_flow:
            plHeating += eta_v_strainlike.Dot(deltaDp_v);
            break;
          default:
            FOUR_C_THROW("unknown plastic dissipation mode: %d", DisMode());
            break;
        }

        // derivative w.r.t. temperature
        double dPlHeatingDT = (0. - Isohard() * HardSoft() * aI +
                                  (Infyield() * (-HardSoft()) - Inityield() * (-YieldSoft())) *
                                      (1. - exp(-Expisohard() * aI))) *
                              delta_alpha_i_[gp];
        switch (DisMode())
        {
          case Inpar::TSI::pl_multiplier:
            dPlHeatingDT += -delta_alpha_i_[gp] *
                            (0. + Inityield() * YieldSoft() + Isohard() * HardSoft() * aI +
                                (Infyield() * HardSoft() - Inityield() * YieldSoft()) *
                                    (1. - exp(-Expisohard() * aI)));
            break;
          case Inpar::TSI::pl_flow:
            // do nothing
            break;
          default:
            FOUR_C_THROW("unknown plastic dissipation mode: %d", DisMode());
            break;
        }

        // derivative w.r.t. Delta alpha i
        double dPlHeatingDdai =
            (*temp) * (0. - Isohard() * HardSoft() * aI +
                          (Infyield() * (-HardSoft()) - Inityield() * (-YieldSoft())) *
                              (1. - exp(-Expisohard() * aI))) +
            (*temp) * delta_alpha_i_[gp] *
                (0. - Isohard() * HardSoft() +
                    (-Infyield() * HardSoft() + Inityield() * YieldSoft()) * Expisohard() *
                        exp(-Expisohard() * aI));
        switch (DisMode())
        {
          case Inpar::TSI::pl_multiplier:
            dPlHeatingDdai +=
                +Inityield() * (1. - YieldSoft() * dT) +
                Isohard() * (1. - HardSoft() * dT) *
                    (last_alpha_isotropic_[gp] + 2. * delta_alpha_i_[gp]) +
                (Infyield() * (1. - HardSoft() * dT) - Inityield() * (1. - YieldSoft() * dT)) *
                    ((1. - exp(-Expisohard() * aI)) +
                        delta_alpha_i_[gp] * Expisohard() * exp(-Expisohard() * aI));
            break;
          case Inpar::TSI::pl_flow:
            // do nothing
            break;
          default:
            FOUR_C_THROW("unknown plastic dissipation mode: %d", DisMode());
            break;
        }

        // this factor is from the evolution equation for delta_alpha_i
        dPlHeatingDdai *= sq;

        // derivative w.r.t. eta
        Core::LinAlg::Matrix<6, 1> dHpDeta(true);
        if (dDpHeta > 0.)
        {
          tmp61.Multiply(PlAniso_full_, eta_v_strainlike);
          dHpDeta.Update(dPlHeatingDdai * dDpHeta / (abseta_H * absHeta * absHeta), tmp61, 1.);
          dHpDeta.Update(dPlHeatingDdai * abseta_H / (absHeta * absHeta), HdDp_strainlike, 1.);
          dHpDeta.Update(
              -2. * dPlHeatingDdai * abseta_H * dDpHeta / (pow(absHeta, 4.)), HetaH_strainlike, 1.);
        }

        if (DisMode() == Inpar::TSI::pl_flow) dHpDeta.Update(1., deltaDp_v_strainlike, 1.);

        // derivative w.r.t. C
        dHdC->MultiplyTN(*dMdC, dHpDeta);

        // derivative w.r.t. Delta Dp
        dHdDp->MultiplyTN(detaddp, dHpDeta);
        if (dDpHeta > 0.)
        {
          tmp61.Multiply(PlAniso_full_, eta_v_strainlike);
          dHdDp->Update(dPlHeatingDdai * abseta_H / (absHeta * absHeta), tmp61, 1.);
        }
        if (DisMode() == Inpar::TSI::pl_flow) dHdDp->Update(1., eta_v_strainlike, 1.);

        // scaling with time step
        plHeating /= dt;
        dPlHeatingDT /= dt;
        dHdC->Scale(1. / dt);
        dHdDp->Scale(1. / dt);

        // communicate to the element via params (not nice)
        HepDiss(gp) += plHeating;
        (*dHepDissdT_)[gp] += dPlHeatingDT;
      }
    }  // TSI

    // plastic gp
    if (*active)
    {
      // this is a plastic gp
      *elast = false;

      // derivative of the complementarity function w.r.t. to the mandel stress tensor
      Core::LinAlg::Matrix<6, 6> dNCPdeta;
      dNCPdeta.Update(1. - ypl / absetatr_H, pdev, 1.);
      tmp61.Multiply(PlAniso_full_, etatr_v_strainlike);
      dNCPdeta.MultiplyNT(1. / (absetatr_H * absetatr_H), eta_s_v, tmp61, 1.);
      if (dDpHeta > 0.)
      {
        tmp61.Multiply(PlAniso_full_, eta_v_strainlike);
        dNCPdeta.MultiplyNT(
            -dYplDai * dDpHeta / (abseta_H * absHeta * absHeta * absetatr_H), etatr_v, tmp61, 1.);
        dNCPdeta.MultiplyNT(
            -dYplDai * abseta_H / (absHeta * absHeta * absetatr_H), etatr_v, HdDp_strainlike, 1.);
        dNCPdeta.MultiplyNT(2. * dYplDai * abseta_H * dDpHeta / (pow(absHeta, 4.) * absetatr_H),
            etatr_v, HetaH_strainlike, 1.);
      }

      // derivative w.r.t. C
      dNCPdC->Multiply(dNCPdeta, *dMdC);

      // derivative w.r.t. deltaDp
      dNCPdDp->Multiply(dNCPdeta, detaddp);
      Core::LinAlg::Matrix<6, 6> dNCPdetatr;
      tmp61.Multiply(PlAniso_full_, etatr_v_strainlike);
      dNCPdetatr.MultiplyNT(cpl() / (absetatr_H * absetatr_H), eta_s_v, tmp61, 1.);
      dNCPdetatr.Update(-cpl() * ypl / absetatr_H, pdev, 1.);
      dNCPdDp->Multiply(1., dNCPdetatr, InvPlAniso_full_, 1.);
      if (dDpHeta > 0.)
      {
        tmp61.Multiply(PlAniso_full_, eta_v_strainlike);
        dNCPdDp->MultiplyNT(
            -dYplDai * abseta_H / (absetatr_H * absHeta * absHeta), etatr_v, tmp61, 1.);
      }

      // residual
      NCP->Update(eta_v);
      NCP->Update(-ypl / absetatr_H, etatr_v, 1.);

      // derivative w.r.t. temperature
      if (dNCPdT != nullptr) dNCPdT->Update(-dYpldT / absetatr_H, etatr_v, 0.);
    }

    // not active but needs condensation due to acitivity in last iteration
    else if (dDpHeta > 0.)
    {
      // this is an elastic gp, but the flag "elast" is reserved for those
      // elastic GP that do not require a condensation
      *elast = false;

      // residual
      NCP->Multiply(-cpl(), InvPlAniso_full_, deltaDp_v, 1.);

      // derivative of the complementarity function w.r.t. to the mandel stress tensor
      Core::LinAlg::Matrix<6, 6> dNCPdeta;

      tmp61.Multiply(PlAniso_full_, eta_v_strainlike);
      dNCPdeta.MultiplyNT(
          -s() * dYplDai * dDpHeta / (abseta_H * absHeta * absHeta * ypl), *NCP, tmp61, 1.);
      dNCPdeta.MultiplyNT(
          -s() * dYplDai * abseta_H / (absHeta * absHeta * ypl), *NCP, HdDp_strainlike, 1.);
      dNCPdeta.MultiplyNT(s() * 2. * dYplDai * abseta_H * dDpHeta / (pow(absHeta, 4.) * ypl), *NCP,
          HetaH_strainlike, 1.);

      // derivative w.r.t. C
      dNCPdC->Multiply(dNCPdeta, *dMdC);

      // derivative w.r.t. deltaDp
      dNCPdDp->Update(-cpl(), InvPlAniso_full_, 1.);
      dNCPdDp->Multiply(1., dNCPdeta, detaddp, 1.);
      tmp61.Multiply(PlAniso_full_, eta_v_strainlike);
      dNCPdDp->MultiplyNT(dYplDai / (ypl * absHeta), *NCP, tmp61, 1.);

      // derivative w.r.t. temperature
      if (dNCPdT != nullptr)
        dNCPdT->Multiply(-s() * cpl() / ypl * dYpldT, InvPlAniso_full_, deltaDp_v, 0.);
    }
    else
    {
      // Cpl = cpl* Delta D^p = 0
      // we don't build the matrix blocks here.
      // The trivial identity is enforced at the element
      *elast = true;
      // todo: currently, the coupling term dypldT is neglected
      // in the elastic realm. this is a slightly inconsistent linearization for s!=0
      // however, that way we can ensure that deltaLp=0 (exactly) at any Newton iteration

      if (dNCPdT != nullptr) dNCPdT->Clear();
    }
  }
  // elastic gp
  else
  {
    *elast = true;
  }

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate plastic stress and stiffness (with pl. spin)   seitz 05/14 |
 *----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::EvaluatePlast(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<3, 3>* deltaLp, const double* temp, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 9>* dPK2dLp, Core::LinAlg::Matrix<9, 1>* NCP,
    Core::LinAlg::Matrix<9, 6>* dNCPdC, Core::LinAlg::Matrix<9, 9>* dNCPdLp, bool* active,
    bool* elast, bool* as_converged, const int gp, Core::LinAlg::Matrix<9, 1>* dNCPdT,
    Core::LinAlg::Matrix<6, 1>* dHdC, Core::LinAlg::Matrix<9, 1>* dHdLp, const double dt,
    const int eleGID, Core::LinAlg::Matrix<6, 1>* cauchy, Core::LinAlg::Matrix<6, 9>* d_cauchy_ddp,
    Core::LinAlg::Matrix<6, 6>* d_cauchy_dC, Core::LinAlg::Matrix<6, 9>* d_cauchy_dF,
    Core::LinAlg::Matrix<6, 1>* d_cauchy_dT)
{
  int check = +(cauchy != nullptr) + (d_cauchy_dC != nullptr) + (d_cauchy_dF != nullptr) +
              (d_cauchy_ddp != nullptr) + (d_cauchy_dT != nullptr);
  if (!(check == 0 || check == 5)) FOUR_C_THROW("some inconsistency with provided variables");

  Core::LinAlg::Matrix<3, 1> dPI(true);
  Core::LinAlg::Matrix<6, 1> ddPII(true);
  Core::LinAlg::Matrix<3, 1> gamma(true);
  Core::LinAlg::Matrix<8, 1> delta(true);

  if (evaluate_kin_quant_plast(defgrd, deltaLp, gp, params)) return;

  ElastHyperEvaluateInvariantDerivatives(
      prinv_, dPI, ddPII, potsum_, summandProperties_, gp, eleGID);
  CalculateGammaDelta(gamma, delta, prinv_, dPI, ddPII);

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  dPK2dLp->Clear();
  NCP->Clear();
  dNCPdC->Clear();
  dNCPdLp->Clear();
  if (dNCPdT != nullptr) dNCPdT->Clear();

  // new temporary matrices
  Core::LinAlg::Matrix<3, 3> mStr;  // Mandel stress tensor
  Core::LinAlg::Matrix<6, 6> dMdC;  // derivative of Mandel stress w.r.t. RCG
  Core::LinAlg::Matrix<6, 9>
      dMdFpinv;  // derivative of Mandel stress w.r.t. inverse plastic deformation gradient
  Core::LinAlg::Matrix<6, 9> dPK2dFpinv;
  Core::LinAlg::Matrix<6, 9> d_cauchy_dFpi;

  // isotropic elasticity in coupled strain energy format
  // isotropic elasticity in decoupled ("mod") format go here as well
  // as the modified gammas and deltas have been converted
  if (summandProperties_.isoprinc || summandProperties_.isomod)
  {
    evaluate_isotropic_princ_plast(dPK2dFpinv, mStr, dMdC, dMdFpinv, gamma, delta);
  }
  else
    FOUR_C_THROW("only isotropic hypereleastic materials");

  if (cauchy)
    EvaluateCauchyPlast(dPI, ddPII, defgrd, *cauchy, d_cauchy_dFpi, *d_cauchy_dC, *d_cauchy_dF);

  evaluate_nc_pand_spin(&mStr, &dMdC, &dMdFpinv, &dPK2dFpinv, deltaLp, gp, NCP, dNCPdC, dNCPdLp,
      dPK2dLp, active, elast, as_converged, dt, &d_cauchy_dFpi, d_cauchy_ddp);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate NCP function and plastic spin equation         seitz 05/14 |
 *----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::evaluate_nc_pand_spin(const Core::LinAlg::Matrix<3, 3>* mStr,
    const Core::LinAlg::Matrix<6, 6>* dMdC, const Core::LinAlg::Matrix<6, 9>* dMdFpinv,
    const Core::LinAlg::Matrix<6, 9>* dPK2dFpinv, const Core::LinAlg::Matrix<3, 3>* deltaLp,
    const int gp, Core::LinAlg::Matrix<9, 1>* NCP, Core::LinAlg::Matrix<9, 6>* dNCPdC,
    Core::LinAlg::Matrix<9, 9>* dNCPdLp, Core::LinAlg::Matrix<6, 9>* dPK2dLp, bool* active,
    bool* elast, bool* as_converged, const double dt,
    const Core::LinAlg::Matrix<6, 9>* d_cauchy_dFpi, Core::LinAlg::Matrix<6, 9>* d_cauchy_ddp)
{
  const double sq = sqrt(2. / 3.);
  Core::LinAlg::Matrix<6, 1> tmp61;

  // deviatoric projection tensor
  Core::LinAlg::Matrix<6, 6> pdev(true);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (i == j)
        pdev(i, j) = +2. / 3.;
      else
        pdev(i, j) = -1. / 3.;
  for (int i = 3; i < 6; i++) pdev(i, i) = 1.;

  // deviatoric symmetric projection tensor (A-->dev(sym(A))
  Core::LinAlg::Matrix<6, 9> psymdev(true);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (i == j)
        psymdev(i, j) = +2. / 3.;
      else
        psymdev(i, j) = -1. / 3.;
  for (int i = 3; i < 6; i++) psymdev(i, i) = psymdev(i, i + 3) = .5;

  // symmetric identity
  Core::LinAlg::Matrix<6, 9> psym(true);
  for (int i = 0; i < 3; i++) psym(i, i) = 1.;
  for (int i = 3; i < 6; i++) psym(i, i) = psym(i, i + 3) = .5;


  // effective stress
  Core::LinAlg::Matrix<3, 3> eta(*mStr);
  for (int i = 0; i < 3; i++)
    eta(i, i) -= 1. / 3. * ((*mStr)(0, 0) + (*mStr)(1, 1) + (*mStr)(2, 2));
  eta.Update(2. / 3. * Kinhard(), last_alpha_kinematic_[gp], 1.);
  eta.Update(-1. / 3. * Kinhard(), *deltaLp, 1.);
  eta.UpdateT(-1. / 3. * Kinhard(), *deltaLp, 1.);

  // in stress-like voigt notation
  Core::LinAlg::Matrix<6, 1> eta_v;      // in stress-like voigt notation
  Core::LinAlg::Matrix<6, 1> etatr_v;    // in stress-like voigt notation
  Core::LinAlg::Matrix<6, 1> eta_s_v;    // in stress-like voigt notation
  Core::LinAlg::Matrix<6, 1> deltaDp_v;  // in stress-like voigt notation
  for (int i = 0; i < 3; i++)
  {
    eta_v(i) = eta(i, i);
    deltaDp_v(i) = (*deltaLp)(i, i);
  }
  eta_v(3) = .5 * (eta(1, 0) + eta(0, 1));
  deltaDp_v(3) = .5 * ((*deltaLp)(0, 1) + (*deltaLp)(1, 0));
  eta_v(4) = .5 * (eta(1, 2) + eta(2, 1));
  deltaDp_v(4) = .5 * ((*deltaLp)(2, 1) + (*deltaLp)(1, 2));
  eta_v(5) = .5 * (eta(2, 0) + eta(0, 2));
  deltaDp_v(5) = .5 * ((*deltaLp)(0, 2) + (*deltaLp)(2, 0));

  // trial effective stress
  etatr_v.Update(eta_v);
  etatr_v.Multiply(cpl(), InvPlAniso_full_, deltaDp_v, 1.);

  // in strain-like voigt notation
  Core::LinAlg::Matrix<6, 1> eta_v_strainlike(eta_v);          // in strain-like voigt notation
  Core::LinAlg::Matrix<6, 1> etatr_v_strainlike(etatr_v);      // in strain-like voigt notation
  Core::LinAlg::Matrix<6, 1> deltaDp_v_strainlike(deltaDp_v);  // in strain-like voigt notation
  for (int i = 3; i < 6; i++)
  {
    eta_v_strainlike(i) *= 2.;
    etatr_v_strainlike(i) *= 2.;
    deltaDp_v_strainlike(i) *= 2.;
  }

  tmp61.Multiply(PlAniso_full_, eta_v);
  double absHeta = norm_stress_like(tmp61);
  double abseta_H = tmp61.Dot(eta_v_strainlike);
  if (abseta_H < -1.e-16)
    FOUR_C_THROW("this should not happen. tmp=%f", abseta_H);
  else if (abseta_H >= 0.)
    abseta_H = sqrt(abseta_H);
  else
    FOUR_C_THROW("this should not happen. tmp=%f", abseta_H);
  double dDpHeta = tmp61.Dot(deltaDp_v_strainlike);
  tmp61.Multiply(PlAniso_full_, etatr_v);
  double absetatr_H = tmp61.Dot(etatr_v_strainlike);
  if (absetatr_H < -1.e-16)
    FOUR_C_THROW("this should not happen. tmp=%f", absetatr_H);
  else if (absetatr_H >= 0.)
    absetatr_H = sqrt(absetatr_H);
  else
    FOUR_C_THROW("this should not happen. tmp=%f", absetatr_H);
  Core::LinAlg::Matrix<6, 1> HdDp;
  HdDp.Multiply(PlAniso_full_, deltaDp_v);
  Core::LinAlg::Matrix<6, 1> HdDp_strainlike;
  HdDp_strainlike.Multiply(PlAniso_full_, deltaDp_v_strainlike);
  Core::LinAlg::Matrix<6, 1> HetaH_strainlike;
  tmp61.Multiply(PlAniso_full_, eta_v_strainlike);
  HetaH_strainlike.Multiply(PlAniso_full_, tmp61);

  // isotropic hardening increment
  delta_alpha_i_[gp] = 0.;
  if (dDpHeta > 0. && absHeta > 0.)
    delta_alpha_i_[gp] = sq * dDpHeta * abseta_H / (absHeta * absHeta);

  // new isotropic hardening value
  const double aI = last_alpha_isotropic_[gp] + delta_alpha_i_[gp];

  // current yield stress equivalent (yield stress scaled by sqrt(2/3))
  double ypl =
      sq *
      ((Infyield() - Inityield()) * (1. - exp(-Expisohard() * aI)) + Isohard() * aI + Inityield()) *
      pow(1. + Visc() * delta_alpha_i_[gp] / dt, ViscRate());

  // check activity state
  if (ypl < absetatr_H)
  {
    if (activity_state_[gp] == false)  // gp switches state
    {
      if (abs(ypl - absetatr_H) > AS_CONVERGENCE_TOL * Inityield() ||
          deltaLp->NormInf() > AS_CONVERGENCE_TOL * Inityield() / cpl())
        *as_converged = false;
    }
    activity_state_[gp] = true;
    *active = true;
  }
  else
  {
    if (activity_state_[gp] == true)  // gp switches state
    {
      if (abs(ypl - absetatr_H) > AS_CONVERGENCE_TOL * Inityield() ||
          deltaLp->NormInf() > AS_CONVERGENCE_TOL * Inityield() / cpl())
        *as_converged = false;
    }
    activity_state_[gp] = false;
    *active = false;
  }

  // check, if gp is at the corner of the root of the NCP function
  if (abs(ypl - absetatr_H) > AS_CONVERGENCE_TOL * Inityield() ||
      deltaLp->NormInf() > AS_CONVERGENCE_TOL * Inityield() / cpl())
  { /* gp not at the corner point */
  }
  else
  {
    /* gp at the corner point --> elastic and plastic branch are equally valid, take the elastic
     * one*/
    *active = false;
    activity_state_[gp] = false;
  }

  // these cases have some terms in common
  if (*active || dDpHeta > 0. || deltaLp->NormInf() > 0.)
  {
    // derivative of the NCP function w.r.t. RCG / Delta Lp
    // without the lines corresponding to the plastic spin
    Core::LinAlg::Matrix<6, 6> dNCPdC_red;
    Core::LinAlg::Matrix<6, 9> dNCPdLp_red;
    Core::LinAlg::Matrix<6, 1> NCP_red;

    // damping parameter apl
    double apl = 1.;
    if (ypl / abseta_H < 1.) apl = ypl / abseta_H;

    eta_s_v.Update((1. - s()) * ypl / absetatr_H, etatr_v, 1.);
    eta_s_v.Update(apl * s(), eta_v, 1.);

    // matrix exponential derivative
    Core::LinAlg::Matrix<9, 9> Dexp(false);
    Core::LinAlg::Matrix<3, 3> tmp(*deltaLp);
    tmp.Scale(-1.);
    matrix_exponential_derivative3x3(tmp, Dexp);

    // Derivative of inverse plastic deformation gradient
    Core::LinAlg::Matrix<9, 9> dFpiDdeltaLp(true);
    for (int A = 0; A < 3; A++)
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
          for (int i = 0; i < 9; i++)
            dFpiDdeltaLp(
                Core::LinAlg::Voigt::IndexMappings::non_symmetric_tensor_to_voigt9_index(A, a),
                i) -=
                last_plastic_defgrd_inverse_[gp](A, b) *
                Dexp(Core::LinAlg::Voigt::IndexMappings::non_symmetric_tensor_to_voigt9_index(b, a),
                    i);

    // derivative of mandel stress
    Core::LinAlg::Matrix<6, 9> dMdLp;
    dMdLp.Multiply(*dMdFpinv, dFpiDdeltaLp);
    // we spare the deviatoric projection of the mandel stress derivative to get the effective
    // stress derivative. It is enforced implicitly as the "detadLp" is always contracted with
    // deviatoric tensors.
    Core::LinAlg::Matrix<6, 9> detadLp(dMdLp);
    detadLp.Update(-2. / 3. * Kinhard(), psymdev, 1.);
    dPK2dLp->Multiply(*dPK2dFpinv, dFpiDdeltaLp);
    if (d_cauchy_ddp) d_cauchy_ddp->Multiply(*d_cauchy_dFpi, dFpiDdeltaLp);

    // Factor of derivative of Y^pl w.r.t. delta alpha ^i
    // we have added the factor sqrt(2/3) from delta_alpha_i=sq*... here
    double dYplDai =
        2. / 3. *
        (+Isohard() + (Infyield() - Inityield()) * Expisohard() * exp(-Expisohard() * aI)) *
        pow(1. + Visc() * delta_alpha_i_[gp] / dt, ViscRate());
    dYplDai += 2. / 3. *
               ((Infyield() - Inityield()) * (1. - exp(-Expisohard() * aI)) + Isohard() * aI +
                   Inityield()) *
               pow(1. + Visc() * delta_alpha_i_[gp] / dt, ViscRate() - 1) * ViscRate() * Visc() /
               dt;

    // plastic gp
    if (*active)
    {
      // this is a plastic gp
      *elast = false;

      // derivative of the complementarity function w.r.t. to the mandel stress tensor
      Core::LinAlg::Matrix<6, 6> dNCPdeta;
      dNCPdeta.Update(1. - ypl / absetatr_H, pdev, 1.);
      tmp61.Multiply(PlAniso_full_, etatr_v_strainlike);
      dNCPdeta.MultiplyNT(1. / (absetatr_H * absetatr_H), eta_s_v, tmp61, 1.);
      if (dDpHeta > 0.)
      {
        tmp61.Multiply(PlAniso_full_, eta_v_strainlike);
        dNCPdeta.MultiplyNT(
            -dYplDai * dDpHeta / (abseta_H * absHeta * absHeta * absetatr_H), etatr_v, tmp61, 1.);
        dNCPdeta.MultiplyNT(
            -dYplDai * abseta_H / (absHeta * absHeta * absetatr_H), etatr_v, HdDp_strainlike, 1.);
        dNCPdeta.MultiplyNT(2. * dYplDai * abseta_H * dDpHeta / (pow(absHeta, 4.) * absetatr_H),
            etatr_v, HetaH_strainlike, 1.);
      }

      // derivative w.r.t. C
      dNCPdC_red.Multiply(dNCPdeta, *dMdC);

      // derivative w.r.t. deltaLp
      dNCPdLp_red.Multiply(dNCPdeta, detadLp);
      Core::LinAlg::Matrix<6, 6> dNCPdetatr;
      tmp61.Multiply(PlAniso_full_, etatr_v_strainlike);
      dNCPdetatr.MultiplyNT(cpl() / (absetatr_H * absetatr_H), eta_s_v, tmp61, 1.);
      dNCPdetatr.Update(-cpl() * ypl / absetatr_H, pdev, 1.);
      Core::LinAlg::Matrix<6, 6> dNCPdDp;
      dNCPdDp.Multiply(dNCPdetatr, InvPlAniso_full_);
      if (dDpHeta > 0.)
      {
        tmp61.Multiply(PlAniso_full_, eta_v_strainlike);
        dNCPdDp.MultiplyNT(
            -dYplDai * abseta_H / (absetatr_H * absHeta * absHeta), etatr_v, tmp61, 1.);
      }
      dNCPdLp_red.Multiply(1., dNCPdDp, psym, 1.);

      // residual
      NCP_red.Update(eta_v);
      NCP_red.Update(-ypl / absetatr_H, etatr_v, 1.);
    }

    // not active but needs condensation due to acitivity in last iteration
    else if (dDpHeta > 0.)
    {
      // this is an elastic gp, but the flag "elast" is reserved for those
      // elastic GP that do not require a condensation
      *elast = false;

      // residual
      NCP_red.Multiply(-cpl(), InvPlAniso_full_, deltaDp_v, 1.);

      // derivative of the complementarity function w.r.t. to the mandel stress tensor
      Core::LinAlg::Matrix<6, 6> dNCPdeta;

      tmp61.Multiply(PlAniso_full_, eta_v_strainlike);
      dNCPdeta.MultiplyNT(
          -s() * dYplDai * dDpHeta / (abseta_H * absHeta * absHeta * ypl), NCP_red, tmp61, 1.);
      dNCPdeta.MultiplyNT(
          -s() * dYplDai * abseta_H / (absHeta * absHeta * ypl), NCP_red, HdDp_strainlike, 1.);
      dNCPdeta.MultiplyNT(s() * 2. * dYplDai * abseta_H * dDpHeta / (pow(absHeta, 4.) * ypl),
          NCP_red, HetaH_strainlike, 1.);

      // derivative w.r.t. C
      dNCPdC_red.Multiply(dNCPdeta, *dMdC);

      Core::LinAlg::Matrix<6, 6> dNCPdDp;
      dNCPdDp.Update(-cpl(), InvPlAniso_full_, 1.);
      tmp61.Multiply(PlAniso_full_, eta_v_strainlike);
      dNCPdDp.MultiplyNT(dYplDai / (ypl * absHeta), NCP_red, tmp61, 1.);
      dNCPdLp_red.Multiply(1., dNCPdDp, psym, 1.);
      dNCPdLp_red.Multiply(1., dNCPdeta, detadLp, 1.);
    }
    else
    {
      // Cpl = cpl* Delta D^p = 0
      // we don't build the matrix blocks here.
      // The trivial identity is enforced at the element
      *elast = true;
    }

    // plastic spin equation
    // the tensor product Sigma.Dp is not made for voigt notation so we do it the hard way
    Core::LinAlg::Matrix<3, 1> spEq;
    spEq(0) = .5 * ((*deltaLp)(0, 1) - (*deltaLp)(1, 0)) -
              PlSpinChi() / Inityield() *
                  ((*mStr)(0, 0) * deltaDp_v(3) + (*mStr)(0, 1) * deltaDp_v(1) +
                      (*mStr)(0, 2) * deltaDp_v(4) - (*mStr)(0, 1) * deltaDp_v(0) -
                      (*mStr)(1, 1) * deltaDp_v(3) - (*mStr)(1, 2) * deltaDp_v(5));
    spEq(1) = .5 * ((*deltaLp)(1, 2) - (*deltaLp)(2, 1)) -
              PlSpinChi() / Inityield() *
                  ((*mStr)(0, 1) * deltaDp_v(5) + (*mStr)(1, 1) * deltaDp_v(4) +
                      (*mStr)(1, 2) * deltaDp_v(2) - (*mStr)(0, 2) * deltaDp_v(3) -
                      (*mStr)(1, 2) * deltaDp_v(1) - (*mStr)(2, 2) * deltaDp_v(4));
    spEq(2) = .5 * ((*deltaLp)(0, 2) - (*deltaLp)(2, 0)) -
              PlSpinChi() / Inityield() *
                  ((*mStr)(0, 0) * deltaDp_v(5) + (*mStr)(0, 1) * deltaDp_v(4) +
                      (*mStr)(0, 2) * deltaDp_v(2) - (*mStr)(0, 2) * deltaDp_v(0) -
                      (*mStr)(1, 2) * deltaDp_v(3) - (*mStr)(2, 2) * deltaDp_v(5));

    // Derivative of plastic spin equation w.r.t. mandel stress
    Core::LinAlg::Matrix<3, 6> dSpdM;
    dSpdM(0, 0) = +deltaDp_v(3);
    dSpdM(0, 1) = -deltaDp_v(3);
    dSpdM(0, 2) = 0.;
    dSpdM(0, 3) = +deltaDp_v(1) - deltaDp_v(0);
    dSpdM(0, 4) = -deltaDp_v(5);
    dSpdM(0, 5) = deltaDp_v(4);
    dSpdM(1, 0) = 0.;
    dSpdM(1, 1) = deltaDp_v(4);
    dSpdM(1, 2) = -deltaDp_v(4);
    dSpdM(1, 3) = deltaDp_v(5);
    dSpdM(1, 4) = deltaDp_v(2) - deltaDp_v(1);
    dSpdM(1, 5) = -deltaDp_v(3);
    dSpdM(2, 0) = deltaDp_v(5);
    dSpdM(2, 1) = 0.;
    dSpdM(2, 2) = -deltaDp_v(5);
    dSpdM(2, 3) = deltaDp_v(4);
    dSpdM(2, 4) = -deltaDp_v(3);
    dSpdM(2, 5) = deltaDp_v(2) - deltaDp_v(0);
    dSpdM.Scale(-PlSpinChi() / Inityield());

    // derivative of plastic spin equation w.r.t. deltaDp
    Core::LinAlg::Matrix<3, 6> dSpdDp;
    dSpdDp(0, 0) = -(*mStr)(0, 1);
    dSpdDp(0, 1) = (*mStr)(0, 1);
    dSpdDp(0, 2) = 0.;
    dSpdDp(0, 3) = (*mStr)(0, 0) - (*mStr)(1, 1);
    dSpdDp(0, 4) = +(*mStr)(0, 2);
    dSpdDp(0, 5) = -(*mStr)(1, 2);
    dSpdDp(1, 0) = 0.;
    dSpdDp(1, 1) = -(*mStr)(1, 2);
    dSpdDp(1, 2) = +(*mStr)(1, 2);
    dSpdDp(1, 3) = -(*mStr)(0, 2);
    dSpdDp(1, 4) = +(*mStr)(1, 1) - (*mStr)(2, 2);
    dSpdDp(1, 5) = (*mStr)(0, 1);
    dSpdDp(2, 0) = -(*mStr)(0, 2);
    dSpdDp(2, 1) = 0.;
    dSpdDp(2, 2) = +(*mStr)(0, 2);
    dSpdDp(2, 3) = -(*mStr)(1, 2);
    dSpdDp(2, 4) = +(*mStr)(0, 1);
    dSpdDp(2, 5) = (*mStr)(0, 0) - (*mStr)(2, 2);
    dSpdDp.Scale(-PlSpinChi() / Inityield());

    // derivative of plastic spin equation w.r.t. RCG
    Core::LinAlg::Matrix<3, 6> dSpdC;
    dSpdC.Multiply(dSpdM, *dMdC);

    // derivative of plastic spin equation w.r.t. deltaLp
    Core::LinAlg::Matrix<3, 9> dSpddLp;
    dSpddLp.Multiply(dSpdDp, psym);
    dSpddLp.Multiply(1., dSpdM, dMdLp, 1.);
    dSpddLp(0, 3) += .5;
    dSpddLp(0, 6) -= .5;
    dSpddLp(1, 4) += .5;
    dSpddLp(1, 7) -= .5;
    dSpddLp(2, 5) += .5;
    dSpddLp(2, 8) -= .5;

    // combine NCP and plastic spin equations
    for (int i = 0; i < 6; i++) (*NCP)(i) = NCP_red(i);
    for (int i = 6; i < 9; i++) (*NCP)(i) = spEq(i - 6);

    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++) (*dNCPdC)(i, j) = dNCPdC_red(i, j);
    for (int i = 6; i < 9; i++)
      for (int j = 0; j < 6; j++) (*dNCPdC)(i, j) = dSpdC(i - 6, j);

    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 9; j++) (*dNCPdLp)(i, j) = dNCPdLp_red(i, j);
    for (int i = 6; i < 9; i++)
      for (int j = 0; j < 9; j++) (*dNCPdLp)(i, j) = dSpddLp(i - 6, j);
  }
  // elastic gp
  else
  {
    *elast = true;
  }

  return;
}

void Mat::PlasticElastHyper::EvaluateCauchyPlast(const Core::LinAlg::Matrix<3, 1>& dPI,
    const Core::LinAlg::Matrix<6, 1>& ddPII, const Core::LinAlg::Matrix<3, 3>* defgrd,
    Core::LinAlg::Matrix<6, 1>& cauchy, Core::LinAlg::Matrix<6, 9>& d_cauchy_dFpi,
    Core::LinAlg::Matrix<6, 6>& d_cauchy_dC, Core::LinAlg::Matrix<6, 9>& d_cauchy_dF,
    Core::LinAlg::Matrix<6, 1>* d_cauchy_dT)
{
  cauchy.Clear();
  d_cauchy_dC.Clear();
  d_cauchy_dF.Clear();
  d_cauchy_dFpi.Clear();

  cauchy.Update(sqrt(prinv_(2)) * dPI(2), id2V_, 1.);
  cauchy.Update((dPI(0) + prinv_(0) * dPI(1)) / sqrt(prinv_(2)), bev_, 1.);
  cauchy.Update(-dPI(1) / sqrt(prinv_(2)), be2v_, 1.);
  cauchy.Scale(2.);

  d_cauchy_dC.MultiplyNT(sqrt(prinv_(2)) * (ddPII(4) + prinv_(0) * ddPII(3)), id2V_, Cpi_, 1.);
  d_cauchy_dC.MultiplyNT(-sqrt(prinv_(2)) * ddPII(3), id2V_, CpiCCpi_, 1.);
  d_cauchy_dC.MultiplyNT(sqrt(prinv_(2)) * (.5 * dPI(2) + prinv_(2) * ddPII(2)), id2V_, ircg_, 1.);

  d_cauchy_dC.MultiplyNT(
      (ddPII(0) + dPI(1) + 2. * prinv_(0) * ddPII(5) + prinv_(0) * prinv_(0) * ddPII(1)) /
          sqrt(prinv_(2)),
      bev_, Cpi_, 1.);
  d_cauchy_dC.MultiplyNT((-ddPII(5) - prinv_(0) * ddPII(1)) / sqrt(prinv_(2)), bev_, CpiCCpi_, 1.);
  d_cauchy_dC.MultiplyNT((-.5 * dPI(0) - .5 * prinv_(0) * dPI(1) + prinv_(2) * ddPII(4) +
                             prinv_(0) * prinv_(2) * ddPII(3)) /
                             sqrt(prinv_(2)),
      bev_, ircg_, 1.);

  d_cauchy_dC.MultiplyNT(-(ddPII(5) - prinv_(0) * ddPII(1)) / sqrt(prinv_(2)), be2v_, Cpi_, 1.);
  d_cauchy_dC.MultiplyNT(ddPII(1) / sqrt(prinv_(2)), be2v_, CpiCCpi_, 1.);
  d_cauchy_dC.MultiplyNT((.5 * dPI(1) / prinv_(2) - ddPII(3)) / sqrt(prinv_(2)), be2v_, ircg_, 1.);
  d_cauchy_dC.Scale(4.);

  d_cauchy_dFpi.MultiplyNT(sqrt(prinv_(2)) * (ddPII(4) + prinv_(0) * ddPII(3)), id2V_, CFpi_, 1.);
  d_cauchy_dFpi.MultiplyNT(-sqrt(prinv_(2)) * ddPII(3), id2V_, CFpiCe_, 1.);
  d_cauchy_dFpi.MultiplyNT(
      sqrt(prinv_(2)) * (.5 * dPI(2) + prinv_(2) * ddPII(2)), id2V_, CFpiCei_, 1.);

  d_cauchy_dFpi.MultiplyNT(
      (ddPII(0) + dPI(1) + 2. * prinv_(0) * ddPII(5) + prinv_(0) * prinv_(0) * ddPII(1)) /
          sqrt(prinv_(2)),
      bev_, CFpi_, 1.);
  d_cauchy_dFpi.MultiplyNT((-ddPII(5) - prinv_(0) * ddPII(1)) / sqrt(prinv_(2)), bev_, CFpiCe_, 1.);
  d_cauchy_dFpi.MultiplyNT((-.5 * dPI(0) - .5 * prinv_(0) * dPI(1) + prinv_(2) * ddPII(4) +
                               prinv_(0) * prinv_(2) * ddPII(3)) /
                               sqrt(prinv_(2)),
      bev_, CFpiCei_, 1.);

  d_cauchy_dFpi.MultiplyNT(-(ddPII(5) - prinv_(0) * ddPII(1)) / sqrt(prinv_(2)), be2v_, CFpi_, 1.);
  d_cauchy_dFpi.MultiplyNT(ddPII(1) / sqrt(prinv_(2)), be2v_, CFpiCe_, 1.);
  d_cauchy_dFpi.MultiplyNT(
      (.5 * dPI(1) / prinv_(2) - ddPII(3)) / sqrt(prinv_(2)), be2v_, CFpiCei_, 1.);
  d_cauchy_dFpi.Scale(2.);

  AddRightNonSymmetricHolzapfelProduct(
      d_cauchy_dFpi, *defgrd, Fe_, (dPI(0) + prinv_(0) * dPI(1)) / sqrt(prinv_(2)));
  AddRightNonSymmetricHolzapfelProduct(d_cauchy_dFpi, *defgrd, beFe_, -dPI(1) / sqrt(prinv_(2)));
  AddRightNonSymmetricHolzapfelProduct(d_cauchy_dFpi, beF_, Fe_, -dPI(1) / sqrt(prinv_(2)));

  AddRightNonSymmetricHolzapfelProduct(
      d_cauchy_dF, id2_, FCpi_, (dPI(0) + prinv_(0) * dPI(1)) / sqrt(prinv_(2)));
  AddRightNonSymmetricHolzapfelProduct(d_cauchy_dF, id2_, beFCpi_, -dPI(1) / sqrt(prinv_(2)));
  AddRightNonSymmetricHolzapfelProduct(d_cauchy_dF, be_, FCpi_, -dPI(1) / sqrt(prinv_(2)));
  d_cauchy_dF.Scale(2.);
  d_cauchy_dFpi.Scale(2.);

  if (d_cauchy_dT)
  {
    d_cauchy_dT->Clear();

    const double j = sqrt(prinv_(2));
    double d_dPI2_dT = 0.;
    if (summandProperties_.isomod) FOUR_C_THROW("only coupled SEF are supposed to end up here");
    for (unsigned int p = 0; p < potsum_.size(); ++p)
      potsum_[p]->AddCoupDerivVol(j, nullptr, &d_dPI2_dT, nullptr, nullptr);

    const double fac = -3. * Cte();
    d_dPI2_dT *= fac * .5 / j;

    d_cauchy_dT->Update(sqrt(prinv_(2)) * d_dPI2_dT, id2V_, 1.);
    d_cauchy_dT->Scale(2.);
  }
}

void Mat::PlasticElastHyper::evaluate_cauchy_derivs(const Core::LinAlg::Matrix<3, 1>& prinv,
    const int gp, int eleGID, Core::LinAlg::Matrix<3, 1>& dPI, Core::LinAlg::Matrix<6, 1>& ddPII,
    Core::LinAlg::Matrix<10, 1>& dddPIII, const double* temp)
{
  for (unsigned i = 0; i < potsum_.size(); ++i)
  {
    if (summandProperties_.isoprinc)
    {
      potsum_[i]->add_derivatives_principal(dPI, ddPII, prinv, gp, eleGID);
      potsum_[i]->add_third_derivatives_principal_iso(dddPIII, prinv, gp, eleGID);
    }
    if (summandProperties_.isomod || summandProperties_.anisomod || summandProperties_.anisoprinc)
      FOUR_C_THROW("not implemented for this form of strain energy function");
  }
  if (temp)
  {
    if (summandProperties_.isomod || summandProperties_.anisomod || summandProperties_.anisoprinc)
      FOUR_C_THROW(
          "thermo-elastic Nitsche contact only with strain energies"
          "\ndepending on unmodified invariants");

    const double dT = *temp - InitTemp();

    const double j = sqrt(prinv(2));
    double dPmodI = 0.;
    double ddPmodII = 0.;
    double dddPmodIII = 0.;
    double ddddPmodIIII = 0.;
    for (unsigned int p = 0; p < potsum_.size(); ++p)
      potsum_[p]->AddCoupDerivVol(j, &dPmodI, &ddPmodII, &dddPmodIII, &ddddPmodIIII);

    const double fac = -3. * Cte() * dT;
    dPI(2) += fac * .5 / j * ddPmodII;
    ddPII(2) += fac * (.25 / (j * j) * dddPmodIII - .25 / (j * j * j) * ddPmodII);
    dddPIII(2) +=
        fac * (+1. / (8. * j * j * j) * ddddPmodIIII - 3. / (8. * j * j * j * j) * dddPmodIII +
                  3. / (8. * j * j * j * j * j) * ddPmodII);
  }
}

void Mat::PlasticElastHyper::evaluate_cauchy_temp_deriv(const Core::LinAlg::Matrix<3, 1>& prinv,
    const double ndt, const double bdndt, const double ibdndt, const double* temp, double* DsntDT,
    const Core::LinAlg::Matrix<9, 1>& iFTV, const Core::LinAlg::Matrix<9, 1>& DbdndtDFV,
    const Core::LinAlg::Matrix<9, 1>& DibdndtDFV, const Core::LinAlg::Matrix<9, 1>& DI1DF,
    const Core::LinAlg::Matrix<9, 1>& DI2DF, const Core::LinAlg::Matrix<9, 1>& DI3DF,
    Core::LinAlg::Matrix<9, 1>* D2sntDFDT)
{
  const double sqI3 = sqrt(prinv(2));
  const double prefac = 2.0 / sqI3;
  double dPmodI = 0.0;
  double ddPmodII = 0.0;
  double dddPmodIII = 0.0;
  for (unsigned int p = 0; p < potsum_.size(); ++p)
    potsum_[p]->AddCoupDerivVol(sqI3, &dPmodI, &ddPmodII, &dddPmodIII, nullptr);

  // those are actually the temperature derivatives of the coefficients
  // then we plug them into the cauchy stress derivative and voila
  // we keep many zero entries here for possible future extension to more
  // than just thermal expansion
  static Core::LinAlg::Matrix<3, 1> dPI(true);
  static Core::LinAlg::Matrix<6, 1> ddPII(true);
  dPI.Clear();
  ddPII.Clear();

  const double fac = -3.0 * Cte();
  dPI(2) += fac * 0.5 / sqI3 * ddPmodII;
  ddPII(2) += fac * (0.25 / (sqI3 * sqI3) * dddPmodIII - 0.25 / (sqI3 * sqI3 * sqI3) * ddPmodII);

  *DsntDT = prefac * (prinv(1) * dPI(1) * ndt + prinv(2) * dPI(2) * ndt + dPI(0) * bdndt -
                         prinv(2) * dPI(1) * ibdndt);

  if (D2sntDFDT)
  {
    D2sntDFDT->Update(-prefac * (prinv(1) * dPI(1) * ndt + prinv(2) * dPI(2) * ndt +
                                    dPI(0) * bdndt - prinv(2) * dPI(1) * ibdndt),
        iFTV, 0.0);  // D2sntDFDT is cleared here
    D2sntDFDT->Update(prefac * dPI(0), DbdndtDFV, 1.0);
    D2sntDFDT->Update(-prefac * prinv(2) * dPI(1), DibdndtDFV, 1.0);
    D2sntDFDT->Update(prefac * (prinv(1) * ddPII(5) * ndt + prinv(2) * ddPII(4) * ndt +
                                   ddPII(0) * bdndt - prinv(2) * ddPII(5) * ibdndt),
        DI1DF, 1.0);
    D2sntDFDT->Update(
        prefac * (dPI(1) * ndt + prinv(1) * ddPII(1) * ndt + prinv(2) * ddPII(3) * ndt +
                     ddPII(5) * bdndt - prinv(2) * ddPII(1) * ibdndt),
        DI2DF, 1.0);
    D2sntDFDT->Update(
        prefac * (prinv(1) * ddPII(3) * ndt + dPI(2) * ndt + prinv(2) * ddPII(2) * ndt +
                     ddPII(4) * bdndt - dPI(1) * ibdndt - prinv(2) * ddPII(3) * ibdndt),
        DI3DF, 1.0);
  }
}


void Mat::PlasticElastHyper::add_thermal_expansion_derivs(const Core::LinAlg::Matrix<3, 1>& prinv,
    Core::LinAlg::Matrix<3, 1>& dPI, Core::LinAlg::Matrix<6, 1>& ddPII, const int gp, int eleGID,
    const double& temp)
{
  const double j = sqrt(prinv(2));
  Core::LinAlg::Matrix<3, 1> modinv(true);
  modinv(2) = j;
  Core::LinAlg::Matrix<3, 1> dPmodI;
  Core::LinAlg::Matrix<6, 1> ddPmodII;
  double dddPmodIII = 0.;
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->add_derivatives_modified(dPmodI, ddPmodII, modinv, gp, eleGID);
    potsum_[p]->Add3rdVolDeriv(modinv, dddPmodIII);
    potsum_[p]->AddCoupDerivVol(j, nullptr, &(ddPmodII(2)), &dddPmodIII, nullptr);
  }

  const double dT = temp - InitTemp();
  const double fac = -3. * Cte() * dT;
  dPI(2) += fac * .5 / j * ddPmodII(2);
  ddPII(2) += fac * (.25 / (j * j) * dddPmodIII - .25 / (j * j * j) * ddPmodII(2));
}

void Mat::PlasticElastHyper::UpdateGP(const int gp, const Core::LinAlg::Matrix<3, 3>* deltaDp)
{
  if (activity_state_[gp] == true)
  {
    // update plastic deformation gradient
    Core::LinAlg::Matrix<3, 3> tmp;
    tmp.Update(-1., *deltaDp);
    matrix_exponential3x3(tmp);
    Core::LinAlg::Matrix<3, 3> fpi_last = last_plastic_defgrd_inverse_[gp];
    last_plastic_defgrd_inverse_[gp].Multiply(fpi_last, tmp);
    // update isotropic hardening
    last_alpha_isotropic_[gp] += delta_alpha_i_[gp];

    // update kinematic hardening
    last_alpha_kinematic_[gp].Update(-.5, *deltaDp, 1.);
    last_alpha_kinematic_[gp].UpdateT(-.5, *deltaDp, 1.);
  }

  return;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::evaluate_kin_quant_elast(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<3, 3>* deltaLp, const int gp)
{
  Core::LinAlg::Matrix<3, 3> tmp;
  Core::LinAlg::Matrix<3, 3> invpldefgrd;
  Core::LinAlg::Matrix<3, 3>& InvPlasticDefgrdLast = last_plastic_defgrd_inverse_[gp];
  tmp.Update(-1., *deltaLp);
  matrix_exponential3x3(tmp);
  invpldefgrd.Multiply(InvPlasticDefgrdLast, tmp);

  // inverse plastic right Cauchy-Green
  Core::LinAlg::Matrix<3, 3> CpiM;
  CpiM.MultiplyNT(invpldefgrd, invpldefgrd);
  // stress-like Voigt notation
  for (int i = 0; i < 3; i++) Cpi_(i) = CpiM(i, i);
  Cpi_(3) = (CpiM(0, 1) + CpiM(1, 0)) / 2.;
  Cpi_(4) = (CpiM(2, 1) + CpiM(1, 2)) / 2.;
  Cpi_(5) = (CpiM(0, 2) + CpiM(2, 0)) / 2.;

  // inverse RCG
  Core::LinAlg::Matrix<3, 3> iRCG;
  Core::LinAlg::Matrix<3, 3> RCG;
  RCG.MultiplyTN(*defgrd, *defgrd);
  iRCG.Invert(RCG);
  // stress-like Voigt notation
  for (int i = 0; i < 3; i++) ircg_(i) = iRCG(i, i);
  ircg_(3) = (iRCG(0, 1) + iRCG(1, 0)) / 2.;
  ircg_(4) = (iRCG(2, 1) + iRCG(1, 2)) / 2.;
  ircg_(5) = (iRCG(0, 2) + iRCG(2, 0)) / 2.;

  // C_p^-1 * C * C_p^-1
  Core::LinAlg::Matrix<3, 3> CpiCCpiM;
  tmp.Multiply(CpiM, RCG);
  CpiCCpiM.Multiply(tmp, CpiM);
  // stress-like Voigt notation
  for (int i = 0; i < 3; i++) CpiCCpi_(i) = CpiCCpiM(i, i);
  CpiCCpi_(3) = (CpiCCpiM(0, 1) + CpiCCpiM(1, 0)) / 2.;
  CpiCCpi_(4) = (CpiCCpiM(2, 1) + CpiCCpiM(1, 2)) / 2.;
  CpiCCpi_(5) = (CpiCCpiM(0, 2) + CpiCCpiM(2, 0)) / 2.;

  tmp.Multiply(*defgrd, invpldefgrd);
  Core::LinAlg::Matrix<3, 3> CeM;
  CeM.MultiplyTN(tmp, tmp);
  // elastic right Cauchy-Green in strain-like Voigt notation.
  Core::LinAlg::Matrix<6, 1> elasticRCGv;
  for (int i = 0; i < 3; i++) elasticRCGv(i) = CeM(i, i);
  elasticRCGv(3) = (CeM(0, 1) + CeM(1, 0));
  elasticRCGv(4) = (CeM(2, 1) + CeM(1, 2));
  elasticRCGv(5) = (CeM(0, 2) + CeM(2, 0));

  // principal invariants of elastic Cauchy-Green strain
  Core::LinAlg::Voigt::Strains::invariants_principal(prinv_, elasticRCGv);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Mat::PlasticElastHyper::evaluate_kin_quant_plast(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<3, 3>* deltaLp, const int gp, Teuchos::ParameterList& params)
{
  id2_.Clear();
  id2V_.Clear();
  for (int i = 0; i < 3; i++)
  {
    id2V_(i) = 1.;
    id2_(i, i) = 1.;
  }
  Core::LinAlg::Matrix<3, 3> tmp;
  Core::LinAlg::Matrix<3, 3> tmp33;
  Core::LinAlg::Matrix<3, 3>& InvPlasticDefgrdLast = last_plastic_defgrd_inverse_[gp];
  tmp.Update(-1., *deltaLp);
  matrix_exponential3x3(tmp);
  invpldefgrd_.Multiply(InvPlasticDefgrdLast, tmp);

  tmp33.Multiply(*defgrd, invpldefgrd_);
  Fe_.Update(tmp33);
  be_.MultiplyNT(Fe_, Fe_);
  for (int i = 0; i < 3; ++i) bev_(i) = be_(i, i);
  bev_(3) = be_(0, 1);
  bev_(4) = be_(1, 2);
  bev_(5) = be_(0, 2);

  be2_.Multiply(be_, be_);
  for (int i = 0; i < 3; ++i) be2v_(i) = be2_(i, i);
  be2v_(3) = be2_(0, 1);
  be2v_(4) = be2_(1, 2);
  be2v_(5) = be2_(0, 2);

  beF_.Multiply(be_, *defgrd);
  beFe_.Multiply(be_, Fe_);

  CeM_.MultiplyTN(tmp33, tmp33);
  // elastic right Cauchy-Green in strain-like Voigt notation.
  Core::LinAlg::Matrix<6, 1> elasticRCGv;
  for (int i = 0; i < 3; i++) elasticRCGv(i) = CeM_(i, i);
  elasticRCGv(3) = (CeM_(0, 1) + CeM_(1, 0));
  elasticRCGv(4) = (CeM_(2, 1) + CeM_(1, 2));
  elasticRCGv(5) = (CeM_(0, 2) + CeM_(2, 0));
  // elastic right Cauchy-Green in stress-like Voigt notation.
  for (int i = 0; i < 3; i++) Ce_(i) = CeM_(i, i);
  Ce_(3) = (CeM_(0, 1) + CeM_(1, 0)) / 2.;
  Ce_(4) = (CeM_(2, 1) + CeM_(1, 2)) / 2.;
  Ce_(5) = (CeM_(0, 2) + CeM_(2, 0)) / 2.;

  // square of elastic right Cauchy-Green in stress-like Voigt notation.
  tmp.Multiply(CeM_, CeM_);
  for (int i = 0; i < 3; i++) Ce2_(i) = tmp(i, i);
  Ce2_(3) = (tmp(0, 1) + tmp(1, 0)) / 2.;
  Ce2_(4) = (tmp(2, 1) + tmp(1, 2)) / 2.;
  Ce2_(5) = (tmp(0, 2) + tmp(2, 0)) / 2.;

  // principal invariants of elastic Cauchy-Green strain
  Core::LinAlg::Voigt::Strains::invariants_principal(prinv_, elasticRCGv);

  // inverse plastic right Cauchy-Green
  Core::LinAlg::Matrix<3, 3> CpiM;
  CpiM.MultiplyNT(invpldefgrd_, invpldefgrd_);
  // stress-like Voigt notation
  for (int i = 0; i < 3; i++) Cpi_(i) = CpiM(i, i);
  Cpi_(3) = (CpiM(0, 1) + CpiM(1, 0)) / 2.;
  Cpi_(4) = (CpiM(2, 1) + CpiM(1, 2)) / 2.;
  Cpi_(5) = (CpiM(0, 2) + CpiM(2, 0)) / 2.;
  FCpi_.Multiply(*defgrd, CpiM);
  beFCpi_.Multiply(be_, FCpi_);

  // inverse RCG
  Core::LinAlg::Matrix<3, 3> iRCG;
  Core::LinAlg::Matrix<3, 3> RCG;
  RCG.MultiplyTN(*defgrd, *defgrd);
  iRCG.Invert(RCG);
  // stress-like Voigt notation
  for (int i = 0; i < 3; i++) ircg_(i) = iRCG(i, i);
  ircg_(3) = (iRCG(0, 1) + iRCG(1, 0)) / 2.;
  ircg_(4) = (iRCG(2, 1) + iRCG(1, 2)) / 2.;
  ircg_(5) = (iRCG(0, 2) + iRCG(2, 0)) / 2.;

  // C_p^-1 * C * C_p^-1
  Core::LinAlg::Matrix<3, 3> CpiCCpiM;
  tmp33.Multiply(CpiM, RCG);
  CpiCCpiM.Multiply(tmp33, CpiM);
  // stress-like Voigt notation
  for (int i = 0; i < 3; i++) CpiCCpi_(i) = CpiCCpiM(i, i);
  CpiCCpi_(3) = (CpiCCpiM(0, 1) + CpiCCpiM(1, 0)) / 2.;
  CpiCCpi_(4) = (CpiCCpiM(2, 1) + CpiCCpiM(1, 2)) / 2.;
  CpiCCpi_(5) = (CpiCCpiM(0, 2) + CpiCCpiM(2, 0)) / 2.;

  CpiC_.Multiply(CpiM, RCG);
  FpiCe_.Multiply(invpldefgrd_, CeM_);

  FpiTC_.MultiplyTN(invpldefgrd_, RCG);
  CeFpiTC_.Multiply(CeM_, FpiTC_);

  tmp.Multiply(RCG, invpldefgrd_);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(tmp, CFpi_);
  tmp33.Multiply(tmp, CeM_);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(tmp33, CFpiCe_);

  double det = CeM_.Determinant();
  if (det > -1e-30 and det < 1e-30)
    if (params.isParameter("tolerate_errors"))
      if (params.get<bool>("tolerate_errors") == true)
      {
        params.get<bool>("eval_error") = true;
        return 1;
      }

  tmp.Invert(CeM_);
  tmp33.Multiply(invpldefgrd_, tmp);
  tmp.Multiply(RCG, tmp33);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(tmp, CFpiCei_);

  return 0;
}

/*----------------------------------------------------------------------*
 |  tensor norm from stress-like Voigt-notation             seitz 05/14 |
 *----------------------------------------------------------------------*/
double Mat::PlasticElastHyper::norm_stress_like(const Core::LinAlg::Matrix<6, 1>& stress)
{
  double norm = 0;
  for (int i = 0; i < 3; ++i) norm += stress(i) * stress(i);
  for (int i = 3; i < 6; ++i) norm += 2. * stress(i) * stress(i);
  norm = sqrt(norm);
  return norm;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::evaluate_isotropic_princ_elast(
    Core::LinAlg::Matrix<6, 1>& stressisoprinc, Core::LinAlg::Matrix<6, 6>& cmatisoprinc,
    Core::LinAlg::Matrix<3, 1> dPI, Core::LinAlg::Matrix<6, 1> ddPII)
{
  // 2nd Piola Kirchhoff stress (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
  // factors
  Core::LinAlg::Matrix<3, 1> gamma(true);
  gamma(0) = 2. * (dPI(0) + prinv_(0) * dPI(1));
  gamma(1) = -2. * dPI(1);
  gamma(2) = 2. * prinv_(2) * dPI(2);

  //  // 2nd Piola Kirchhoff stresses
  stressisoprinc.Update(gamma(0), Cpi_, 1.0);
  stressisoprinc.Update(gamma(1), CpiCCpi_, 1.0);
  stressisoprinc.Update(gamma(2), ircg_, 1.0);

  // constitutive tensor according to Holzapfel-Nonlinear Solid Mechanics p. 261)
  // factors
  Core::LinAlg::Matrix<8, 1> delta(true);
  delta(0) =
      4. * (ddPII(0) + 2. * prinv_(0) * ddPII(5) + dPI(1) + prinv_(0) * prinv_(0) * ddPII(1));
  delta(1) = -4. * (ddPII(5) + prinv_(0) * ddPII(1));
  delta(2) = 4. * (prinv_(2) * ddPII(4) + prinv_(0) * prinv_(2) * ddPII(3));
  delta(3) = 4. * ddPII(1);
  delta(4) = -4. * prinv_(2) * ddPII(3);
  delta(5) = 4. * (prinv_(2) * dPI(2) + prinv_(2) * prinv_(2) * ddPII(2));
  delta(6) = -4. * prinv_(2) * dPI(2);
  delta(7) = -4. * dPI(1);

  // constitutive tensor
  cmatisoprinc.MultiplyNT(delta(0), Cpi_, Cpi_, 1.);
  cmatisoprinc.MultiplyNT(delta(1), CpiCCpi_, Cpi_, 1.);
  cmatisoprinc.MultiplyNT(delta(1), Cpi_, CpiCCpi_, 1.);
  cmatisoprinc.MultiplyNT(delta(2), Cpi_, ircg_, 1.);
  cmatisoprinc.MultiplyNT(delta(2), ircg_, Cpi_, 1.);
  cmatisoprinc.MultiplyNT(delta(3), CpiCCpi_, CpiCCpi_, 1.);
  cmatisoprinc.MultiplyNT(delta(4), CpiCCpi_, ircg_, 1.);
  cmatisoprinc.MultiplyNT(delta(4), ircg_, CpiCCpi_, 1.);
  cmatisoprinc.MultiplyNT(delta(5), ircg_, ircg_, 1.);
  AddtoCmatHolzapfelProduct(cmatisoprinc, ircg_, delta(6));
  AddtoCmatHolzapfelProduct(cmatisoprinc, Cpi_, delta(7));

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::evaluate_isotropic_princ_plast(
    Core::LinAlg::Matrix<6, 9>& dPK2dFpinvIsoprinc,
    Core::LinAlg::Matrix<3, 3>& MandelStressIsoprinc, Core::LinAlg::Matrix<6, 6>& dMdCisoprinc,
    Core::LinAlg::Matrix<6, 9>& dMdFpinvIsoprinc, const Core::LinAlg::Matrix<3, 1>& gamma,
    const Core::LinAlg::Matrix<8, 1>& delta)
{
  // derivative of PK2 w.r.t. inverse plastic deformation gradient
  AddRightNonSymmetricHolzapfelProduct(dPK2dFpinvIsoprinc, id2_, invpldefgrd_, gamma(0));
  AddRightNonSymmetricHolzapfelProduct(dPK2dFpinvIsoprinc, CpiC_, invpldefgrd_, gamma(1));
  dPK2dFpinvIsoprinc.MultiplyNT(delta(0), Cpi_, CFpi_, 1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(1), Cpi_, CFpiCe_, 1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(1), CpiCCpi_, CFpi_, 1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(2), Cpi_, CFpiCei_, 1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(2), ircg_, CFpi_, 1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(3), CpiCCpi_, CFpiCe_, 1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(4), CpiCCpi_, CFpiCei_, 1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(4), ircg_, CFpiCe_, 1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(5), ircg_, CFpiCei_, 1.);
  AddRightNonSymmetricHolzapfelProduct(dPK2dFpinvIsoprinc, id2_, FpiCe_, 0.5 * delta(7));

  // Mandel stress
  Core::LinAlg::Matrix<6, 1> Mv;
  Mv.Update(gamma(0), Ce_);
  Mv.Update(gamma(1), Ce2_, 1.);
  Mv.Update(gamma(2), id2V_, 1.);
  for (int i = 0; i < 3; i++) MandelStressIsoprinc(i, i) += Mv(i);
  MandelStressIsoprinc(0, 1) += Mv(3);
  MandelStressIsoprinc(1, 0) += Mv(3);
  MandelStressIsoprinc(1, 2) += Mv(4);
  MandelStressIsoprinc(2, 1) += Mv(4);
  MandelStressIsoprinc(0, 2) += Mv(5);
  MandelStressIsoprinc(2, 0) += Mv(5);

  // derivative of Mandel stress w.r.t. GL
  AddSymmetricHolzapfelProduct(dMdCisoprinc, invpldefgrd_, invpldefgrd_, .5 * gamma(0));
  AddSymmetricHolzapfelProduct(dMdCisoprinc, invpldefgrd_, FpiCe_, gamma(1));
  dMdCisoprinc.MultiplyNT(delta(0), Ce_, Cpi_, 1.);
  dMdCisoprinc.MultiplyNT(delta(1), Ce_, CpiCCpi_, 1.);
  dMdCisoprinc.MultiplyNT(delta(1), Ce2_, Cpi_, 1.);
  dMdCisoprinc.MultiplyNT(delta(2), Ce_, ircg_, 1.);
  dMdCisoprinc.MultiplyNT(delta(2), id2V_, Cpi_, 1.);
  dMdCisoprinc.MultiplyNT(delta(3), Ce2_, CpiCCpi_, 1.);
  dMdCisoprinc.MultiplyNT(delta(4), Ce2_, ircg_, 1.);
  dMdCisoprinc.MultiplyNT(delta(4), id2V_, CpiCCpi_, 1.);
  dMdCisoprinc.MultiplyNT(delta(5), id2V_, ircg_, 1.);

  // derivative of Mandel stress w.r.t. inverse plastic deformation gradient
  AddRightNonSymmetricHolzapfelProduct(dMdFpinvIsoprinc, FpiTC_, id2_, gamma(0));
  AddRightNonSymmetricHolzapfelProduct(dMdFpinvIsoprinc, FpiTC_, CeM_, gamma(1));
  AddRightNonSymmetricHolzapfelProduct(dMdFpinvIsoprinc, CeFpiTC_, id2_, gamma(1));
  dMdFpinvIsoprinc.MultiplyNT(delta(0), Ce_, CFpi_, 1.);
  dMdFpinvIsoprinc.MultiplyNT(delta(1), Ce_, CFpiCe_, 1.);
  dMdFpinvIsoprinc.MultiplyNT(delta(1), Ce2_, CFpi_, 1.);
  dMdFpinvIsoprinc.MultiplyNT(delta(2), Ce_, CFpiCei_, 1.);
  dMdFpinvIsoprinc.MultiplyNT(delta(2), id2V_, CFpi_, 1.);
  dMdFpinvIsoprinc.MultiplyNT(delta(3), Ce2_, CFpiCe_, 1.);
  dMdFpinvIsoprinc.MultiplyNT(delta(4), Ce2_, CFpiCei_, 1.);
  dMdFpinvIsoprinc.MultiplyNT(delta(4), id2V_, CFpiCe_, 1.);
  dMdFpinvIsoprinc.MultiplyNT(delta(5), id2V_, CFpiCei_, 1.);

  return;
}

/*---------------------------------------------------------------------*
 | return names of visualization data (public)                         |
 *---------------------------------------------------------------------*/
void Mat::PlasticElastHyper::VisNames(std::map<std::string, int>& names) const
{
  std::string accumulatedstrain = "accumulatedstrain";
  names[accumulatedstrain] = 1;  // scalar
  std::string plastic_strain_incr = "plastic_strain_incr";
  names[plastic_strain_incr] = 1;  // scalar
  std::string plastic_zone = "plastic_zone";
  names[plastic_zone] = 1;  // scalar
  std::string kinematic_plastic_strain = "kinematic_plastic_strain";
  names[kinematic_plastic_strain] = 9;  // tensor

}  // VisNames()


/*---------------------------------------------------------------------*
 | return visualization data (public)                                  |
 *---------------------------------------------------------------------*/
bool Mat::PlasticElastHyper::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (name == "accumulatedstrain")
  {
    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");
    double tmp = 0.;
    for (unsigned gp = 0; gp < last_alpha_isotropic_.size(); gp++) tmp += AccumulatedStrain(gp);
    data[0] = tmp / last_alpha_isotropic_.size();
  }
  else if (name == "plastic_strain_incr")
  {
    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");
    double tmp = 0.;
    for (unsigned gp = 0; gp < delta_alpha_i_.size(); gp++) tmp += delta_alpha_i_.at(gp);
    data[0] = tmp / delta_alpha_i_.size();
  }
  else if (name == "plastic_zone")
  {
    bool plastic_history = false;
    bool curr_active = false;
    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");
    for (unsigned gp = 0; gp < last_alpha_isotropic_.size(); gp++)
    {
      if (AccumulatedStrain(gp) != 0.) plastic_history = true;
      if (Active(gp)) curr_active = true;
    }
    data[0] = plastic_history + curr_active;
  }
  else if (name == "kinematic_plastic_strain")
  {
    if ((int)data.size() != 9) FOUR_C_THROW("size mismatch");
    std::vector<double> tmp(9, 0.);
    for (std::size_t gp = 0; gp < last_alpha_kinematic_.size(); ++gp)
    {
      const double* values = last_alpha_kinematic_[gp].A();
      for (std::size_t i = 0; i < 9; ++i)
      {
        tmp[i] += values[i];
      }
    }
    for (std::size_t i = 0; i < 9; ++i) data[i] = tmp[i] / last_alpha_kinematic_.size();
  }
  else
  {
    return false;
  }
  return true;

}  // VisData()

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
void Mat::PlasticElastHyper::register_output_data_names(
    std::unordered_map<std::string, int>& names_and_size) const
{
  names_and_size["accumulated_plastic_strain"] = 1;  // scalar
  names_and_size["plastic_strain_incr"] = 1;         // scalar
  names_and_size["plastic_zone"] = 1;                // scalar
  names_and_size["kinematic_plastic_strain"] = 9;    // tensor
}

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
bool Mat::PlasticElastHyper::EvaluateOutputData(
    const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const
{
  if (name == "accumulated_plastic_strain")
  {
    for (std::size_t gp = 0; gp < last_alpha_isotropic_.size(); ++gp)
    {
      data(gp, 0) = last_alpha_isotropic_[gp];
    }
    return true;
  }
  if (name == "plastic_strain_incr")
  {
    for (std::size_t gp = 0; gp < delta_alpha_i_.size(); ++gp)
    {
      data(gp, 0) = delta_alpha_i_[gp];
    }
    return true;
  }
  if (name == "plastic_zone")
  {
    bool plastic_history = false;
    bool curr_active = false;
    for (std::size_t gp = 0; gp < last_alpha_isotropic_.size(); ++gp)
    {
      if (AccumulatedStrain(gp) != 0.) plastic_history = true;
      if (Active(gp)) curr_active = true;
      data(gp, 0) = plastic_history + curr_active;
    }
    return true;
  }
  if (name == "kinematic_plastic_strain")
  {
    for (std::size_t gp = 0; gp < last_alpha_kinematic_.size(); ++gp)
    {
      const double* values = last_alpha_kinematic_[gp].A();
      for (std::size_t i = 0; i < 9; ++i)
      {
        data(gp, i) = values[i];
      }
    }
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------*
 |  matrix exponential                                      seitz 07/13 |
 *----------------------------------------------------------------------*/
void Mat::PlasticElastHyper::matrix_exponential3x3(Core::LinAlg::Matrix<3, 3>& MatrixInOut)
{
  double Norm = MatrixInOut.Norm2();
  // direct calculation for zero-matrix
  if (Norm == 0.)
  {
    MatrixInOut.Clear();
    for (int i = 0; i < 3; i++) MatrixInOut(i, i) = 1.;
    return;
  }

  // Calculation of matrix exponential via power series. This is usually
  // faster than by polar decomposition for matrices are close to zero.
  // For small plastic increments this is the case
  Core::LinAlg::Matrix<3, 3> In(MatrixInOut);
  int n = 0;
  int facn = 1;
  MatrixInOut.Clear();
  for (int i = 0; i < 3; i++) MatrixInOut(i, i) = 1.;
  Core::LinAlg::Matrix<3, 3> tmp(MatrixInOut);
  Core::LinAlg::Matrix<3, 3> tmp2(MatrixInOut);
  while (n < 50 && tmp.Norm2() / facn > 1.e-32)
  {
    n++;
    facn *= n;
    tmp.Multiply(tmp2, In);
    tmp2 = tmp;
    MatrixInOut.Update(1. / facn, tmp, 1.);
  }
  if (n == 50) FOUR_C_THROW("matrix exponential unconverged in %i steps", n);

  return;
}

/*---------------------------------------------------------------------------*
 |  matrix exponential derivative of a symmetric matrix          seitz 07/13 |
 *---------------------------------------------------------------------------*/
void Mat::PlasticElastHyper::matrix_exponential_derivative_sym3x3(
    const Core::LinAlg::Matrix<3, 3> MatrixIn, Core::LinAlg::Matrix<6, 6>& MatrixExpDeriv)
{
  double norm = MatrixIn.Norm2();

  Core::LinAlg::Matrix<6, 6> id4sharp(true);
  for (int i = 0; i < 3; i++) id4sharp(i, i) = 1.0;
  for (int i = 3; i < 6; i++) id4sharp(i, i) = 0.5;

  // direct calculation for zero-matrix
  if (norm == 0.)
  {
    MatrixExpDeriv = id4sharp;
    return;
  }

  if (norm < 0.3)
  {
    // see Souza-Neto: Computational Methods for plasticity, Box B.2.
    int nmax = 0;
    int nIter = 0;
    int nfac = 1;
    Core::LinAlg::Matrix<3, 3> tmp1;
    Core::LinAlg::Matrix<3, 3> tmp2(true);
    for (int i = 0; i < 3; i++) tmp2(i, i) = 1.;

    // all needed powers of X
    std::vector<Core::LinAlg::Matrix<3, 3>> Xn;
    Xn.resize(0);
    Xn.push_back(tmp2);

    // all needed factorials
    std::vector<int> fac;
    fac.resize(0);
    fac.push_back(nfac);

    // compute nmax and Xn
    while (nIter < 50 && tmp2.Norm2() / nfac > 1.e-32)
    {
      nIter++;
      nfac *= nIter;
      fac.push_back(nfac);
      tmp1.Multiply(tmp2, MatrixIn);
      Xn.push_back(tmp1);
      tmp2 = tmp1;
    }
    if (nIter == 50) FOUR_C_THROW("matrix exponential unconverged in %i steps", nIter);
    nmax = nIter;

    // compose derivative of matrix exponential (symmetric Voigt-notation)
    MatrixExpDeriv.Clear();
    for (int n = 1; n <= nmax; n++)
    {
      for (int m = 1; m <= n / 2; m++)
        AddSymmetricHolzapfelProduct(MatrixExpDeriv, Xn.at(m - 1), Xn.at(n - m), .5 / fac[n]);
      if (n % 2 == 1)
        AddSymmetricHolzapfelProduct(
            MatrixExpDeriv, Xn.at((n - 1) / 2), Xn.at((n - 1) / 2), .25 / fac[n]);
    }
  }
  else
  {
    double EWtolerance = 1.e-12;

    Core::LinAlg::Matrix<3, 3> EV(MatrixIn);
    Core::LinAlg::Matrix<3, 3> EW;
    Core::LinAlg::SYEV(EV, EW, EV);

    Core::LinAlg::Matrix<3, 1> vec1;
    Core::LinAlg::Matrix<3, 1> vec2;
    Core::LinAlg::Matrix<3, 3> tmp1;
    Core::LinAlg::Matrix<3, 3> tmp2;

    MatrixExpDeriv.Clear();
    // souza eq. (A.52)
    // note: EW stored in ascending order

    //  d X^2 / d X  =  1/2 * (  delta_jk X_lj + delta_il X_kj
    //                         + delta_jl X_ik + delta_kj X_il )
    //
    // y_i = log(x_i)
    // dy_i / dx_j = delta_ij 1/x_i

    Core::LinAlg::Matrix<3, 3> id2(true);
    for (int i = 0; i < 3; i++) id2(i, i) = 1.0;
    //  // --------------------------------- switch by number of equal eigenvalues

    if (abs(EW(0, 0) - EW(1, 1)) < EWtolerance &&
        abs(EW(1, 1) - EW(2, 2)) < EWtolerance)  // ------------------ x_a == x_b == x_c
    {
      // calculate derivative
      MatrixExpDeriv = id4sharp;
      MatrixExpDeriv.Scale(exp(EW(0, 0)));
    }

    else if ((abs(EW(0, 0) - EW(1, 1)) < EWtolerance && abs(EW(1, 1) - EW(2, 2)) > EWtolerance) ||
             (abs(EW(0, 0) - EW(1, 1)) > EWtolerance &&
                 abs(EW(1, 1) - EW(2, 2)) <
                     EWtolerance))  // ---- x_a != x_b == x_c or x_a == x_b != x_c
    {
      // factors
      double s1 = 0.0;
      double s2 = 0.0;
      double s3 = 0.0;
      double s4 = 0.0;
      double s5 = 0.0;
      double s6 = 0.0;

      int a = 0;
      int c = 0;

      // switch which two EW are equal
      if (abs(EW(0, 0) - EW(1, 1)) < EWtolerance &&
          abs(EW(1, 1) - EW(2, 2)) > EWtolerance)  // ----------------------- x_a == x_b != x_c
      {
        a = 2;
        c = 0;
      }
      else if (abs(EW(0, 0) - EW(1, 1)) > EWtolerance &&
               abs(EW(1, 1) - EW(2, 2)) < EWtolerance)  // ------------------ x_a != x_b == x_c
      {
        a = 0;
        c = 2;
      }
      else
        FOUR_C_THROW("you should not be here");

      // in souza eq. (A.53):
      s1 = (exp(EW(a, a)) - exp(EW(c, c))) / (pow(EW(a, a) - EW(c, c), 2.0)) -
           exp(EW(c, c)) / (EW(a, a) - EW(c, c));
      s2 = 2.0 * EW(c, c) * (exp(EW(a, a)) - exp(EW(c, c))) / (pow(EW(a, a) - EW(c, c), 2.0)) -
           (EW(a, a) + EW(c, c)) / (EW(a, a) - EW(c, c)) * exp(EW(c, c));
      s3 = 2.0 * (exp(EW(a, a)) - exp(EW(c, c))) / (pow(EW(a, a) - EW(c, c), 3.0)) -
           (exp(EW(a, a)) + exp(EW(c, c))) / (pow(EW(a, a) - EW(c, c), 2.0));
      s4 = EW(c, c) * s3;
      s5 = s4;
      s6 = EW(c, c) * EW(c, c) * s3;

      // calculate derivative
      Mat::AddToCmatDerivTensorSquare(MatrixExpDeriv, s1, MatrixIn, 1.);
      MatrixExpDeriv.Update(-s2, id4sharp, 1.);
      Mat::ElastSymTensorMultiply(MatrixExpDeriv, -1. * s3, MatrixIn, MatrixIn, 1.);
      Mat::ElastSymTensorMultiply(MatrixExpDeriv, s4, MatrixIn, id2, 1.);
      Mat::ElastSymTensorMultiply(MatrixExpDeriv, s5, id2, MatrixIn, 1.);
      Mat::ElastSymTensorMultiply(MatrixExpDeriv, -s6, id2, id2, 1.);
    }

    else if (abs(EW(0, 0) - EW(1, 1)) > EWtolerance &&
             abs(EW(1, 1) - EW(2, 2)) > EWtolerance)  // ----------------- x_a != x_b != x_c
    {
      for (int a = 0; a < 3; a++)  // loop over all eigenvalues
      {
        int b = (a + 1) % 3;
        int c = (a + 2) % 3;

        Core::LinAlg::Matrix<3, 1> ea;
        Core::LinAlg::Matrix<3, 1> eb;
        Core::LinAlg::Matrix<3, 1> ec;
        for (int i = 0; i < 3; i++)
        {
          ea(i) = EV(i, a);
          eb(i) = EV(i, b);
          ec(i) = EV(i, c);
        }
        Core::LinAlg::Matrix<3, 3> Ea;
        Ea.MultiplyNT(ea, ea);
        Core::LinAlg::Matrix<3, 3> Eb;
        Eb.MultiplyNT(eb, eb);
        Core::LinAlg::Matrix<3, 3> Ec;
        Ec.MultiplyNT(ec, ec);

        double fac = exp(EW(a, a)) / ((EW(a, a) - EW(b, b)) * (EW(a, a) - EW(c, c)));

        // + d X^2 / d X
        Mat::AddToCmatDerivTensorSquare(MatrixExpDeriv, fac, MatrixIn, 1.);

        // - (x_b + x_c) I_s
        MatrixExpDeriv.Update(-1. * (EW(b, b) + EW(c, c)) * fac, id4sharp, 1.);

        // - [(x_a - x_b) + (x_a - x_c)] E_a \dyad E_a
        Mat::ElastSymTensorMultiply(MatrixExpDeriv,
            -1. * fac * ((EW(a, a) - EW(b, b)) + (EW(a, a) - EW(c, c))), Ea, Ea, 1.);


        // - (x_b - x_c) (E_b \dyad E_b)
        Mat::ElastSymTensorMultiply(MatrixExpDeriv, -1. * fac * (EW(b, b) - EW(c, c)), Eb, Eb, 1.);

        // + (x_b - x_c) (E_c \dyad E_c)
        Mat::ElastSymTensorMultiply(MatrixExpDeriv, fac * (EW(b, b) - EW(c, c)), Ec, Ec, 1.);

        // dy / dx_a E_a \dyad E_a
        Mat::ElastSymTensorMultiply(MatrixExpDeriv, exp(EW(a, a)), Ea, Ea, 1.);
      }  // end loop over all eigenvalues
    }

    else
      FOUR_C_THROW("you should not be here.");
  }
  return;
}

/*---------------------------------------------------------------------------*
 |  matrix exponential derivative of a symmetric matrix          seitz 09/13 |
 *---------------------------------------------------------------------------*/
void Mat::PlasticElastHyper::matrix_exponential_derivative3x3(
    const Core::LinAlg::Matrix<3, 3> MatrixIn, Core::LinAlg::Matrix<9, 9>& MatrixExpDeriv)
{
  // see Souza-Neto: Computational Methods for plasticity, Box B.2.
  int nmax = 0;
  int nIter = 0;
  int nfac = 1;
  Core::LinAlg::Matrix<3, 3> tmp1;
  Core::LinAlg::Matrix<3, 3> tmp2(true);
  for (int i = 0; i < 3; i++) tmp2(i, i) = 1.;

  // all needed powers of X
  std::vector<Core::LinAlg::Matrix<3, 3>> Xn;
  Xn.resize(0);
  Xn.push_back(tmp2);

  // all needed factorials
  std::vector<int> fac;
  fac.resize(0);
  fac.push_back(nfac);

  // compute nmax and Xn
  while (nIter < 50 && tmp2.Norm2() / nfac > 1.e-32)
  {
    nIter++;
    nfac *= nIter;
    fac.push_back(nfac);
    tmp1.Multiply(tmp2, MatrixIn);
    Xn.push_back(tmp1);
    tmp2 = tmp1;
  }
  if (nIter == 50) FOUR_C_THROW("matrix exponential unconverged in %i steps", nIter);
  nmax = nIter;

  // compose derivative of matrix exponential (non-symmetric Voigt-notation)
  MatrixExpDeriv.Clear();
  for (int n = 1; n <= nmax; n++)
    for (int m = 1; m <= n; m++)
      AddNonSymmetricProduct(1. / fac[n], Xn.at(m - 1), Xn.at(n - m), MatrixExpDeriv);

  return;
}

FOUR_C_NAMESPACE_CLOSE
