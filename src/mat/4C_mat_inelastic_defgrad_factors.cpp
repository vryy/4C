/*----------------------------------------------------------------------*/
/*! \file
\brief evaluation of inelastic deformation gradients and their derivatives

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_inelastic_defgrad_factors.hpp"

#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_electrode.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_utils_function_of_time.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticDefgradNoGrowth::InelasticDefgradNoGrowth(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata)
{
  // do nothing here
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticDefgradScalar::InelasticDefgradScalar(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      scalar1_(matdata.parameters.get<int>("SCALAR1")),
      scalar1_ref_conc_(matdata.parameters.get<double>("SCALAR1_RefConc"))
{
  // safety checks
  // in case not all scatra dofs are transported scalars, the last scatra dof is a potential and can
  // not be treated as a concentration but it is treated like that in so3_scatra_evaluate.cpp in the
  // pre_evaluate method!
  if (scalar1_ != 1) FOUR_C_THROW("At the moment it is only possible that SCALAR1 induces growth");
  if (matdata.parameters.get<double>("SCALAR1_RefConc") < 0.0)
    FOUR_C_THROW("The reference concentration of SCALAR1 can't be negative");
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticDefgradLinScalar::InelasticDefgradLinScalar(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : InelasticDefgradScalar(matdata),
      scalar1_molar_growth_fac_(matdata.parameters.get<double>("SCALAR1_MolarGrowthFac"))
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticDefgradIntercalFrac::InelasticDefgradIntercalFrac(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : InelasticDefgradScalar(matdata)
{
  // get matid
  const int matid = matdata.parameters.get<int>("MATID");

  // Check if the material specified by user with MATID is an electrode material
  if (matid > 0)
  {
    // retrieve problem instance to read from
    const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
    // retrieve validated input line of material ID in question
    auto* curmat = Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    switch (curmat->Type())
    {
      case Core::Materials::m_electrode:
      {
        // Get C_max and Chi_max of electrode material
        c_max_ = curmat->raw_parameters().get<double>("C_MAX");
        chi_max_ = curmat->raw_parameters().get<double>("CHI_MAX");
        break;
      }
      default:
        FOUR_C_THROW("The material you have specified by MATID has to be an electrode material!");
    }
  }
  else
  {
    FOUR_C_THROW("You have to enter a valid MATID for the corresponding electrode material!");
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticDefgradPolyIntercalFrac::InelasticDefgradPolyIntercalFrac(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : InelasticDefgradIntercalFrac(matdata),
      poly_coeffs_(matdata.parameters.get<std::vector<double>>("POLY_PARAMS")),
      x_max_(matdata.parameters.get<double>("X_max")),
      x_min_(matdata.parameters.get<double>("X_min"))
{
  // safety check
  if (poly_coeffs_.size() !=
      static_cast<unsigned int>(matdata.parameters.get<int>("POLY_PARA_NUM")))
  {
    FOUR_C_THROW(
        "Number of coefficients POLY_PARA_NUM you entered in input file has to match the size "
        "of coefficient vector POLY_PARAMS");
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticDefgradLinScalarAniso::InelasticDefgradLinScalarAniso(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : InelasticDefgradLinScalar(matdata),
      growth_dir_(Teuchos::rcp(new InelasticDeformationDirection(
          matdata.parameters.get<std::vector<double>>("GrowthDirection"))))
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticDefgradPolyIntercalFracAniso::InelasticDefgradPolyIntercalFracAniso(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : InelasticDefgradPolyIntercalFrac(matdata),
      growth_dir_(Teuchos::rcp(new InelasticDeformationDirection(
          matdata.parameters.get<std::vector<double>>("GrowthDirection"))))
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticDeformationDirection::InelasticDeformationDirection(
    std::vector<double> growthdirection)
    : growth_dir_mat_(true)
{
  if (growthdirection.size() != 3)
  {
    FOUR_C_THROW(
        "Since we have a 3D problem here, vector that defines the growth direction also needs to "
        "have the size 3!");
  }

  // fill matrix that determines the growth direction
  const double growthdirvecnorm =
      std::sqrt(std::pow(growthdirection[0], 2.0) + std::pow(growthdirection[1], 2.0) +
                std::pow(growthdirection[2], 2.0));
  const double invquadrgrowthdirvecnorm = 1.0 / (growthdirvecnorm * growthdirvecnorm);

  // loop over all rows and columns to fill the matrix and scale it correctly on the fly
  for (unsigned i = 0; i < growthdirection.size(); ++i)
  {
    for (unsigned j = 0; j < growthdirection.size(); ++j)
    {
      growth_dir_mat_(i, j) = invquadrgrowthdirvecnorm * growthdirection[i] * growthdirection[j];
    }
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticDefgradLinTempIso::InelasticDefgradLinTempIso(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      ref_temp_(matdata.parameters.get<double>("RefTemp")),
      temp_growth_fac_(matdata.parameters.get<double>("Temp_GrowthFac"))

{
  // safety checks
  if (ref_temp_ < 0.0) FOUR_C_THROW("Avoid negative reference temperatures");
  if (temp_growth_fac_ == 0.0)
  {
    FOUR_C_THROW(
        "Do not use 'MAT_InelasticDefgradLinTempIso' with a growth factor of 0.0. Use "
        "'MAT_InelasticDefgradNoGrowth' instead!");
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticDefgradTimeFunct::InelasticDefgradTimeFunct(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata), funct_num_(matdata.parameters.get<int>("FUNCT_NUM"))
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::InelasticDefgradFactors::InelasticDefgradFactors(Core::Mat::PAR::Parameter* params)
    : gp_(-1), params_(params)
{
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Teuchos::RCP<Mat::InelasticDefgradFactors> Mat::InelasticDefgradFactors::Factory(int matnum)
{
  // for the sake of safety
  if (Global::Problem::Instance()->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");

  // another safety check
  if (Global::Problem::Instance()->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // check correct masslin type
  const Teuchos::ParameterList& sdyn = Global::Problem::Instance()->structural_dynamic_params();
  if (Core::UTILS::IntegralValue<Inpar::Solid::MassLin>(sdyn, "MASSLIN") != Inpar::Solid::ml_none)
  {
    FOUR_C_THROW(
        "If you use the material 'InelasticDefgradFactors' please set 'MASSLIN' in the "
        "STRUCTURAL DYNAMIC Section to 'None', or feel free to implement other possibility!");
  }

  // retrieve problem instance to read from
  const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  auto* curmat = Global::Problem::Instance(probinst)->Materials()->ParameterById(matnum);

  // get material type and call corresponding constructors
  const Core::Materials::MaterialType currentMaterialType = curmat->Type();
  switch (currentMaterialType)
  {
    case Core::Materials::mfi_no_growth:
    {
      auto* params = dynamic_cast<Mat::PAR::InelasticDefgradNoGrowth*>(curmat);

      return Teuchos::rcp(new InelasticDefgradNoGrowth(params));
    }
    case Core::Materials::mfi_lin_scalar_aniso:
    {
      // get pointer to parameter class
      auto* params = dynamic_cast<Mat::PAR::InelasticDefgradLinScalarAniso*>(curmat);

      // return pointer to inelastic deformation gradient object
      return Teuchos::rcp(new InelasticDefgradLinScalarAniso(params));
    }
    case Core::Materials::mfi_lin_scalar_iso:
    {
      // get pointer to parameter class
      auto* params = dynamic_cast<Mat::PAR::InelasticDefgradScalar*>(curmat);

      // return pointer to inelastic deformation gradient object
      return Teuchos::rcp(new InelasticDefgradLinScalarIso(params));
    }
    case Core::Materials::mfi_poly_intercal_frac_aniso:
    {
      // get pointer to parameter class
      auto* params = dynamic_cast<Mat::PAR::InelasticDefgradPolyIntercalFracAniso*>(curmat);

      // return pointer to inelastic deformation gradient object
      return Teuchos::rcp(new InelasticDefgradPolyIntercalFracAniso(params));
    }
    case Core::Materials::mfi_poly_intercal_frac_iso:
    {
      // get pointer to parameter class
      auto* params = dynamic_cast<Mat::PAR::InelasticDefgradPolyIntercalFrac*>(curmat);

      // return pointer to inelastic deformation gradient object
      return Teuchos::rcp(new InelasticDefgradPolyIntercalFracIso(params));
    }

    case Core::Materials::mfi_lin_temp_iso:
    {
      auto* params = dynamic_cast<Mat::PAR::InelasticDefgradLinTempIso*>(curmat);
      return Teuchos::rcp(new InelasticDefgradLinTempIso(params));
    }
    case Core::Materials::mfi_time_funct:
    {
      auto* params = dynamic_cast<Mat::PAR::InelasticDefgradTimeFunct*>(curmat);
      return Teuchos::rcp(new InelasticDefgradTimeFunct(params));
    }

    default:
      FOUR_C_THROW("cannot deal with type %d", curmat->Type());
  }
  // dummy return
  return Teuchos::null;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::InelasticDefgradScalar::InelasticDefgradScalar(Core::Mat::PAR::Parameter* params)
    : InelasticDefgradFactors(params), concentrations_(Teuchos::null)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradScalar::pre_evaluate(Teuchos::ParameterList& params, const int gp)
{
  set_gp(gp);

  // set pointer to vector of gp_conc, only if gp is 0, because this is the first gp
  if (gp == 0)
    concentrations_ = params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc");
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradScalar::SetConcentrationGP(const double concentration)
{
  // this method is only called for a certain gauss point whose id is not accessible, thus we set a
  // dummy id here and set the corresponding concentration value afterwards
  const int dummy_gp(0);
  set_gp(dummy_gp);

  const int scalar1 = Parameter()->Scalar1();

  get_concentration_gp().at(scalar1 - 1) = concentration;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::InelasticDefgradPolyIntercalFrac::InelasticDefgradPolyIntercalFrac(
    Core::Mat::PAR::Parameter* params)
    : InelasticDefgradScalar(params)
{
  polynomial_growth_ = Teuchos::rcp(new InelasticDefgradPolynomialShape(
      Parameter()->PolyCoeffs(), Parameter()->XMin(), Parameter()->XMax()));

  // get reference intercalation fraction
  const double x_ref = Mat::Electrode::compute_intercalation_fraction(
      Parameter()->Scalar1RefConc(), Parameter()->Chimax(), Parameter()->Cmax(), 1.0);

  // set the polynomial value in the reference configuration
  Parameter()->set_polynom_reference_value(polynomial_growth_->ComputePolynomial(x_ref));
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double Mat::InelasticDefgradPolyIntercalFrac::EvaluatePolynomial(
    const double concentration, const double detjacobian)
{
  // get intercalation fraction
  const double x = Mat::Electrode::compute_intercalation_fraction(
      concentration, Parameter()->Chimax(), Parameter()->Cmax(), detjacobian);

  // check bounds of validity of polynomial
  polynomial_growth_->check_polynomial_bounds(x);

  // calculate and return the value of the polynomial
  return polynomial_growth_->ComputePolynomial(x);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double Mat::InelasticDefgradPolyIntercalFrac::evaluate_polynomial_derivative(
    const double concentration, const double detjacobian)
{
  // get intercalation fraction
  const double x = Mat::Electrode::compute_intercalation_fraction(
      concentration, Parameter()->Chimax(), Parameter()->Cmax(), detjacobian);

  return polynomial_growth_->compute_polynomial_derivative(x);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticSource Mat::InelasticDefgradPolyIntercalFrac::GetInelasticSource()
{
  return Mat::PAR::InelasticSource::concentration;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::InelasticDefgradLinScalarIso::InelasticDefgradLinScalarIso(Core::Mat::PAR::Parameter* params)
    : InelasticDefgradScalar(params)
{
  linear_growth_ = Teuchos::rcp(new InelasticDefgradLinearShape(
      Parameter()->scalar1_molar_growth_fac(), Parameter()->Scalar1RefConc()));
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticSource Mat::InelasticDefgradLinScalarIso::GetInelasticSource()
{
  return Mat::PAR::InelasticSource::concentration;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinScalarIso::evaluate_inverse_inelastic_def_grad(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, Core::LinAlg::Matrix<3, 3>& iFinM)
{
  // get parameter
  const int sc1 = Parameter()->Scalar1();
  const double material_concentration = get_concentration_gp().at(sc1 - 1) * defgrad->determinant();

  // get growth factor
  const double growth_factor = linear_growth_->evaluate_linear_growth(material_concentration);

  const double isoinelasticdefo = std::pow(1.0 + growth_factor, (1.0 / 3.0));

  // calculate inverse inelastic deformation gradient
  for (int i = 0; i < 3; ++i) iFinM(i, i) = 1.0 / isoinelasticdefo;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinScalarIso::evaluate_additional_cmat(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<3, 3>& iFinjM,
    const Core::LinAlg::Matrix<6, 1>& iCV, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
    Core::LinAlg::Matrix<6, 6>& cmatadd)
{
  // static variables
  static Core::LinAlg::Matrix<9, 6> diFinjdC(true);
  static Core::LinAlg::Matrix<9, 1> id9x1(true);

  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double sc1GrowthFac = linear_growth_->GrowthFac();
  const double concentration = get_concentration_gp().at(sc1 - 1);
  const double detjacobian = defgrad->determinant();

  // get growth factor
  const double growth_factor = linear_growth_->evaluate_linear_growth(concentration * detjacobian);

  // evaluate scaling factor
  const double scalefac =
      -sc1GrowthFac * concentration * detjacobian / 6.0 * std::pow(1 + growth_factor, -4.0 / 3.0);

  // calculate diFindC
  diFinjdC.multiply_nt(scalefac, id9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.multiply_nn(2.0, dSdiFinj, diFinjdC, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinScalarIso::EvaluateODStiffMat(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<3, 3>& iFinjM,
    const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 1>& dstressdc)
{
  static Core::LinAlg::Matrix<9, 1> id9x1(true);
  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double sc1GrowthFac = linear_growth_->GrowthFac();
  const double detjacobian = defgrad->determinant();
  const double material_concentration = get_concentration_gp().at(sc1 - 1) * detjacobian;

  // get growth factor
  const double growth_factor = linear_growth_->evaluate_linear_growth(material_concentration);

  // calculate scalefac
  const double scalefac =
      -sc1GrowthFac / 3.0 * detjacobian * std::pow(1 + growth_factor, -4.0 / 3.0);

  // calculate diFindc and add contribution to dstressdc = dSdiFinj : diFinjdc
  dstressdc.multiply_nn(scalefac, dSdiFinj, id9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinScalarIso::evaluate_inelastic_def_grad_derivative(
    const double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx)
{
  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double material_concentration = get_concentration_gp().at(sc1 - 1) * detjacobian;

  // get growth factor
  const double growth_factor = linear_growth_->evaluate_linear_growth(material_concentration);
  // calculate the scale factor needed to calculate the derivative below
  const double scalefac = 1.0 / 3.0 * std::pow(1 + growth_factor, -2.0 / 3.0) *
                          linear_growth_->GrowthFac() * detjacobian;

  // prepare identity tensor as 9x1 vector
  static Core::LinAlg::Matrix<9, 1> id9x1(true);
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // here dFindc is zeroed out and filled with the current value
  dFindx.update(scalefac, id9x1, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::InelasticDefgradLinScalarAniso::InelasticDefgradLinScalarAniso(
    Core::Mat::PAR::Parameter* params)
    : InelasticDefgradScalar(params)
{
  linear_growth_ = Teuchos::rcp(new InelasticDefgradLinearShape(
      Parameter()->scalar1_molar_growth_fac(), Parameter()->Scalar1RefConc()));
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticSource Mat::InelasticDefgradLinScalarAniso::GetInelasticSource()
{
  return Mat::PAR::InelasticSource::concentration;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinScalarAniso::evaluate_inverse_inelastic_def_grad(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, Core::LinAlg::Matrix<3, 3>& iFinM)
{
  // init and clear variable
  static Core::LinAlg::Matrix<3, 3> FinM(true);
  FinM.clear();

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double material_concentration = get_concentration_gp().at(sc1 - 1) * defgrad->determinant();

  // get growth factor
  const double growth_factor = linear_growth_->evaluate_linear_growth(material_concentration);

  // calculate inelastic deformation gradient and its inverse
  for (int i = 0; i < 3; ++i) FinM(i, i) = 1.0;

  // finalize inelastic deformation gradient matrix (FinM is calculated, such that the volume change
  // is a linear function of the scalar (mapped to reference frame) that causes it)
  FinM.update(growth_factor, Parameter()->GrowthDirMat(), 1.0);

  // calculate inverse of inelastic deformation gradient matrix
  iFinM.invert(FinM);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinScalarAniso::evaluate_additional_cmat(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<3, 3>& iFinjM,
    const Core::LinAlg::Matrix<6, 1>& iCV, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
    Core::LinAlg::Matrix<6, 6>& cmatadd)
{
  static Core::LinAlg::Matrix<3, 3> temp(true);
  static Core::LinAlg::Matrix<3, 3> iFinjGiFinj(true);
  static Core::LinAlg::Matrix<9, 1> iFinjGiFinj9x1(true);
  static Core::LinAlg::Matrix<9, 6> diFinjdC(true);

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double sc1GrowthFac = linear_growth_->GrowthFac();
  const double concentration = get_concentration_gp().at(sc1 - 1);
  const double detjacobian = defgrad->determinant();

  // prepare scalefac
  const double scalefac = -sc1GrowthFac * concentration * detjacobian / 2.0;

  // calculate F_{in,j}^{-1} . G . F_{in,j}^{-1} with F_{in,j}, the j-th factor of F_{in}
  temp.multiply_nn(1.0, iFinjM, Parameter()->GrowthDirMat(), 0.0);
  iFinjGiFinj.multiply_nn(1.0, temp, iFinjM, 0.0);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(iFinjGiFinj, iFinjGiFinj9x1);

  // calculate diFinjdC
  diFinjdC.multiply_nt(scalefac, iFinjGiFinj9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.multiply_nn(2.0, dSdiFinj, diFinjdC, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinScalarAniso::EvaluateODStiffMat(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<3, 3>& iFinjM,
    const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 1>& dstressdc)
{
  // static variables
  static Core::LinAlg::Matrix<3, 3> tmp(true);
  static Core::LinAlg::Matrix<3, 3> diFinjdcM(true);
  static Core::LinAlg::Matrix<9, 1> diFinjdc9x1(true);

  // get parameters
  const double sc1GrowthFac = linear_growth_->GrowthFac();
  const double detjacobian = defgrad->determinant();

  // prepare scalefac
  const double scalefac = -sc1GrowthFac * detjacobian;

  // calculate diFinjdc
  tmp.multiply_nn(1.0, iFinjM, Parameter()->GrowthDirMat(), 0.0);
  diFinjdcM.multiply_nn(scalefac, tmp, iFinjM, 0.0);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(diFinjdcM, diFinjdc9x1);

  // dstressdc = dSdiFinj : diFinjdc
  dstressdc.multiply_nn(1.0, dSdiFinj, diFinjdc9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinScalarAniso::evaluate_inelastic_def_grad_derivative(
    const double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx)
{
  const double scalefac = linear_growth_->GrowthFac() * detjacobian;

  // get the growth direction matrix as a 9x1 vector
  static Core::LinAlg::Matrix<9, 1> growthdirmat9x1(true);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(Parameter()->GrowthDirMat(), growthdirmat9x1);

  // here dFindc is zeroed out and filled with the current value
  dFindx.update(scalefac, growthdirmat9x1, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::InelasticDefgradPolyIntercalFracIso::InelasticDefgradPolyIntercalFracIso(
    Core::Mat::PAR::Parameter* params)
    : InelasticDefgradPolyIntercalFrac(params)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradPolyIntercalFracIso::evaluate_inverse_inelastic_def_grad(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, Core::LinAlg::Matrix<3, 3>& iFinM)
{
  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double polynomReferenceValue = Parameter()->get_polynom_reference_value();

  // get polynomial
  const double polynomValue =
      EvaluatePolynomial(get_concentration_gp().at(sc1 - 1), defgrad->determinant());

  // calculate growth
  const double isoInelasticDefo =
      std::pow((1.0 + polynomValue) / (1.0 + polynomReferenceValue), (1.0 / 3.0));
  // calculate inverse inelastic deformation gradient
  for (int i = 0; i < 3; ++i) iFinM(i, i) = 1.0 / isoInelasticDefo;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradPolyIntercalFracIso::evaluate_additional_cmat(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<3, 3>& iFinjM,
    const Core::LinAlg::Matrix<6, 1>& iCV, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
    Core::LinAlg::Matrix<6, 6>& cmatadd)
{
  // static variables
  static Core::LinAlg::Matrix<9, 6> diFinjdC(true);
  static Core::LinAlg::Matrix<9, 1> id9x1(true);

  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double chi_max = Parameter()->Chimax();
  const double c_max = Parameter()->Cmax();
  const double detjacobian = defgrad->determinant();
  const double concentration = get_concentration_gp().at(sc1 - 1);
  const double polynomReferenceValue = Parameter()->get_polynom_reference_value();

  // get polynomials
  const double polynomValue = EvaluatePolynomial(concentration, detjacobian);
  const double polynomDerivativeValue = evaluate_polynomial_derivative(concentration, detjacobian);

  // prepare scalefac
  const double scalefac = -1.0 / (6.0 * c_max) * concentration * chi_max * detjacobian *
                          std::pow(1.0 + polynomValue, -4.0 / 3.0) * polynomDerivativeValue *
                          std::pow(1.0 + polynomReferenceValue, 1.0 / 3.0);

  // calculate diFinjdC
  diFinjdC.multiply_nt(scalefac, id9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.multiply_nn(2.0, dSdiFinj, diFinjdC, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradPolyIntercalFracIso::EvaluateODStiffMat(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<3, 3>& iFinjM,
    const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 1>& dstressdc)
{
  static Core::LinAlg::Matrix<9, 1> id9x1(true);
  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double concentration = get_concentration_gp().at(sc1 - 1);
  const double detjacobian = defgrad->determinant();
  const double polynomReferenceValue = Parameter()->get_polynom_reference_value();

  // get polynomial and derivatives
  const double polynomValue = EvaluatePolynomial(concentration, detjacobian);
  const double polynomDerivativeValue = evaluate_polynomial_derivative(concentration, detjacobian);
  const double dChidc = Mat::Electrode::compute_d_intercalation_fraction_d_concentration(
      Parameter()->Chimax(), Parameter()->Cmax(), detjacobian);

  // prepare scalefac
  const double scalefac = -1.0 / 3.0 * std::pow(1.0 + polynomValue, -4.0 / 3.0) *
                          std::pow(1.0 + polynomReferenceValue, 1.0 / 3.0) *
                          polynomDerivativeValue * dChidc;

  // calculate diFinjdc and add contribution to dstressdc = dSdiFinj : diFinjdc
  dstressdc.multiply_nn(scalefac, dSdiFinj, id9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradPolyIntercalFracIso::evaluate_inelastic_def_grad_derivative(
    const double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx)
{
  static Core::LinAlg::Matrix<9, 1> id9x1(true);
  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double concentration = get_concentration_gp().at(sc1 - 1);
  const double polynomReferenceValue = Parameter()->get_polynom_reference_value();

  // get polynomial and its derivative
  const double polynomValue = EvaluatePolynomial(concentration, detjacobian);
  const double polynomDerivativeValue = evaluate_polynomial_derivative(concentration, detjacobian);

  // calculate the scale factor needed to get the derivative later
  const double denominator = 1.0 / (polynomReferenceValue + 1.0);
  const double base = (polynomValue + 1.0) * denominator;
  const double dChidc = Mat::Electrode::compute_d_intercalation_fraction_d_concentration(
      Parameter()->Chimax(), Parameter()->Cmax(), detjacobian);
  const double scalefac =
      1.0 / 3.0 * std::pow(base, -2.0 / 3.0) * polynomDerivativeValue * denominator * dChidc;

  // here dFindc is zeroed out and filled with the current value
  dFindx.update(scalefac, id9x1, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::InelasticDefgradPolyIntercalFracAniso::InelasticDefgradPolyIntercalFracAniso(
    Core::Mat::PAR::Parameter* params)
    : InelasticDefgradPolyIntercalFrac(params)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradPolyIntercalFracAniso::evaluate_inverse_inelastic_def_grad(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, Core::LinAlg::Matrix<3, 3>& iFinM)
{
  // init and clear variable
  static Core::LinAlg::Matrix<3, 3> FinM(true);
  FinM.clear();

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double polynomReferenceValue = Parameter()->get_polynom_reference_value();

  // get polynomials
  const double polynomValue =
      EvaluatePolynomial(get_concentration_gp().at(sc1 - 1), defgrad->determinant());

  // calculate growth factor
  const double growth_factor =
      (polynomValue - polynomReferenceValue) / (polynomReferenceValue + 1.0);

  // calculate inelastic deformation gradient and its inverse
  for (int i = 0; i < 3; ++i) FinM(i, i) = 1.0;

  // add the growth part
  FinM.update(growth_factor, Parameter()->GrowthDirMat(), 1.0);

  // calculate inverse of inelastic deformation gradient matrix
  iFinM.invert(FinM);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradPolyIntercalFracAniso::evaluate_additional_cmat(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<3, 3>& iFinjM,
    const Core::LinAlg::Matrix<6, 1>& iCV, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
    Core::LinAlg::Matrix<6, 6>& cmatadd)
{
  static Core::LinAlg::Matrix<3, 3> temp(true);
  static Core::LinAlg::Matrix<3, 3> iFinjGiFinj(true);
  static Core::LinAlg::Matrix<9, 1> iFinjGiFinj9x1(true);
  static Core::LinAlg::Matrix<9, 6> diFinjdC(true);

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double chi_max = Parameter()->Chimax();
  const double c_max = Parameter()->Cmax();
  const double concentration = get_concentration_gp().at(sc1 - 1);
  const double detjacobian = defgrad->determinant();
  const double polynomReferenceValue = Parameter()->get_polynom_reference_value();

  // get first derivative of polynomial
  const double polynomDerivativeValue = evaluate_polynomial_derivative(concentration, detjacobian);

  // prepare scalefac
  const double scalefac = -detjacobian * concentration * chi_max * polynomDerivativeValue /
                          (2.0 * c_max * (polynomReferenceValue + 1.0));

  // calculate F_{in,j}^{-1} . G . F_{in,j}^{-1} with F_{in,j}, the j-th factor of F_{in}
  temp.multiply_nn(1.0, iFinjM, Parameter()->GrowthDirMat(), 0.0);
  iFinjGiFinj.multiply_nn(1.0, temp, iFinjM, 0.0);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(iFinjGiFinj, iFinjGiFinj9x1);

  // calculate diFinjdC
  diFinjdC.multiply_nt(scalefac, iFinjGiFinj9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.multiply_nn(2.0, dSdiFinj, diFinjdC, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradPolyIntercalFracAniso::EvaluateODStiffMat(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<3, 3>& iFinjM,
    const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 1>& dstressdc)
{
  // static variables
  static Core::LinAlg::Matrix<3, 3> tmp(true);
  static Core::LinAlg::Matrix<3, 3> diFinjdcM(true);
  static Core::LinAlg::Matrix<9, 1> diFinjdc9x1(true);

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double detjacobian = defgrad->determinant();
  const double polynomReferenceValue = Parameter()->get_polynom_reference_value();

  // get derivatives
  const double polynomDerivativeValue =
      evaluate_polynomial_derivative(get_concentration_gp().at(sc1 - 1), detjacobian);
  const double dChidc = Mat::Electrode::compute_d_intercalation_fraction_d_concentration(
      Parameter()->Chimax(), Parameter()->Cmax(), detjacobian);

  // prepare scalefac
  const double scalefac = -polynomDerivativeValue / (polynomReferenceValue + 1.0) * dChidc;

  // calculate diFinjdc
  tmp.multiply_nn(1.0, iFinjM, Parameter()->GrowthDirMat(), 0.0);
  diFinjdcM.multiply_nn(scalefac, tmp, iFinjM, 0.0);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(diFinjdcM, diFinjdc9x1);

  // dstressdc = dSdiFinj : diFinjdc
  dstressdc.multiply_nn(1.0, dSdiFinj, diFinjdc9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradPolyIntercalFracAniso::evaluate_inelastic_def_grad_derivative(
    const double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx)
{
  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double concentration = get_concentration_gp().at(sc1 - 1);
  const double polynomReferenceValue = Parameter()->get_polynom_reference_value();

  // get polynomial derivative
  const double polynomDerivativeValue = evaluate_polynomial_derivative(concentration, detjacobian);

  const double dChidc = Mat::Electrode::compute_d_intercalation_fraction_d_concentration(
      Parameter()->Chimax(), Parameter()->Cmax(), detjacobian);
  const double scalefac = polynomDerivativeValue / (polynomReferenceValue + 1.0) * dChidc;

  // get the growth direction matrix as a 9x1 vector
  static Core::LinAlg::Matrix<9, 1> growthdirmat9x1(true);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(Parameter()->GrowthDirMat(), growthdirmat9x1);

  // here dFindc is zeroed out and filled with the current value
  dFindx.update(scalefac, growthdirmat9x1, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::InelasticDefgradLinearShape::InelasticDefgradLinearShape(
    const double growth_fac, const double reference_value)
    : growth_fac_(growth_fac), reference_value_(reference_value)
{
  // safety checks
  if (growth_fac < 0.0)
    FOUR_C_THROW("Growth factor can not be negative, please check your input file!");
  if (growth_fac == 0.0)
  {
    FOUR_C_THROW(
        "Do not use linear growth laws with a growth factor of 0.0. Use "
        "'MAT_InelasticDefgradNoGrowth' instead!");
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double Mat::InelasticDefgradLinearShape::evaluate_linear_growth(const double value) const
{
  // calculate and return the linear growth factor
  return growth_fac_ * (value - reference_value_);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::InelasticDefgradPolynomialShape::InelasticDefgradPolynomialShape(
    std::vector<double> poly_coeffs, const double x_min, const double x_max)
    : poly_coeffs_(std::move(poly_coeffs)), x_min_(x_min), x_max_(x_max)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double Mat::InelasticDefgradPolynomialShape::ComputePolynomial(const double x)
{
  // initialize the variable for the evaluation of the polynomial
  double polynom(0.0);

  // compute polynomial
  for (unsigned i = 0; i < poly_coeffs_.size(); ++i) polynom += poly_coeffs_[i] * std::pow(x, i);

  return polynom;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double Mat::InelasticDefgradPolynomialShape::compute_polynomial_derivative(const double x)
{
  // initialize the variable for the derivative of the polynomial
  double polynomDerivative(0.0);

  // compute first derivative of polynomial
  for (unsigned i = 1; i < poly_coeffs_.size(); ++i)
    polynomDerivative += i * poly_coeffs_[i] * std::pow(x, i - 1);

  return polynomDerivative;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradPolynomialShape::check_polynomial_bounds(const double x) const
{
  // safety check for validity of polynomial
  if ((x < x_min_) or (x > x_max_))
  {
    std::cout << "WARNING: Polynomial is evaluated outside its range of validity!" << std::endl;
    std::cout << "Evaluation at: " << x << " Lower bound is " << x_min_ << " Upper bound is "
              << x_max_ << std::endl;
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::InelasticDefgradLinTempIso::InelasticDefgradLinTempIso(Core::Mat::PAR::Parameter* params)
    : InelasticDefgradFactors(params), temperatures_(Teuchos::null)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinTempIso::pre_evaluate(Teuchos::ParameterList& params, int gp)
{
  // get Gauss point number
  set_gp(gp);

  // set pointer to vector of gp_temp, only if gp is 0, because this is the first gp
  if (gp == 0) temperatures_ = params.get<Teuchos::RCP<std::vector<double>>>("gp_temp");
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinTempIso::evaluate_inverse_inelastic_def_grad(
    const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM)
{
  // get parameters
  const double tempgrowthfac = Parameter()->GetTempGrowthFac();
  const double reftemp = Parameter()->RefTemp();

  const double growthfactor = 1.0 + tempgrowthfac * (get_temperature_gp() - reftemp);
  if (growthfactor <= 0.0) FOUR_C_THROW("Determinante of growth must not become negative");
  const double isoinelasticdefo = std::pow(growthfactor, (1.0 / 3.0));

  // calculate inverse inelastic deformation gradient
  for (int i = 0; i < 3; ++i) iFinM(i, i) = 1.0 / isoinelasticdefo;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinTempIso::evaluate_inelastic_def_grad_derivative(
    double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx)
{
  // get parameters
  const double tempgrowthfac = Parameter()->GetTempGrowthFac();
  const double reftemp = Parameter()->RefTemp();

  const double growthfactor = 1.0 + tempgrowthfac * (get_temperature_gp() - reftemp);
  const double scalefac = tempgrowthfac / 3.0 * std::pow(growthfactor, -2.0 / 3.0);

  // prepare identity tensor as 9x1 vector
  static Core::LinAlg::Matrix<9, 1> id9x1(true);
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // here dFindT is zeroed out and filled with the current value
  dFindx.update(scalefac, id9x1, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinTempIso::evaluate_additional_cmat(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<3, 3>& iFinjM,
    const Core::LinAlg::Matrix<6, 1>& iCV, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
    Core::LinAlg::Matrix<6, 6>& cmatadd)
{
  // nothing to do so far, as current growth model is not a function of displacements (and thus C)
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradLinTempIso::EvaluateODStiffMat(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<3, 3>& iFinjM,
    const Core::LinAlg::Matrix<6, 9>& dSdiFinj, Core::LinAlg::Matrix<6, 1>& dstressdT)
{
  static Core::LinAlg::Matrix<9, 1> id9x1(true);
  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters from parameter class
  const double tempgrowthfac = Parameter()->GetTempGrowthFac();
  const double reftemp = Parameter()->RefTemp();

  const double growthfactor = 1.0 + tempgrowthfac * (get_temperature_gp() - reftemp);
  if (growthfactor <= 0.0) FOUR_C_THROW("Determinante of growth must not become negative");

  const double scalefac = -tempgrowthfac / (3.0 * std::pow(growthfactor, 4.0 / 3.0));

  // dstressdT = dSdiFinj : diFinjdT
  // diFinjdT = - growthfac/(3*[1 + growthfac*(T-T_{ref})]^(4/3)) * I
  dstressdT.multiply_nn(scalefac, dSdiFinj, id9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticSource Mat::InelasticDefgradLinTempIso::GetInelasticSource()
{
  return PAR::InelasticSource::temperature;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradNoGrowth::evaluate_additional_cmat(
    const Core::LinAlg::Matrix<3, 3>* defgrad, const Core::LinAlg::Matrix<3, 3>& iFinjM,
    const Core::LinAlg::Matrix<6, 1>& iCV, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
    Core::LinAlg::Matrix<6, 6>& cmatadd)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradNoGrowth::evaluate_inelastic_def_grad_derivative(
    double detjacobian, Core::LinAlg::Matrix<9, 1>& dFindx)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradNoGrowth::evaluate_inverse_inelastic_def_grad(
    const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM)
{
  iFinM = identity_;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradNoGrowth::EvaluateODStiffMat(const Core::LinAlg::Matrix<3, 3>* defgrad,
    const Core::LinAlg::Matrix<3, 3>& iFinjM, const Core::LinAlg::Matrix<6, 9>& dSdiFinj,
    Core::LinAlg::Matrix<6, 1>& dstressdx)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticSource Mat::InelasticDefgradNoGrowth::GetInelasticSource()
{
  return PAR::InelasticSource::none;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::InelasticDefgradNoGrowth::InelasticDefgradNoGrowth(Core::Mat::PAR::Parameter* params)
    : InelasticDefgradFactors(params), identity_(true)
{
  // add 1.0 to main diagonal
  identity_(0, 0) = identity_(1, 1) = identity_(2, 2) = 1.0;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradNoGrowth::pre_evaluate(Teuchos::ParameterList& params, int gp) {}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradTimeFunct::evaluate_inverse_inelastic_def_grad(
    const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM)
{
  const double idetFin = std::pow(funct_value_, -1.0 / 3.0);
  iFinM.update(idetFin, identity_, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::InelasticSource Mat::InelasticDefgradTimeFunct::GetInelasticSource()
{
  return PAR::InelasticSource::none;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::InelasticDefgradTimeFunct::InelasticDefgradTimeFunct(Core::Mat::PAR::Parameter* params)
    : InelasticDefgradFactors(params), funct_value_(0.0), identity_(true)
{
  // add 1.0 to main diagonal
  identity_(0, 0) = identity_(1, 1) = identity_(2, 2) = 1.0;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticDefgradTimeFunct::pre_evaluate(Teuchos::ParameterList& params, int gp)
{
  // evaluate function value for current time step.
  auto& funct = Global::Problem::Instance()->FunctionById<Core::UTILS::FunctionOfTime>(
      Parameter()->FunctNum() - 1);
  const double time = params.get<double>("total time");
  funct_value_ = funct.evaluate(time);
}
FOUR_C_NAMESPACE_CLOSE
