/*----------------------------------------------------------------------*/
/*! \file
\brief evaluation of inelastic deformation gradients and their derivatives

\level 3

*/
/*----------------------------------------------------------------------*/

#include "inelastic_defgrad_factors.H"

#include "electrode.H"
#include "matpar_bundle.H"
#include "material_service.H"
#include "multiplicative_split_defgrad_elasthyper.H"
#include <utility>

#include "drt_globalproblem.H"
#include "voigt_notation.H"

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradNoGrowth::InelasticDefgradNoGrowth(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata)
{
  // do nothing here
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradScalar::InelasticDefgradScalar(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), scalar1_(matdata->GetInt("SCALAR1"))
{
  // safety checks
  // in case not all scatra dofs are transported scalars, the last scatra dof is a potential and can
  // not be treated as a concentration but it is treated like that in so3_scatra_evaluate.cpp in the
  // PreEvaluate method!
  if (scalar1_ != 1) dserror("At the moment it is only possible that SCALAR1 induces growth");
  if (matdata->GetDouble("SCALAR1_RefConc") < 0.0)
    dserror("The reference concentration of SCALAR1 can't be negative");
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradIntercalFrac::InelasticDefgradIntercalFrac(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : InelasticDefgradScalar(matdata), c_max_(-1.0), chi_max_(-1.0)
{
  // get matid
  const int matid = matdata->GetInt("MATID");

  // Check if the material specified by user with MATID is an electrode material
  if (matid > 0)
  {
    // retrieve problem instance to read from
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    // retrieve validated input line of material ID in question
    auto curmat = DRT::Problem::Instance(probinst)->Materials()->ById(matid);
    switch (curmat->Type())
    {
      case INPAR::MAT::m_electrode:
      {
        // Get C_max and Chi_max of electrode material
        c_max_ = curmat->GetDouble("C_MAX");
        chi_max_ = curmat->GetDouble("CHI_MAX");
        break;
      }
      default:
        dserror("The material you have specified by MATID has to be an electrode material!");
    }
  }
  else
  {
    dserror("You have to enter a valid MATID for the corresponding electrode material!");
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradPolyIntercalFrac::InelasticDefgradPolyIntercalFrac(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : InelasticDefgradIntercalFrac(matdata), polynomReferenceValue_(-1.0)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradLinScalarAniso::InelasticDefgradLinScalarAniso(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : InelasticDefgradScalar(matdata),
      growthdir_(Teuchos::rcp(
          new InelasticDeformationDirection(*matdata->Get<std::vector<double>>("GrowthDirection"))))
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradPolyIntercalFracAniso::InelasticDefgradPolyIntercalFracAniso(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : InelasticDefgradPolyIntercalFrac(matdata),
      growthdir_(Teuchos::rcp(
          new InelasticDeformationDirection(*matdata->Get<std::vector<double>>("GrowthDirection"))))
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDeformationDirection::InelasticDeformationDirection(
    std::vector<double> growthdirection)
    : growthdirmat_(true)
{
  if (growthdirection.size() != 3)
  {
    dserror(
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
      growthdirmat_(i, j) = invquadrgrowthdirvecnorm * growthdirection[i] * growthdirection[j];
    }
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradLinTempIso::InelasticDefgradLinTempIso(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      reftemp_(matdata->GetDouble("RefTemp")),
      tempgrowthfac_(matdata->GetDouble("Temp_GrowthFac"))

{
  // safety checks
  if (reftemp_ < 0.0) dserror("Avoid negative reference temperatures");
  if (tempgrowthfac_ == 0.0)
  {
    dserror(
        "Do not use 'MAT_InelasticDefgradLinTempIso' with a growth factor of 0.0. Use "
        "'MAT_InelasticDefgradNoGrowth' instead!");
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradFactors::InelasticDefgradFactors(MAT::PAR::Parameter* params)
    : gp_(-1), params_(params)
{
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Teuchos::RCP<MAT::InelasticDefgradFactors> MAT::InelasticDefgradFactors::Factory(int matnum)
{
  // for the sake of safety
  if (DRT::Problem::Instance()->Materials() == Teuchos::null)
    dserror("List of materials cannot be accessed in the global problem instance.");

  // another safety check
  if (DRT::Problem::Instance()->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // check correct masslin type
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  if (DRT::INPUT::IntegralValue<INPAR::STR::MassLin>(sdyn, "MASSLIN") != INPAR::STR::ml_none)
  {
    dserror(
        "If you use the material 'InelasticDefgradFactors' please set 'MASSLIN' in the "
        "STRUCTURAL DYNAMIC Section to 'None', or feel free to implement other possibility!");
  }

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matnum);

  // get material type and call corresponding constructors
  const INPAR::MAT::MaterialType currentMaterialType = curmat->Type();
  switch (currentMaterialType)
  {
    case INPAR::MAT::mfi_no_growth:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::InelasticDefgradNoGrowth(curmat));

      auto* params = dynamic_cast<MAT::PAR::InelasticDefgradNoGrowth*>(curmat->Parameter());

      return Teuchos::rcp(new InelasticDefgradNoGrowth(params));
    }
    case INPAR::MAT::mfi_lin_scalar_aniso:
    case INPAR::MAT::mfi_lin_scalar_iso:
    {
      const double scalar1MolarGrowthFac = curmat->GetDouble("SCALAR1_MolarGrowthFac");

      // get pointer to linear growth object
      auto linearGrowth = Teuchos::rcp(new InelasticDefgradLinearShape(
          scalar1MolarGrowthFac, curmat->GetDouble("SCALAR1_RefConc")));

      // construct and return pointer to anisotropic version
      if (currentMaterialType == INPAR::MAT::mfi_lin_scalar_aniso)
      {
        if (curmat->Parameter() == nullptr)
          curmat->SetParameter(new MAT::PAR::InelasticDefgradLinScalarAniso(curmat));

        // get pointer to parameter class
        auto* params = dynamic_cast<MAT::PAR::InelasticDefgradLinScalarAniso*>(curmat->Parameter());

        // return pointer to inelastic deformation gradient object
        return Teuchos::rcp(new InelasticDefgradLinScalarAniso(params, linearGrowth));
      }
      // construct and return pointer to isotropic version
      else
      {
        if (curmat->Parameter() == nullptr)
          curmat->SetParameter(new MAT::PAR::InelasticDefgradScalar(curmat));

        // get pointer to parameter class
        auto* params = dynamic_cast<MAT::PAR::InelasticDefgradScalar*>(curmat->Parameter());

        // return pointer to inelastic deformation gradient object
        return Teuchos::rcp(new InelasticDefgradLinScalarIso(params, linearGrowth));
      }
    }
    case INPAR::MAT::mfi_poly_intercal_frac_aniso:
    case INPAR::MAT::mfi_poly_intercal_frac_iso:
    {
      // safety check
      std::vector<double> polyCoeffs(*curmat->Get<std::vector<double>>("POLY_PARAMS"));
      if (polyCoeffs.size() != static_cast<unsigned int>(curmat->GetInt("POLY_PARA_NUM")))
      {
        dserror(
            "Number of coefficients POLY_PARA_NUM you entered in input file has to match the size "
            "of coefficient vector POLY_PARAMS");
      }

      // get pointer to polynomial growth object
      auto polynomialGrowth = Teuchos::rcp(new InelasticDefgradPolynomialShape(
          polyCoeffs, curmat->GetDouble("X_min"), curmat->GetDouble("X_max")));

      // construct and return pointer to anisotropic version
      if (currentMaterialType == INPAR::MAT::mfi_poly_intercal_frac_aniso)
      {
        if (curmat->Parameter() == nullptr)
          curmat->SetParameter(new MAT::PAR::InelasticDefgradPolyIntercalFracAniso(curmat));

        // get pointer to parameter class
        auto* params =
            dynamic_cast<MAT::PAR::InelasticDefgradPolyIntercalFracAniso*>(curmat->Parameter());

        // get reference intercalation fraction
        const double x_ref = MAT::Electrode::ComputeIntercalationFraction(
            curmat->GetDouble("SCALAR1_RefConc"), params->Chimax(), params->Cmax(), 1.0);

        // set the polynomial value in the reference configuration
        params->SetPolynomReferenceValue(polynomialGrowth->ComputePolynomial(x_ref));

        // return pointer to inelastic deformation gradient object
        return Teuchos::rcp(new InelasticDefgradPolyIntercalFracAniso(params, polynomialGrowth));
      }
      // construct and return pointer to isotropic version
      else
      {
        if (curmat->Parameter() == nullptr)
          curmat->SetParameter(new MAT::PAR::InelasticDefgradPolyIntercalFrac(curmat));

        // get pointer to parameter class
        auto* params =
            dynamic_cast<MAT::PAR::InelasticDefgradPolyIntercalFrac*>(curmat->Parameter());

        // get reference intercalation fraction
        const double x_ref = MAT::Electrode::ComputeIntercalationFraction(
            curmat->GetDouble("SCALAR1_RefConc"), params->Chimax(), params->Cmax(), 1.0);

        // set the polynomial value in the reference configuration
        params->SetPolynomReferenceValue(polynomialGrowth->ComputePolynomial(x_ref));

        // return pointer to inelastic deformation gradient object
        return Teuchos::rcp(new InelasticDefgradPolyIntercalFracIso(params, polynomialGrowth));
      }
    }
    case INPAR::MAT::mfi_lin_temp_iso:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::InelasticDefgradLinTempIso(curmat));

      auto* params = dynamic_cast<MAT::PAR::InelasticDefgradLinTempIso*>(curmat->Parameter());
      return Teuchos::rcp(new InelasticDefgradLinTempIso(params));
    }
    default:
      dserror("cannot deal with type %d", curmat->Type());
  }
  // dummy return
  return Teuchos::null;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradScalar::InelasticDefgradScalar(MAT::PAR::Parameter* params)
    : InelasticDefgradFactors(params), concentrations_(Teuchos::null)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradScalar::PreEvaluate(Teuchos::ParameterList& params, const int gp)
{
  SetGP(gp);

  // set pointer to vector of gp_conc, only if gp is 0, because this is the first gp
  if (gp == 0)
    concentrations_ = params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc");
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradScalar::SetConcentrationGP(const double concentration)
{
  // this method is only called for a certain gauss point whose id is not accessible, thus we set a
  // dummy id here and set the corresponding concentration value afterwards
  const int dummy_gp(0);
  SetGP(dummy_gp);

  const int scalar1 = Parameter()->Scalar1();

  GetConcentrationGP().at(scalar1 - 1) = concentration;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolyIntercalFrac::InelasticDefgradPolyIntercalFrac(
    MAT::PAR::Parameter* params, Teuchos::RCP<InelasticDefgradPolynomialShape> polynomial_growth)
    : InelasticDefgradScalar(params), polynomial_growth_(std::move(polynomial_growth))
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradPolyIntercalFrac::EvaluatePolynomial(
    const double concentration, const double detjacobian)
{
  // get intercalation fraction
  const double x = MAT::Electrode::ComputeIntercalationFraction(
      concentration, Parameter()->Chimax(), Parameter()->Cmax(), detjacobian);

  // check bounds of validity of polynomial
  polynomial_growth_->CheckPolynomialBounds(x);

  // calculate and return the value of the polynomial
  return polynomial_growth_->ComputePolynomial(x);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradPolyIntercalFrac::EvaluatePolynomialDerivative(
    const double concentration, const double detjacobian)
{
  // get intercalation fraction
  const double x = MAT::Electrode::ComputeIntercalationFraction(
      concentration, Parameter()->Chimax(), Parameter()->Cmax(), detjacobian);

  return polynomial_growth_->ComputePolynomialDerivative(x);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticSource MAT::InelasticDefgradPolyIntercalFrac::GetInelasticSource()
{
  return MAT::PAR::InelasticSource::concentration;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarIso::InelasticDefgradLinScalarIso(
    MAT::PAR::Parameter* params, Teuchos::RCP<InelasticDefgradLinearShape> linear_growth)
    : InelasticDefgradScalar(params), linear_growth_(std::move(linear_growth))
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticSource MAT::InelasticDefgradLinScalarIso::GetInelasticSource()
{
  return MAT::PAR::InelasticSource::concentration;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarIso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* const defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  // get parameter
  const int sc1 = Parameter()->Scalar1();
  const double material_concentration = GetConcentrationGP().at(sc1 - 1) * defgrad->Determinant();

  // get growth factor
  const double growth_factor = linear_growth_->EvaluateLinearGrowth(material_concentration);

  const double isoinelasticdefo = std::pow(1.0 + growth_factor, (1.0 / 3.0));

  // calculate inverse inelastic deformation gradient
  for (int i = 0; i < 3; ++i) iFinM(i, i) = 1.0 / isoinelasticdefo;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarIso::EvaluateAdditionalCmat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 1>& iCV, const LINALG::Matrix<6, 9>& dSdiFinj,
    LINALG::Matrix<6, 6>& cmatadd)
{
  // static variables
  static LINALG::Matrix<9, 6> diFinjdC(true);
  static LINALG::Matrix<9, 1> id9x1(true);

  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double sc1GrowthFac = linear_growth_->GrowthFac();
  const double concentration = GetConcentrationGP().at(sc1 - 1);
  const double detjacobian = defgrad->Determinant();

  // get growth factor
  const double growth_factor = linear_growth_->EvaluateLinearGrowth(concentration * detjacobian);

  // evaluate scaling factor
  const double scalefac =
      -sc1GrowthFac * concentration * detjacobian / 6.0 * std::pow(1 + growth_factor, -4.0 / 3.0);

  // calculate diFindC
  diFinjdC.MultiplyNT(scalefac, id9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.MultiplyNN(2.0, dSdiFinj, diFinjdC, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarIso::EvaluateODStiffMat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc)
{
  static LINALG::Matrix<9, 1> id9x1(true);
  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double sc1GrowthFac = linear_growth_->GrowthFac();
  const double detjacobian = defgrad->Determinant();
  const double material_concentration = GetConcentrationGP().at(sc1 - 1) * detjacobian;

  // get growth factor
  const double growth_factor = linear_growth_->EvaluateLinearGrowth(material_concentration);

  // calculate scalefac
  const double scalefac =
      -sc1GrowthFac / 3.0 * detjacobian * std::pow(1 + growth_factor, -4.0 / 3.0);

  // calculate diFindc and add contribution to dstressdc = dSdiFinj : diFinjdc
  dstressdc.MultiplyNN(scalefac, dSdiFinj, id9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarIso::EvaluateInelasticDefGradDerivative(
    const double detjacobian, LINALG::Matrix<9, 1>& dFindx)
{
  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double material_concentration = GetConcentrationGP().at(sc1 - 1) * detjacobian;

  // get growth factor
  const double growth_factor = linear_growth_->EvaluateLinearGrowth(material_concentration);
  // calculate the scale factor needed to calculate the derivative below
  const double scalefac = 1.0 / 3.0 * std::pow(1 + growth_factor, -2.0 / 3.0) *
                          linear_growth_->GrowthFac() * detjacobian;

  // prepare identity tensor as 9x1 vector
  static LINALG::Matrix<9, 1> id9x1(true);
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // here dFindc is zeroed out and filled with the current value
  dFindx.Update(scalefac, id9x1, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarAniso::InelasticDefgradLinScalarAniso(
    MAT::PAR::Parameter* params, Teuchos::RCP<InelasticDefgradLinearShape> linear_growth)
    : InelasticDefgradScalar(params), linear_growth_(std::move(linear_growth))
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticSource MAT::InelasticDefgradLinScalarAniso::GetInelasticSource()
{
  return MAT::PAR::InelasticSource::concentration;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarAniso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* const defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  // init and clear variable
  static LINALG::Matrix<3, 3> FinM(true);
  FinM.Clear();

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double material_concentration = GetConcentrationGP().at(sc1 - 1) * defgrad->Determinant();

  // get growth factor
  const double growth_factor = linear_growth_->EvaluateLinearGrowth(material_concentration);

  // calculate inelastic deformation gradient and its inverse
  for (int i = 0; i < 3; ++i) FinM(i, i) = 1.0;

  // finalize inelastic deformation gradient matrix (FinM is calculated, such that the volume change
  // is a linear function of the scalar (mapped to reference frame) that causes it)
  FinM.Update(growth_factor, Parameter()->Growthdirmat(), 1.0);

  // calculate inverse of inelastic deformation gradient matrix
  iFinM.Invert(FinM);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarAniso::EvaluateAdditionalCmat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 1>& iCV, const LINALG::Matrix<6, 9>& dSdiFinj,
    LINALG::Matrix<6, 6>& cmatadd)
{
  static LINALG::Matrix<3, 3> temp(true);
  static LINALG::Matrix<3, 3> iFinjGiFinj(true);
  static LINALG::Matrix<9, 1> iFinjGiFinj9x1(true);
  static LINALG::Matrix<9, 6> diFinjdC(true);

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double sc1GrowthFac = linear_growth_->GrowthFac();
  const double concentration = GetConcentrationGP().at(sc1 - 1);
  const double detjacobian = defgrad->Determinant();

  // prepare scalefac
  const double scalefac = -sc1GrowthFac * concentration * detjacobian / 2.0;

  // calculate F_{in,j}^{-1} . G . F_{in,j}^{-1} with F_{in,j}, the j-th factor of F_{in}
  temp.MultiplyNN(1.0, iFinjM, Parameter()->Growthdirmat(), 0.0);
  iFinjGiFinj.MultiplyNN(1.0, temp, iFinjM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(iFinjGiFinj, iFinjGiFinj9x1);

  // calculate diFinjdC
  diFinjdC.MultiplyNT(scalefac, iFinjGiFinj9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.MultiplyNN(2.0, dSdiFinj, diFinjdC, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarAniso::EvaluateODStiffMat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc)
{
  // static variables
  static LINALG::Matrix<3, 3> tmp(true);
  static LINALG::Matrix<3, 3> diFinjdcM(true);
  static LINALG::Matrix<9, 1> diFinjdc9x1(true);

  // get parameters
  const double sc1GrowthFac = linear_growth_->GrowthFac();
  const double detjacobian = defgrad->Determinant();

  // prepare scalefac
  const double scalefac = -sc1GrowthFac * detjacobian;

  // calculate diFinjdc
  tmp.MultiplyNN(1.0, iFinjM, Parameter()->Growthdirmat(), 0.0);
  diFinjdcM.MultiplyNN(scalefac, tmp, iFinjM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(diFinjdcM, diFinjdc9x1);

  // dstressdc = dSdiFinj : diFinjdc
  dstressdc.MultiplyNN(1.0, dSdiFinj, diFinjdc9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarAniso::EvaluateInelasticDefGradDerivative(
    const double detjacobian, LINALG::Matrix<9, 1>& dFindx)
{
  const double scalefac = linear_growth_->GrowthFac() * detjacobian;

  // get the growth direction matrix as a 9x1 vector
  static LINALG::Matrix<9, 1> growthdirmat9x1(true);
  UTILS::VOIGT::Matrix3x3to9x1(Parameter()->Growthdirmat(), growthdirmat9x1);

  // here dFindc is zeroed out and filled with the current value
  dFindx.Update(scalefac, growthdirmat9x1, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolyIntercalFracIso::InelasticDefgradPolyIntercalFracIso(
    MAT::PAR::Parameter* params,
    const Teuchos::RCP<InelasticDefgradPolynomialShape>& polynomial_growth)
    : InelasticDefgradPolyIntercalFrac(params, polynomial_growth)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyIntercalFracIso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* const defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double polynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get polynomial
  const double polynomValue =
      EvaluatePolynomial(GetConcentrationGP().at(sc1 - 1), defgrad->Determinant());

  // calculate growth
  const double isoInelasticDefo =
      std::pow((1.0 + polynomValue) / (1.0 + polynomReferenceValue), (1.0 / 3.0));
  // calculate inverse inelastic deformation gradient
  for (int i = 0; i < 3; ++i) iFinM(i, i) = 1.0 / isoInelasticDefo;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyIntercalFracIso::EvaluateAdditionalCmat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 1>& iCV, const LINALG::Matrix<6, 9>& dSdiFinj,
    LINALG::Matrix<6, 6>& cmatadd)
{
  // static variables
  static LINALG::Matrix<9, 6> diFinjdC(true);
  static LINALG::Matrix<9, 1> id9x1(true);

  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double chi_max = Parameter()->Chimax();
  const double c_max = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();
  const double concentration = GetConcentrationGP().at(sc1 - 1);
  const double polynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get polynomials
  const double polynomValue = EvaluatePolynomial(concentration, detjacobian);
  const double polynomDerivativeValue = EvaluatePolynomialDerivative(concentration, detjacobian);

  // prepare scalefac
  const double scalefac = -1.0 / (6.0 * c_max) * concentration * chi_max * detjacobian *
                          std::pow(1.0 + polynomValue, -4.0 / 3.0) * polynomDerivativeValue *
                          std::pow(1.0 + polynomReferenceValue, 1.0 / 3.0);

  // calculate diFinjdC
  diFinjdC.MultiplyNT(scalefac, id9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.MultiplyNN(2.0, dSdiFinj, diFinjdC, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyIntercalFracIso::EvaluateODStiffMat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc)
{
  static LINALG::Matrix<9, 1> id9x1(true);
  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double concentration = GetConcentrationGP().at(sc1 - 1);
  const double detjacobian = defgrad->Determinant();
  const double polynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get polynomial and derivatives
  const double polynomValue = EvaluatePolynomial(concentration, detjacobian);
  const double polynomDerivativeValue = EvaluatePolynomialDerivative(concentration, detjacobian);
  const double dChidc = MAT::Electrode::ComputeIntercalationFractionDerivative(
      Parameter()->Chimax(), Parameter()->Cmax(), detjacobian);

  // prepare scalefac
  const double scalefac = -1.0 / 3.0 * std::pow(1.0 + polynomValue, -4.0 / 3.0) *
                          std::pow(1.0 + polynomReferenceValue, 1.0 / 3.0) *
                          polynomDerivativeValue * dChidc;

  // calculate diFinjdc and add contribution to dstressdc = dSdiFinj : diFinjdc
  dstressdc.MultiplyNN(scalefac, dSdiFinj, id9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyIntercalFracIso::EvaluateInelasticDefGradDerivative(
    const double detjacobian, LINALG::Matrix<9, 1>& dFindx)
{
  static LINALG::Matrix<9, 1> id9x1(true);
  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double concentration = GetConcentrationGP().at(sc1 - 1);
  const double polynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get polynomial and its derivative
  const double polynomValue = EvaluatePolynomial(concentration, detjacobian);
  const double polynomDerivativeValue = EvaluatePolynomialDerivative(concentration, detjacobian);

  // calculate the scale factor needed to get the derivative later
  const double denominator = 1.0 / (polynomReferenceValue + 1.0);
  const double base = (polynomValue + 1.0) * denominator;
  const double dChidc = MAT::Electrode::ComputeIntercalationFractionDerivative(
      Parameter()->Chimax(), Parameter()->Cmax(), detjacobian);
  const double scalefac =
      1.0 / 3.0 * std::pow(base, -2.0 / 3.0) * polynomDerivativeValue * denominator * dChidc;

  // here dFindc is zeroed out and filled with the current value
  dFindx.Update(scalefac, id9x1, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolyIntercalFracAniso::InelasticDefgradPolyIntercalFracAniso(
    MAT::PAR::Parameter* params,
    const Teuchos::RCP<InelasticDefgradPolynomialShape>& polynomial_growth)
    : InelasticDefgradPolyIntercalFrac(params, polynomial_growth)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyIntercalFracAniso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* const defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  // init and clear variable
  static LINALG::Matrix<3, 3> FinM(true);
  FinM.Clear();

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double polynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get polynomials
  const double polynomValue =
      EvaluatePolynomial(GetConcentrationGP().at(sc1 - 1), defgrad->Determinant());

  // calculate growth factor
  const double growth_factor =
      (polynomValue - polynomReferenceValue) / (polynomReferenceValue + 1.0);

  // calculate inelastic deformation gradient and its inverse
  for (int i = 0; i < 3; ++i) FinM(i, i) = 1.0;

  // add the growth part
  FinM.Update(growth_factor, Parameter()->Growthdirmat(), 1.0);

  // calculate inverse of inelastic deformation gradient matrix
  iFinM.Invert(FinM);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyIntercalFracAniso::EvaluateAdditionalCmat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 1>& iCV, const LINALG::Matrix<6, 9>& dSdiFinj,
    LINALG::Matrix<6, 6>& cmatadd)
{
  static LINALG::Matrix<3, 3> temp(true);
  static LINALG::Matrix<3, 3> iFinjGiFinj(true);
  static LINALG::Matrix<9, 1> iFinjGiFinj9x1(true);
  static LINALG::Matrix<9, 6> diFinjdC(true);

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double chi_max = Parameter()->Chimax();
  const double c_max = Parameter()->Cmax();
  const double concentration = GetConcentrationGP().at(sc1 - 1);
  const double detjacobian = defgrad->Determinant();
  const double polynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get first derivative of polynomial
  const double polynomDerivativeValue = EvaluatePolynomialDerivative(concentration, detjacobian);

  // prepare scalefac
  const double scalefac = -detjacobian * concentration * chi_max * polynomDerivativeValue /
                          (2.0 * c_max * (polynomReferenceValue + 1.0));

  // calculate F_{in,j}^{-1} . G . F_{in,j}^{-1} with F_{in,j}, the j-th factor of F_{in}
  temp.MultiplyNN(1.0, iFinjM, Parameter()->Growthdirmat(), 0.0);
  iFinjGiFinj.MultiplyNN(1.0, temp, iFinjM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(iFinjGiFinj, iFinjGiFinj9x1);

  // calculate diFinjdC
  diFinjdC.MultiplyNT(scalefac, iFinjGiFinj9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.MultiplyNN(2.0, dSdiFinj, diFinjdC, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyIntercalFracAniso::EvaluateODStiffMat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc)
{
  // static variables
  static LINALG::Matrix<3, 3> tmp(true);
  static LINALG::Matrix<3, 3> diFinjdcM(true);
  static LINALG::Matrix<9, 1> diFinjdc9x1(true);

  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double detjacobian = defgrad->Determinant();
  const double polynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get derivatives
  const double polynomDerivativeValue =
      EvaluatePolynomialDerivative(GetConcentrationGP().at(sc1 - 1), detjacobian);
  const double dChidc = MAT::Electrode::ComputeIntercalationFractionDerivative(
      Parameter()->Chimax(), Parameter()->Cmax(), detjacobian);

  // prepare scalefac
  const double scalefac = -polynomDerivativeValue / (polynomReferenceValue + 1.0) * dChidc;

  // calculate diFinjdc
  tmp.MultiplyNN(1.0, iFinjM, Parameter()->Growthdirmat(), 0.0);
  diFinjdcM.MultiplyNN(scalefac, tmp, iFinjM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(diFinjdcM, diFinjdc9x1);

  // dstressdc = dSdiFinj : diFinjdc
  dstressdc.MultiplyNN(1.0, dSdiFinj, diFinjdc9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyIntercalFracAniso::EvaluateInelasticDefGradDerivative(
    const double detjacobian, LINALG::Matrix<9, 1>& dFindx)
{
  // get parameters
  const int sc1 = Parameter()->Scalar1();
  const double concentration = GetConcentrationGP().at(sc1 - 1);
  const double polynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get polynomial derivative
  const double polynomDerivativeValue = EvaluatePolynomialDerivative(concentration, detjacobian);

  const double dChidc = MAT::Electrode::ComputeIntercalationFractionDerivative(
      Parameter()->Chimax(), Parameter()->Cmax(), detjacobian);
  const double scalefac = polynomDerivativeValue / (polynomReferenceValue + 1.0) * dChidc;

  // get the growth direction matrix as a 9x1 vector
  static LINALG::Matrix<9, 1> growthdirmat9x1(true);
  UTILS::VOIGT::Matrix3x3to9x1(Parameter()->Growthdirmat(), growthdirmat9x1);

  // here dFindc is zeroed out and filled with the current value
  dFindx.Update(scalefac, growthdirmat9x1, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinearShape::InelasticDefgradLinearShape(
    const double growthFac, const double referenceValue)
    : growthFac_(growthFac), referenceValue_(referenceValue)
{
  // safety checks
  if (growthFac < 0.0) dserror("Growth factor can not be negative, please check your input file!");
  if (growthFac == 0.0)
  {
    dserror(
        "Do not use linear growth laws with a growth factor of 0.0. Use "
        "'MAT_InelasticDefgradNoGrowth' instead!");
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradLinearShape::EvaluateLinearGrowth(const double value) const
{
  // calculate and return the linear growth factor
  return growthFac_ * (value - referenceValue_);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolynomialShape::InelasticDefgradPolynomialShape(
    std::vector<double> poly_coeffs, const double x_min, const double x_max)
    : poly_coeffs_(std::move(poly_coeffs)), x_min_(x_min), x_max_(x_max)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradPolynomialShape::ComputePolynomial(const double x)
{
  // initialize the variable for the evaluation of the polynomial
  double polynom(0.0);

  // compute polynomial
  for (unsigned i = 0; i < poly_coeffs_.size(); ++i) polynom += poly_coeffs_[i] * std::pow(x, i);

  return polynom;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradPolynomialShape::ComputePolynomialDerivative(const double x)
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
void MAT::InelasticDefgradPolynomialShape::CheckPolynomialBounds(const double x) const
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
MAT::InelasticDefgradLinTempIso::InelasticDefgradLinTempIso(MAT::PAR::Parameter* params)
    : InelasticDefgradFactors(params), temperatures_(Teuchos::null)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinTempIso::PreEvaluate(Teuchos::ParameterList& params, int gp)
{
  // get Gauss point number
  SetGP(gp);

  // set pointer to vector of gp_temp, only if gp is 0, because this is the first gp
  if (gp == 0) temperatures_ = params.get<Teuchos::RCP<std::vector<double>>>("gp_temp");
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinTempIso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  // get parameters
  const double tempgrowthfac = Parameter()->GetTempGrowthFac();
  const double reftemp = Parameter()->RefTemp();

  const double growthfactor = 1.0 + tempgrowthfac * (GetTemperatureGP() - reftemp);
  if (growthfactor <= 0.0) dserror("Determinante of growth must not become negative");
  const double isoinelasticdefo = std::pow(growthfactor, (1.0 / 3.0));

  // calculate inverse inelastic deformation gradient
  for (int i = 0; i < 3; ++i) iFinM(i, i) = 1.0 / isoinelasticdefo;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinTempIso::EvaluateInelasticDefGradDerivative(
    double detjacobian, LINALG::Matrix<9, 1>& dFindx)
{
  // get parameters
  const double tempgrowthfac = Parameter()->GetTempGrowthFac();
  const double reftemp = Parameter()->RefTemp();

  const double growthfactor = 1.0 + tempgrowthfac * (GetTemperatureGP() - reftemp);
  const double scalefac = tempgrowthfac / 3.0 * std::pow(growthfactor, -2.0 / 3.0);

  // prepare identity tensor as 9x1 vector
  static LINALG::Matrix<9, 1> id9x1(true);
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // here dFindT is zeroed out and filled with the current value
  dFindx.Update(scalefac, id9x1, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinTempIso::EvaluateAdditionalCmat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 1>& iCV, const LINALG::Matrix<6, 9>& dSdiFinj,
    LINALG::Matrix<6, 6>& cmatadd)
{
  // nothing to do so far, as current growth model is not a function of displacements (and thus C)
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinTempIso::EvaluateODStiffMat(const LINALG::Matrix<3, 3>* const defgrad,
    const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 9>& dSdiFinj,
    LINALG::Matrix<6, 1>& dstressdT)
{
  static LINALG::Matrix<9, 1> id9x1(true);
  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters from parameter class
  const double tempgrowthfac = Parameter()->GetTempGrowthFac();
  const double reftemp = Parameter()->RefTemp();

  const double growthfactor = 1.0 + tempgrowthfac * (GetTemperatureGP() - reftemp);
  if (growthfactor <= 0.0) dserror("Determinante of growth must not become negative");

  const double scalefac = -tempgrowthfac / (3.0 * std::pow(growthfactor, 4.0 / 3.0));

  // dstressdT = dSdiFinj : diFinjdT
  // diFinjdT = - growthfac/(3*[1 + growthfac*(T-T_{ref})]^(4/3)) * I
  dstressdT.MultiplyNN(scalefac, dSdiFinj, id9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticSource MAT::InelasticDefgradLinTempIso::GetInelasticSource()
{
  return PAR::InelasticSource::temperature;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradNoGrowth::EvaluateAdditionalCmat(const LINALG::Matrix<3, 3>* defgrad,
    const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 1>& iCV,
    const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 6>& cmatadd)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradNoGrowth::EvaluateInelasticDefGradDerivative(
    double detjacobian, LINALG::Matrix<9, 1>& dFindx)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradNoGrowth::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  iFinM = identity_;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradNoGrowth::EvaluateODStiffMat(const LINALG::Matrix<3, 3>* defgrad,
    const LINALG::Matrix<3, 3>& iFinjM, const LINALG::Matrix<6, 9>& dSdiFinj,
    LINALG::Matrix<6, 1>& dstressdx)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticSource MAT::InelasticDefgradNoGrowth::GetInelasticSource()
{
  return PAR::InelasticSource::none;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradNoGrowth::InelasticDefgradNoGrowth(MAT::PAR::Parameter* params)
    : InelasticDefgradFactors(params), identity_(true)
{
  // add 1.0 to main diagonal
  identity_(0, 0) = identity_(1, 1) = identity_(2, 2) = 1.0;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradNoGrowth::PreEvaluate(Teuchos::ParameterList& params, int gp) {}