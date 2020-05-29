/*----------------------------------------------------------------------*/
/*! \file
\brief evaluation of inelastic deformation gradients and their derivatives

\level 3

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/

#include "inelastic_defgrad_factors.H"
#include "electrode.H"
#include "matpar_bundle.H"
#include "material_service.H"
#include "multiplicative_split_defgrad_elasthyper.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/voigt_notation.H"

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradScalar::InelasticDefgradScalar(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), Scalar1_(matdata->GetInt("SCALAR1"))
{
  // safety checks
  // in case not all scatra dofs are transported scalars, the last scatra dof is a potential and can
  // not be treated as a concentration but it is treated like that in so3_scatra_evaluate.cpp in the
  // PreEvaluate method!
  if (Scalar1_ != 1) dserror("At the moment it is only possible that SCALAR1 induces growth");
  if (matdata->GetDouble("SCALAR1_RefConc") < 0.0)
    dserror("The reference concentration of SCALAR1 can't be negative");
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradIntercalFrac::InelasticDefgradIntercalFrac(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : InelasticDefgradScalar(matdata), Cmax_(-1.0), Chimax_(-1.0)
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
        Cmax_ = curmat->GetDouble("C_MAX");
        Chimax_ = curmat->GetDouble("CHI_MAX");
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
    : InelasticDefgradIntercalFrac(matdata), PolynomReferenceValue_(-1.0)
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
    dserror(
        "Since we have a 3D problem here, vector that defines the growth direction also needs to "
        "have the size 3!");

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
MAT::InelasticDefgradFactors::InelasticDefgradFactors() : params_(nullptr) {}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradFactors::InelasticDefgradFactors(MAT::PAR::Parameter* params) : params_(params)
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
    dserror(
        "If you use the material 'InelasticDefgradFactors' please set 'MASSLIN' in the "
        "STRUCTURAL DYNAMIC Section to 'None', or feel free to implement other possibility!");

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matnum);

  // get material type and call corresponding constructors
  const INPAR::MAT::MaterialType CurrentMaterialType = curmat->Type();
  switch (CurrentMaterialType)
  {
    case INPAR::MAT::mfi_lin_scalar_aniso:
    case INPAR::MAT::mfi_lin_scalar_iso:
    {
      // safety check
      const double Scalar1MolarGrwothFac = curmat->GetDouble("SCALAR1_MolarGrowthFac");
      if (Scalar1MolarGrwothFac < 0.0)
        dserror("Growth factor can not be negative, please check your input file!");

      // get pointer to linear growth object
      auto LinearGrowth = Teuchos::rcp(new InelasticDefgradLinearShape(
          Scalar1MolarGrwothFac, curmat->GetDouble("SCALAR1_RefConc")));

      // construct and return pointer to anisotropic version
      if (CurrentMaterialType == INPAR::MAT::mfi_lin_scalar_aniso)
      {
        if (curmat->Parameter() == nullptr)
          curmat->SetParameter(new MAT::PAR::InelasticDefgradLinScalarAniso(curmat));

        // get pointer to parameter class
        auto* params = dynamic_cast<MAT::PAR::InelasticDefgradLinScalarAniso*>(curmat->Parameter());

        // return pointer to inelastic deformation gradient object
        return Teuchos::rcp(new InelasticDefgradLinScalarAniso(params, LinearGrowth));
      }
      // construct and return pointer to isotropic version
      else
      {
        if (curmat->Parameter() == nullptr)
          curmat->SetParameter(new MAT::PAR::InelasticDefgradScalar(curmat));

        // get pointer to parameter class
        auto* params = dynamic_cast<MAT::PAR::InelasticDefgradScalar*>(curmat->Parameter());

        // return pointer to inelastic deformation gradient object
        return Teuchos::rcp(new InelasticDefgradLinScalarIso(params, LinearGrowth));
      }
    }
    case INPAR::MAT::mfi_poly_intercal_frac_aniso:
    case INPAR::MAT::mfi_poly_intercal_frac_iso:
    {
      // safety check
      std::vector<double> PolyCoeffs(*curmat->Get<std::vector<double>>("POLY_PARAMS"));
      if (PolyCoeffs.size() != static_cast<unsigned int>(curmat->GetInt("POLY_PARA_NUM")))
        dserror(
            "Number of coefficients POLY_PARA_NUM you entered in input file has to match the size "
            "of coefficient vector POLY_PARAMS");

      // get pointer to polynomial growth object
      auto PolynomialGrowth = Teuchos::rcp(new InelasticDefgradPolynomialShape(
          PolyCoeffs, curmat->GetDouble("X_min"), curmat->GetDouble("X_max")));

      // construct and return pointer to anisotropic version
      if (CurrentMaterialType == INPAR::MAT::mfi_poly_intercal_frac_aniso)
      {
        if (curmat->Parameter() == nullptr)
          curmat->SetParameter(new MAT::PAR::InelasticDefgradPolyIntercalFracAniso(curmat));

        // get pointer to parameter class
        auto* params =
            dynamic_cast<MAT::PAR::InelasticDefgradPolyIntercalFracAniso*>(curmat->Parameter());

        // get reference intercalation fraction
        const double X_ref = MAT::Electrode::ComputeIntercalationFraction(
            curmat->GetDouble("SCALAR1_RefConc"), params->Chimax(), params->Cmax(), 1.0);

        // set the polynomial value in the reference configuration
        params->SetPolynomReferenceValue(PolynomialGrowth->ComputePolynomial(X_ref));

        // return pointer to inelastic deformation gradient object
        return Teuchos::rcp(new InelasticDefgradPolyIntercalFracAniso(params, PolynomialGrowth));
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
        const double X_ref = MAT::Electrode::ComputeIntercalationFraction(
            curmat->GetDouble("SCALAR1_RefConc"), params->Chimax(), params->Cmax(), 1.0);

        // set the polynomial value in the reference configuration
        params->SetPolynomReferenceValue(PolynomialGrowth->ComputePolynomial(X_ref));

        // return pointer to inelastic deformation gradient object
        return Teuchos::rcp(new InelasticDefgradPolyIntercalFracIso(params, PolynomialGrowth));
      }
    }
    default:
      dserror("cannot deal with type %d", curmat->Type());
  }
  // dummy return
  return Teuchos::null;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolyIntercalFrac::InelasticDefgradPolyIntercalFrac(MAT::PAR::Parameter* params,
    const Teuchos::RCP<InelasticDefgradPolynomialShape>& PolynomialGrowth)
    : InelasticDefgradFactors(params), polynomial_growth_(PolynomialGrowth)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradPolyIntercalFrac::EvaluatePolynomial(
    const double concentration, const double DetJacobian)
{
  // get intercalation fraction
  const double X = MAT::Electrode::ComputeIntercalationFraction(
      concentration, Parameter()->Chimax(), Parameter()->Cmax(), DetJacobian);

  // check bounds of validity of polynomial
  polynomial_growth_->CheckPolynomialBounds(X);

  // calculate and return the value of the polynomial
  return polynomial_growth_->ComputePolynomial(X);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradPolyIntercalFrac::EvaluatePolynomialDerivative(
    const double concentration, const double DetJacobian)
{
  // get intercalation fraction
  const double X = MAT::Electrode::ComputeIntercalationFraction(
      concentration, Parameter()->Chimax(), Parameter()->Cmax(), DetJacobian);

  return polynomial_growth_->ComputePolynomialDerivative(X);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticSource MAT::InelasticDefgradPolyIntercalFrac::GetInelasticSource()
{
  return MAT::PAR::InelasticSource::inelastic_concentration;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarIso::InelasticDefgradLinScalarIso()
    : InelasticDefgradFactors(),
      gp_(-1),
      concentrations_(Teuchos::null),
      linear_growth_(Teuchos::null)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarIso::InelasticDefgradLinScalarIso(
    MAT::PAR::Parameter* params, const Teuchos::RCP<InelasticDefgradLinearShape>& LinearGrowth)
    : InelasticDefgradFactors(params),
      gp_(-1),
      concentrations_(Teuchos::null),
      linear_growth_(LinearGrowth)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticSource MAT::InelasticDefgradLinScalarIso::GetInelasticSource()
{
  return MAT::PAR::InelasticSource::inelastic_concentration;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarIso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* const defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  // get parameter
  const int Sc1 = Parameter()->Scalar1();
  const double material_concentration =
      concentrations_->at(gp_).at(Sc1 - 1) * defgrad->Determinant();

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
  const int Sc1 = Parameter()->Scalar1();
  const double Sc1GrowthFac = linear_growth_->GrowthFac();
  const double concentration = concentrations_->at(gp_).at(Sc1 - 1);
  const double detjacobian = defgrad->Determinant();

  // get growth factor
  const double growth_factor = linear_growth_->EvaluateLinearGrowth(concentration * detjacobian);

  // evaluate scaling factor
  const double scalefac =
      -Sc1GrowthFac * concentration * detjacobian / 6.0 * std::pow(1 + growth_factor, -4.0 / 3.0);

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
  const int Sc1 = Parameter()->Scalar1();
  const double Sc1GrowthFac = linear_growth_->GrowthFac();
  const double detjacobian = defgrad->Determinant();
  const double material_concentration = concentrations_->at(gp_).at(Sc1 - 1) * detjacobian;

  // get growth factor
  const double growth_factor = linear_growth_->EvaluateLinearGrowth(material_concentration);

  // calculate scalefac
  const double scalefac =
      -Sc1GrowthFac / 3.0 * detjacobian * std::pow(1 + growth_factor, -4.0 / 3.0);

  // calculate diFindc and add contribution to dstressdc = dSdiFinj : diFinjdc
  dstressdc.MultiplyNN(scalefac, dSdiFinj, id9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarIso::PreEvaluate(Teuchos::ParameterList& params)
{
  // get Gauss point number
  gp_ = params.get<int>("gp", -1);
  if (gp_ == -1) dserror("No Gauss point number provided in material.");

  // set pointer to vector of gp_conc, only if gp is 0, because this is the first gp
  if (gp_ == 0)
    concentrations_ = params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc",
        Teuchos::rcp(new std::vector<std::vector<double>>(30, std::vector<double>(20, 0.0))));
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarAniso::InelasticDefgradLinScalarAniso()
    : InelasticDefgradFactors(),
      gp_(-1),
      concentrations_(Teuchos::null),
      linear_growth_(Teuchos::null)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarAniso::InelasticDefgradLinScalarAniso(
    MAT::PAR::Parameter* params, const Teuchos::RCP<InelasticDefgradLinearShape>& LinearGrowth)
    : InelasticDefgradFactors(params),
      gp_(-1),
      concentrations_(Teuchos::null),
      linear_growth_(LinearGrowth)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticSource MAT::InelasticDefgradLinScalarAniso::GetInelasticSource()
{
  return MAT::PAR::InelasticSource::inelastic_concentration;
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
  const int Sc1 = Parameter()->Scalar1();
  const double material_concentration =
      concentrations_->at(gp_).at(Sc1 - 1) * defgrad->Determinant();

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
  const int Sc1 = Parameter()->Scalar1();
  const double Sc1GrowthFac = linear_growth_->GrowthFac();
  const double concentration = concentrations_->at(gp_).at(Sc1 - 1);
  const double detjacobian = defgrad->Determinant();

  // prepare scalefac
  const double scalefac = -Sc1GrowthFac * concentration * detjacobian / 2.0;

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
  const double Sc1GrowthFac = linear_growth_->GrowthFac();
  const double detjacobian = defgrad->Determinant();

  // prepare scalefac
  const double scalefac = -Sc1GrowthFac * detjacobian;

  // calculate diFinjdc
  tmp.MultiplyNN(1.0, iFinjM, Parameter()->Growthdirmat(), 0.0);
  diFinjdcM.MultiplyNN(scalefac, tmp, iFinjM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(diFinjdcM, diFinjdc9x1);

  // dstressdc = dSdiFinj : diFinjdc
  dstressdc.MultiplyNN(1.0, dSdiFinj, diFinjdc9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarAniso::PreEvaluate(Teuchos::ParameterList& params)
{
  // get Gauss point number
  gp_ = params.get<int>("gp", -1);
  if (gp_ == -1) dserror("No Gauss point number provided in material.");

  // set pointer to vector of gp_conc, only if gp is 0, because this is the first gp
  if (gp_ == 0)
    concentrations_ = params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc",
        Teuchos::rcp(new std::vector<std::vector<double>>(30, std::vector<double>(20, 0.0))));
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolyIntercalFracIso::InelasticDefgradPolyIntercalFracIso()
    : InelasticDefgradPolyIntercalFrac(), gp_(-1), concentrations_(Teuchos::null)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolyIntercalFracIso::InelasticDefgradPolyIntercalFracIso(
    MAT::PAR::Parameter* params,
    const Teuchos::RCP<InelasticDefgradPolynomialShape>& PolynomialGrowth)
    : InelasticDefgradPolyIntercalFrac(params, PolynomialGrowth),
      gp_(-1),
      concentrations_(Teuchos::null)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyIntercalFracIso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* const defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  // get parameters
  const int Sc1 = Parameter()->Scalar1();
  const double PolynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get polynomial
  const double PolynomValue =
      EvaluatePolynomial(concentrations_->at(gp_).at(Sc1 - 1), defgrad->Determinant());

  // calculate growth
  const double isoinelasticdefo =
      std::pow((1.0 + PolynomValue) / (1.0 + PolynomReferenceValue), (1.0 / 3.0));
  // calculate inverse inelastic deformation gradient
  for (int i = 0; i < 3; ++i) iFinM(i, i) = 1.0 / isoinelasticdefo;
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
  const int Sc1 = Parameter()->Scalar1();
  const double Chimax = Parameter()->Chimax();
  const double Cmax = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();
  const double concentration = concentrations_->at(gp_).at(Sc1 - 1);
  const double PolynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get polynomials
  const double PolynomValue = EvaluatePolynomial(concentration, detjacobian);
  const double PolynomDerivativeValue = EvaluatePolynomialDerivative(concentration, detjacobian);

  // prepare scalefac
  const double scalefac = -1.0 / (6.0 * Cmax) * concentration * Chimax * detjacobian *
                          std::pow(1.0 + PolynomValue, -4.0 / 3.0) * PolynomDerivativeValue *
                          std::pow(1.0 + PolynomReferenceValue, 1.0 / 3.0);

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
  const int Sc1 = Parameter()->Scalar1();
  const double concentration = concentrations_->at(gp_).at(Sc1 - 1);
  const double Chimax = Parameter()->Chimax();
  const double Cmax = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();
  const double PolynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get polynomial
  const double PolynomValue = EvaluatePolynomial(concentration, detjacobian);
  const double PolynomDerivativeValue = EvaluatePolynomialDerivative(concentration, detjacobian);

  // prepare scalefac
  const double scalefac = -1.0 / (3.0 * Cmax) * Chimax * detjacobian * PolynomDerivativeValue *
                          std::pow(1.0 + PolynomReferenceValue, 1.0 / 3.0) *
                          std::pow(1.0 + PolynomValue, -4.0 / 3.0);

  // calculate diFinjdc and add contribution to dstressdc = dSdiFinj : diFinjdc
  dstressdc.MultiplyNN(scalefac, dSdiFinj, id9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyIntercalFracIso::PreEvaluate(Teuchos::ParameterList& params)
{
  // get Gauss point number
  gp_ = params.get<int>("gp", -1);
  if (gp_ == -1) dserror("No Gauss point number provided in material.");

  // set pointer to vector of gp_conc, only if gp is 0, because this is the first gp
  if (gp_ == 0)
    concentrations_ = params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc",
        Teuchos::rcp(new std::vector<std::vector<double>>(30, std::vector<double>(20, 0.0))));
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolyIntercalFracAniso::InelasticDefgradPolyIntercalFracAniso()
    : InelasticDefgradPolyIntercalFrac(), gp_(-1), concentrations_(Teuchos::null)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolyIntercalFracAniso::InelasticDefgradPolyIntercalFracAniso(
    MAT::PAR::Parameter* params,
    const Teuchos::RCP<InelasticDefgradPolynomialShape>& PolynomialGrowth)
    : InelasticDefgradPolyIntercalFrac(params, PolynomialGrowth),
      gp_(-1),
      concentrations_(Teuchos::null)
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
  const int Sc1 = Parameter()->Scalar1();
  const double PolynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get polynomials
  const double PolynomValue =
      EvaluatePolynomial(concentrations_->at(gp_).at(Sc1 - 1), defgrad->Determinant());

  // calculate growth factor
  const double growth_factor =
      (PolynomValue - PolynomReferenceValue) / (PolynomReferenceValue + 1.0);

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
  const int Sc1 = Parameter()->Scalar1();
  const double Chimax = Parameter()->Chimax();
  const double Cmax = Parameter()->Cmax();
  const double concentration = concentrations_->at(gp_).at(Sc1 - 1);
  const double detjacobian = defgrad->Determinant();
  const double PolynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get first derivative of polynomial
  const double PolynomDerivativeValue = EvaluatePolynomialDerivative(concentration, detjacobian);

  // prepare scalefac
  const double scalefac = -detjacobian * concentration * Chimax * PolynomDerivativeValue /
                          (2.0 * Cmax * (PolynomReferenceValue + 1.0));

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
  const int Sc1 = Parameter()->Scalar1();
  const double Chimax = Parameter()->Chimax();
  const double Cmax = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();
  const double PolynomReferenceValue = Parameter()->GetPolynomReferenceValue();

  // get polynomials
  const double PolynomDerivativeValue =
      EvaluatePolynomialDerivative(concentrations_->at(gp_).at(Sc1 - 1), detjacobian);

  // prepare scalefac
  const double scalefac =
      -detjacobian * Chimax * PolynomDerivativeValue / (Cmax * (PolynomReferenceValue + 1.0));

  // calculate diFinjdc
  tmp.MultiplyNN(1.0, iFinjM, Parameter()->Growthdirmat(), 0.0);
  diFinjdcM.MultiplyNN(scalefac, tmp, iFinjM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(diFinjdcM, diFinjdc9x1);

  // dstressdc = dSdiFinj : diFinjdc
  dstressdc.MultiplyNN(1.0, dSdiFinj, diFinjdc9x1, 1.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyIntercalFracAniso::PreEvaluate(Teuchos::ParameterList& params)
{
  // get Gauss point number
  gp_ = params.get<int>("gp", -1);
  if (gp_ == -1) dserror("No Gauss point number provided in material.");

  // set pointer to vector of gp_conc, only if gp is 0, because this is the first gp
  if (gp_ == 0)
    concentrations_ = params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc",
        Teuchos::rcp(new std::vector<std::vector<double>>(30, std::vector<double>(20, 0.0))));
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinearShape::InelasticDefgradLinearShape(
    const double GrowthFac, const double ReferenceValue)
    : GrowthFac_(GrowthFac), ReferenceValue_(ReferenceValue)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradLinearShape::EvaluateLinearGrowth(const double Value) const
{
  // calculate and return the linear growth factor
  return GrowthFac_ * (Value - ReferenceValue_);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolynomialShape::InelasticDefgradPolynomialShape(
    std::vector<double> PolyCoeffs, const double Xmin, const double Xmax)
    : PolyCoeffs_(std::move(PolyCoeffs)), Xmin_(Xmin), Xmax_(Xmax)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradPolynomialShape::ComputePolynomial(const double X)
{
  // initialize the variable for the evaluation of the polynomial
  double Polynom(0.0);

  // compute polynomial
  for (unsigned i = 0; i < PolyCoeffs_.size(); ++i) Polynom += PolyCoeffs_[i] * std::pow(X, i);

  return Polynom;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradPolynomialShape::ComputePolynomialDerivative(double X)
{
  // initialize the variable for the derivative of the polynomial
  double PolynomDerivative(0.0);

  // compute first derivative of polynomial
  for (unsigned i = 1; i < PolyCoeffs_.size(); ++i)
    PolynomDerivative += i * PolyCoeffs_[i] * std::pow(X, i - 1);

  return PolynomDerivative;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolynomialShape::CheckPolynomialBounds(const double X) const
{
  // safety check for validity of polynomial
  if ((X < Xmin_) or (X > Xmax_))
  {
    std::cout << "WARNING: Polynomial is evaluated outside its range of validity!" << std::endl;
    std::cout << "Evaluation at: " << X << " Lower bound is " << Xmin_ << " Upper bound is "
              << Xmax_ << std::endl;
  }
}