/*----------------------------------------------------------------------*/
/*! \file
\brief evaluation of inelastic deformation gradients and their derivatives

\level 3

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/

/* headers */
#include "inelastic_defgrad_factors.H"
#include "matpar_bundle.H"
#include "material_service.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/voigt_notation.H"

#include "multiplicative_split_defgrad_elasthyper.H"

/*--------------------------------------------------------------------*
 | constructors of parameter classes                    schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradScalar::InelasticDefgradScalar(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      Scalar1_(matdata->GetInt("SCALAR1")),
      Scalar1refconc_(matdata->GetDouble("SCALAR1_RefConc")),
      Matid_(matdata->GetInt("MATID")),
      Cmax_(-1.0)
{
  // safety checks
  // in case not all scatra dofs are transported scalars, the last scatra dof is a potential and can
  // not be treated as a concentration but it is treated like that in so3_scatra_evaluate.cpp in the
  // PreEvaluate method!
  if (Scalar1_ != 1) dserror("At the moment it is only possible that SCALAR1 induces growth");
  if (Scalar1refconc_ < 0.0) dserror("The reference concentration of SCALAR1 can't be negative");

  // check correct masslin type
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  if (DRT::INPUT::IntegralValue<INPAR::STR::MassLin>(sdyn, "MASSLIN") != INPAR::STR::ml_none)
    dserror(
        "If you use the material 'InelasticDefgradLinScalarIso' please set 'MASSLIN' in the "
        "STRUCTURAL DYNAMIC Section to 'None', or feel free to implement other possibility!");

  // Check if the material specified by user with MATID is one of allowed materials which are
  // currently electrode, elchmat and scatra
  if (Matid_ > 0)
  {
    // retrieve problem instance to read from
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    // retrieve validated input line of material ID in question
    Teuchos::RCP<MAT::PAR::Material> curmat =
        DRT::Problem::Instance(probinst)->Materials()->ById(Matid_);
    switch (curmat->Type())
    {
      case INPAR::MAT::m_electrode:
      {
        // Get C_max of electrode material
        Cmax_ = curmat->GetDouble("C_MAX");
        break;
      }
      case INPAR::MAT::m_elchmat:
      {
        Cmax_ = 1.0;
        break;
      }
      case INPAR::MAT::m_scatra:
      {
        Cmax_ = 1.0;
        break;
      }
      default:
        dserror(
            "The material you have specified by MATID has to be an electrode, "
            "elchmat or scatra material!");
    }
  }
  else
  {
    dserror("You have to enter a valid MATID for the corresponding scatra material!");
  }

  return;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradLinScalar::InelasticDefgradLinScalar(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : InelasticDefgradScalar(matdata), Scalar1growthfac_(matdata->GetDouble("SCALAR1_GrowthFac"))
{
  if (Scalar1growthfac_ < 0.0)
    dserror("The influence of scalar field SCALAR1 to growth can't be negativ");

  return;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradPolyScalar::InelasticDefgradPolyScalar(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : InelasticDefgradScalar(matdata),
      Polyparanum_(matdata->GetInt("POLY_PARA_NUM")),
      Polyparams_(*matdata->Get<std::vector<double>>("POLY_PARAMS")),
      Xmin_(matdata->GetDouble("X_min")),
      Xmax_(matdata->GetDouble("X_max")),
      MATPolynomReference_(0.0)
{
  // safety checks
  // in case not all scatra dofs are transported scalars, the last scatra dof is a potential and can
  // not be treated as a concentration but it is treated like that in so3_scatra_evaluate.cpp in the
  // PreEvaluate method!
  if (Polyparams_.size() != Polyparanum_)
    dserror(
        "Number of coefficients POLY_PARA_NUM you entered in input file has to match the size of "
        "coefficient vector POLY_PARAMS");

  // evaluate polynomial at reference X
  MAT::InelasticDefgradPolynomial Polynomial;
  MATPolynomReference_ = Polynomial.EvaluateMatPolynomial(
      Scalar1(), Polyparanum_, Polyparams_, Xmin_, Xmax_, Scalar1refconc(), Cmax(), 1.0);
  return;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradLinScalarAniso::InelasticDefgradLinScalarAniso(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : InelasticDefgradLinScalar(matdata),
      growthdir_(*matdata->Get<std::vector<double>>("GrowthDirection"))
{
  return;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradPolyScalarAniso::InelasticDefgradPolyScalarAniso(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : InelasticDefgradPolyScalar(matdata),
      growthdir_(*matdata->Get<std::vector<double>>("GrowthDirection"))
{
  return;
}

/*--------------------------------------------------------------------*
 | constructors of Growth direction                                   |
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
  const double growthdirvecnorm = std::sqrt(
      pow(growthdirection[0], 2.0) + pow(growthdirection[1], 2.0) + pow(growthdirection[2], 2.0));
  const double invquadrgrowthdirvecnorm = 1.0 / (growthdirvecnorm * growthdirvecnorm);

  // loop over all rows and colomns to fill the matrix and scale it correctly on the fly
  for (unsigned i = 0; i < growthdirection.size(); ++i)
  {
    for (unsigned j = 0; j < growthdirection.size(); ++j)
    {
      growthdirmat_(i, j) = invquadrgrowthdirvecnorm * growthdirection[i] * growthdirection[j];
    }
  }
}

/*--------------------------------------------------------------------*
 | construct empty material                             schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradFactors::InelasticDefgradFactors() : params_(NULL) {}


/*--------------------------------------------------------------------*
 | construct material with specific material params     schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradFactors::InelasticDefgradFactors(MAT::PAR::Parameter* params) : params_(params)
{
}


/*--------------------------------------------------------------------*
 | factory                                              schmidt 03/18 |
 *--------------------------------------------------------------------*/
Teuchos::RCP<MAT::InelasticDefgradFactors> MAT::InelasticDefgradFactors::Factory(int matnum)
{
  // for the sake of safety
  if (DRT::Problem::Instance()->Materials() == Teuchos::null)
    dserror("List of materials cannot be accessed in the global problem instance.");

  // another safety check
  if (DRT::Problem::Instance()->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matnum);

  switch (curmat->Type())
  {
    case INPAR::MAT::mfi_lin_scalar_iso:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::InelasticDefgradLinScalarIso(curmat));
      MAT::PAR::InelasticDefgradLinScalarIso* params =
          dynamic_cast<MAT::PAR::InelasticDefgradLinScalarIso*>(curmat->Parameter());
      return Teuchos::rcp(new InelasticDefgradLinScalarIso(params));
    }
    case INPAR::MAT::mfi_lin_scalar_aniso:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::InelasticDefgradLinScalarAniso(curmat));
      MAT::PAR::InelasticDefgradLinScalarAniso* params =
          dynamic_cast<MAT::PAR::InelasticDefgradLinScalarAniso*>(curmat->Parameter());
      return Teuchos::rcp(new InelasticDefgradLinScalarAniso(params));
    }
    case INPAR::MAT::mfi_poly_scalar_aniso:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::InelasticDefgradPolyScalarAniso(curmat));
      MAT::PAR::InelasticDefgradPolyScalarAniso* params =
          dynamic_cast<MAT::PAR::InelasticDefgradPolyScalarAniso*>(curmat->Parameter());
      return Teuchos::rcp(new InelasticDefgradPolyScalarAniso(params));
    }
    case INPAR::MAT::mfi_poly_scalar_iso:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::InelasticDefgradPolyScalarIso(curmat));
      MAT::PAR::InelasticDefgradPolyScalarIso* params =
          dynamic_cast<MAT::PAR::InelasticDefgradPolyScalarIso*>(curmat->Parameter());
      return Teuchos::rcp(new InelasticDefgradPolyScalarIso(params));
    }
    default:
      dserror("cannot deal with type %d", curmat->Type());
  }
  // dummy return
  return Teuchos::null;
}


/*--------------------------------------------------------------------*
 | construct empty material                             civaner 07/19 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinear::InelasticDefgradLinear() : InelasticDefgradFactors() {}


/*--------------------------------------------------------------------*
 | construct material with specific material params     schmidt 07/19 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinear::InelasticDefgradLinear(MAT::PAR::Parameter* params)
    : InelasticDefgradFactors(params)
{
}

/*--------------------------------------------------------------------*
 | evaluate scaling factor                              civaner 07/19 |
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradLinear::EvaluateDeltaGrowth(const double Growthfac,
    const double DetJacobian, const double Referenceconc, const double concentration,
    const double Cmax)
{
  // compute delta intercalation fraction between reference and current state
  const double deltaX = (concentration - Referenceconc) / Cmax;
  // calculate growth
  const double scalingfactor = Growthfac * DetJacobian * deltaX;
  return scalingfactor;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticSource MAT::InelasticDefgradLinear::GetInelasticSource()
{
  return PAR::InelasticSource::inelastic_concentration;
};

/*--------------------------------------------------------------------*
 | standard constructor                                 civaner 07/19 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolynomial::InelasticDefgradPolynomial() {}

/*--------------------------------------------------------------------*
 | construct material with specific material params     schmidt 07/19 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolynomial::InelasticDefgradPolynomial(MAT::PAR::Parameter* params)
    : InelasticDefgradFactors(params)
{
}


/*--------------------------------------------------------------------*
 | evaluate polynomial                                   civaner 07/19 |
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradPolynomial::EvaluateMatPolynomial(const double Scalar,
    const int NumberPolyparams, const std::vector<double> ParamsPoly, const double X_min,
    const double X_max, double concentration, const double Cmax, double DetJacobian)
{
  // initialize the variables for polynomial
  double MATPolynom(0.0);
  // safety check for validity of polynomial
  const double X = (concentration * DetJacobian) / Cmax;
  // print warning to screen if prescribed interval of validity for polynomial
  // calculation model is not satisfied
  if ((X < X_min) or (X > X_max))
  {
    std::cout
        << "WARNING: intercalation fraction X = c / c_max * detJacobian is violating prescribed "
           "bounds of polynomial calculation model. Calculated values might therefore not be "
           "reasonable!"
        << std::endl;
    std::cout << "Current X: " << X << " Lower bound is " << X_min << " Upper bound is " << X_max
              << std::endl;
  }
  // compute polynomial
  for (int i = 0; i < NumberPolyparams; ++i) MATPolynom += ParamsPoly[i] * pow(X, i);
  return MATPolynom;
}


/*--------------------------------------------------------------------*
 | evaluate first derivative of polynomial              civaner 07/19 |
 *--------------------------------------------------------------------*/
double MAT::InelasticDefgradPolynomial::EvaluateMatPolynomialDerivative(const double Scalar,
    const int NumberPolyparams, const std::vector<double> ParamsPoly, double concentration,
    const double Cmax, double DetJacobian)
{
  // initialize the variables for polynomial
  double DerivativeMATPolynom(0.0);
  // intercalation fraction
  const double X = (concentration * DetJacobian) / Cmax;
  // compute first derivative of polynomial
  for (int i = 1; i < NumberPolyparams; ++i)
    DerivativeMATPolynom += i * ParamsPoly[i] * pow(X, i - 1);
  return DerivativeMATPolynom;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticSource MAT::InelasticDefgradPolynomial::GetInelasticSource()
{
  return PAR::InelasticSource::inelastic_concentration;
};


/*--------------------------------------------------------------------*
 | construct empty material                             schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarIso::InelasticDefgradLinScalarIso()
    : InelasticDefgradLinear(), gp_(-1), concentrations_(NULL)
{
}


/*--------------------------------------------------------------------*
 | construct material with specific material params     schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarIso::InelasticDefgradLinScalarIso(MAT::PAR::Parameter* params)
    : InelasticDefgradLinear(params), gp_(-1), concentrations_(NULL)
{
}


/*--------------------------------------------------------------------*
 | evaluate the inverse of the inelastic deformation                  |
 | gradient                                             schmidt 03/18 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarIso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const double Sc1growthfac = Parameter()->Scalar1growthfac();
  const double Sc1refconc = Parameter()->Scalar1refconc();
  const double Cmax = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();

  // get scaling factor
  const double deltagrowth = EvaluateDeltaGrowth(
      Sc1growthfac, detjacobian, Sc1refconc, concentrations_->at(gp_).at(Sc1 - 1), Cmax);
  const double isoinelasticdefo = pow(1.0 + deltagrowth, (1.0 / 3.0));

  // calculate inverse inelastic deformation gradient (FinM is modeled, such that the volume change
  // is a linear function of the scalar (mapped to reference frame) that causes it)
  for (int i = 0; i < 3; ++i) iFinM(i, i) = 1.0 / isoinelasticdefo;

  return;
}


/*--------------------------------------------------------------------*
 | evaluate additional contribution to cmat             schmidt 03/18 |
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
  const double Sc1 = Parameter()->Scalar1();
  const double Sc1growthfac = Parameter()->Scalar1growthfac();
  const double Sc1refconc = Parameter()->Scalar1refconc();
  const double Cmax = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();

  // get delta growth
  const double deltagrowth = EvaluateDeltaGrowth(
      Sc1growthfac, detjacobian, Sc1refconc, concentrations_->at(gp_).at(Sc1 - 1), Cmax);
  // evaluate scaling factor
  const double scalefac = -deltagrowth / (6.0 * pow(1.0 + deltagrowth, 4.0 / 3.0));

  // diFinjdC = - growthfac*det(F)*(c - c_{ref})/(6*[1 + growthfac*det(F)*(c - c_{ref})]^(4/3)) * I
  diFinjdC.MultiplyNT(scalefac, id9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.MultiplyNN(2.0, dSdiFinj, diFinjdC, 1.0);

  return;
}


/*--------------------------------------------------------------------*
 | evaluate the off-diagonal contribution to the                      |
 | stiffness matrix (for monolithic calculation)        schmidt 03/18 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarIso::EvaluateODStiffMat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc)
{
  static LINALG::Matrix<9, 1> id9x1(true);
  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const double Sc1growthfac = Parameter()->Scalar1growthfac();
  const double Sc1refconc = Parameter()->Scalar1refconc();
  const double Cmax = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();

  // get delta growth
  const double deltagrowth = EvaluateDeltaGrowth(
      Sc1growthfac, detjacobian, Sc1refconc, concentrations_->at(gp_).at(Sc1 - 1), Cmax);

  // calculate scalefac
  const double temp = Sc1growthfac * detjacobian / Cmax;
  const double scalefac = -temp / (3.0 * pow(1.0 + deltagrowth, 4.0 / 3.0));

  // dstressdc = dSdiFinj : diFinjdc
  // diFinjdc = - growthfac*det(F)/(3*[1 + growthfac*det(F)*(c-c_{ref})]^(4/3)) I
  dstressdc.MultiplyNN(scalefac, dSdiFinj, id9x1, 1.0);

  return;
}


/*--------------------------------------------------------------------*
 | pre evaluate                                         schmidt 03/18 |
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

  return;
}


/*--------------------------------------------------------------------*
 | construct empty material                             schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarAniso::InelasticDefgradLinScalarAniso()
    : InelasticDefgradLinear(), gp_(-1), concentrations_(NULL)
{
}


/*--------------------------------------------------------------------*
 | construct material with specific material params     schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarAniso::InelasticDefgradLinScalarAniso(MAT::PAR::Parameter* params)
    : InelasticDefgradLinear(params), gp_(-1), concentrations_(NULL)
{
}


/*--------------------------------------------------------------------*
 | evaluate the inverse of the inelastic deformation                  |
 | gradient                                             schmidt 03/18 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarAniso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  // init and clear variable
  static LINALG::Matrix<3, 3> FinM(true);
  FinM.Clear();

  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const double Sc1growthfac = Parameter()->Scalar1growthfac();
  const double Sc1refconc = Parameter()->Scalar1refconc();
  const double detjacobian = defgrad->Determinant();
  const double Cmax = Parameter()->Cmax();

  // get delta growth
  const double deltagrowth = EvaluateDeltaGrowth(
      Sc1growthfac, detjacobian, Sc1refconc, concentrations_->at(gp_).at(Sc1 - 1), Cmax);

  // calculate inelastic deformation gradient and its inverse
  for (int i = 0; i < 3; ++i) FinM(i, i) = 1.0;

  // finalize inelastic deformation gradient matrix (FinM is calculated, such that the volume change
  // is a linear function of the scalar (mapped to reference frame) that causes it)
  FinM.Update(deltagrowth, Parameter()->Growthdirmat(), 1.0);
  // calculate inverse of inelastic deformation gradient matrix
  iFinM.Invert(FinM);

  return;
}


/*--------------------------------------------------------------------*
 | evaluate additional contribution to cmat             schmidt 03/18 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarAniso::EvaluateAdditionalCmat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 1>& iCV, const LINALG::Matrix<6, 9>& dSdiFinj,
    LINALG::Matrix<6, 6>& cmatadd)
{
  static LINALG::Matrix<3, 3> temp(true);
  static LINALG::Matrix<3, 3> iFinjGiFinj(true);
  static LINALG::Matrix<9, 1> diFinjddetJ9x1(true);
  static LINALG::Matrix<9, 6> diFinjdC(true);

  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const double Sc1growthfac = Parameter()->Scalar1growthfac();
  const double Sc1refconc = Parameter()->Scalar1refconc();
  const double Cmax = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();

  // get delta growth
  const double deltagrowth = EvaluateDeltaGrowth(
      Sc1growthfac, detjacobian, Sc1refconc, concentrations_->at(gp_).at(Sc1 - 1), Cmax);

  const double scalefac = -deltagrowth / 2.0;

  // calculate F_{in,j}^{-1} . G . F_{in,j}^{-1} with F_{in,j}, the j-th factor of F_{in}
  temp.MultiplyNN(1.0, iFinjM, Parameter()->Growthdirmat(), 0.0);
  iFinjGiFinj.MultiplyNN(1.0, temp, iFinjM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(iFinjGiFinj, diFinjddetJ9x1);

  // diFinjdC = - growthfac*det(F)*(c-c_{ref})/2.0 (F_{in,j}^{-1} . G . F_{in,j}^{-1}) \otimes
  // C^{-1}
  diFinjdC.MultiplyNT(scalefac, diFinjddetJ9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.MultiplyNN(2.0, dSdiFinj, diFinjdC, 1.0);

  return;
}


/*--------------------------------------------------------------------*
 | evaluate the off-diagonal contribution to the                      |
 | stiffness matrix (for monolithic calculation)        schmidt 03/18 |
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
  const double Sc1growthfac = Parameter()->Scalar1growthfac();
  const double Cmax = Parameter()->Cmax();
  // determinant of jacobian
  const double detjacobian = defgrad->Determinant();
  // scalefac
  const double scalefac = -Sc1growthfac * detjacobian / Cmax;

  // diFinjdc = - growthfac*det(F) F_{in,j}^{-1} . G . F_{in,j}^{-1} with F_{in,j}, the j-th factor
  // of F_{in}
  tmp.MultiplyNN(1.0, iFinjM, Parameter()->Growthdirmat(), 0.0);
  diFinjdcM.MultiplyNN(scalefac, tmp, iFinjM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(diFinjdcM, diFinjdc9x1);

  // dstressdc = dSdiFinj : diFinjdc
  dstressdc.MultiplyNN(1.0, dSdiFinj, diFinjdc9x1, 1.0);

  return;
}


/*--------------------------------------------------------------------*
 | pre evaluate                                         schmidt 03/18 |
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

  return;
}

/*--------------------------------------------------------------------*
 | construct empty material                             civaner 03/19 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolyScalarIso::InelasticDefgradPolyScalarIso()
    : InelasticDefgradPolynomial(), gp_(-1), concentrations_(NULL)
{
}


/*--------------------------------------------------------------------*
 | construct material with specific material params     civaner 03/19 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolyScalarIso::InelasticDefgradPolyScalarIso(MAT::PAR::Parameter* params)
    : InelasticDefgradPolynomial(params), gp_(-1), concentrations_(NULL)
{
}

/*--------------------------------------------------------------------*
 | evaluate the inverse of the inelastic deformation                  |
 | gradient                                             civaner 03/19 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyScalarIso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const int Polyparanum = Parameter()->Polyparanum();
  const std::vector<double> Polyparams = Parameter()->Polyparams();
  const double Xmin = Parameter()->Xmin();
  const double Xmax = Parameter()->Xmax();
  const double Cmax = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();
  const double MATPolynomReference = Parameter()->MATPolynomReference();

  // get polynomial
  const double MATPolynom = EvaluateMatPolynomial(Sc1, Polyparanum, Polyparams, Xmin, Xmax,
      concentrations_->at(gp_).at(Sc1 - 1), Cmax, detjacobian);
  // calculate growth
  const double isoinelasticdefo = pow((1.0 + MATPolynom) / (1 + MATPolynomReference), (1.0 / 3.0));
  // calculate inverse inelastic deformation gradient (FinM is modeled, such that the volume change
  // is a polynomial function of the scalar (mapped to reference frame) that causes it)
  for (int i = 0; i < 3; ++i) iFinM(i, i) = 1.0 / isoinelasticdefo;
  return;
}


/*--------------------------------------------------------------------*
 | evaluate additional contribution to cmat             civaner 03/19 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyScalarIso::EvaluateAdditionalCmat(
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
  const double Sc1 = Parameter()->Scalar1();
  const int Polyparanum = Parameter()->Polyparanum();
  const std::vector<double> Polyparams = Parameter()->Polyparams();
  const double Xmin = Parameter()->Xmin();
  const double Xmax = Parameter()->Xmax();
  const double Cmax = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();
  const double MATPolynomReference = Parameter()->MATPolynomReference();

  // calculate scaling factor
  // get polynomial
  const double MATPolynom = EvaluateMatPolynomial(Sc1, Polyparanum, Polyparams, Xmin, Xmax,
      concentrations_->at(gp_).at(Sc1 - 1), Cmax, detjacobian);
  const double DerMATPolynom = EvaluateMatPolynomialDerivative(
      Sc1, Polyparanum, Polyparams, concentrations_->at(gp_).at(Sc1 - 1), Cmax, detjacobian);
  // calculate growth
  const double scalefac = -1.0 / 6.0 * (concentrations_->at(gp_).at(Sc1 - 1)) / Cmax * detjacobian *
                          pow(1.0 + MATPolynom, -4.0 / 3.0) * DerMATPolynom *
                          pow(1.0 + MATPolynomReference, 1.0 / 3.0);

  // diFinjdC = - growthfac*det(F)*(c - c_{ref})/(6*[1 + growthfac*det(F)*(c - c_{ref})]^(4/3)) * I
  diFinjdC.MultiplyNT(scalefac, id9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.MultiplyNN(2.0, dSdiFinj, diFinjdC, 1.0);

  return;
}


/*--------------------------------------------------------------------*
 | evaluate the off-diagonal contribution to the                      |
 | stiffness matrix (for monolithic calculation)        civaner 03/19 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyScalarIso::EvaluateODStiffMat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc)
{
  static LINALG::Matrix<9, 1> id9x1(true);
  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const int Polyparanum = Parameter()->Polyparanum();
  const std::vector<double> Polyparams = Parameter()->Polyparams();
  const double Xmin = Parameter()->Xmin();
  const double Xmax = Parameter()->Xmax();
  const double Cmax = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();
  const double MATPolynomReference = Parameter()->MATPolynomReference();

  // calculate scalefac
  // get polynomial
  const double MATPolynom = EvaluateMatPolynomial(Sc1, Polyparanum, Polyparams, Xmin, Xmax,
      concentrations_->at(gp_).at(Sc1 - 1), Cmax, detjacobian);
  const double DerMATPolynom = EvaluateMatPolynomialDerivative(
      Sc1, Polyparanum, Polyparams, concentrations_->at(gp_).at(Sc1 - 1), Cmax, detjacobian);
  // calculate scalefac
  const double scalefac = -1.0 / 3.0 * pow(1.0 + MATPolynom, -4.0 / 3.0) *
                          pow(1.0 + MATPolynomReference, 1.0 / 3.0) * DerMATPolynom * detjacobian /
                          Cmax;
  // dstressdc = dSdiFinj : diFinjdc
  // diFinjdc = - growthfac*det(F)/(3*[1 + growthfac*det(F)*(c-c_{ref})]^(4/3)) I
  dstressdc.MultiplyNN(scalefac, dSdiFinj, id9x1, 1.0);

  return;
}


/*--------------------------------------------------------------------*
 | pre evaluate                                         civaner 03/19 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyScalarIso::PreEvaluate(Teuchos::ParameterList& params)
{
  // get Gauss point number
  gp_ = params.get<int>("gp", -1);
  if (gp_ == -1) dserror("No Gauss point number provided in material.");

  // set pointer to vector of gp_conc, only if gp is 0, because this is the first gp
  if (gp_ == 0)
    concentrations_ = params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc",
        Teuchos::rcp(new std::vector<std::vector<double>>(30, std::vector<double>(20, 0.0))));

  return;
}



/*--------------------------------------------------------------------*
 | construct empty material                             civaner 03/19 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolyScalarAniso::InelasticDefgradPolyScalarAniso()
    : InelasticDefgradPolynomial(), gp_(-1), concentrations_(NULL)
{
}


/*--------------------------------------------------------------------*
 | construct material with specific material params     civaner 03/19 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradPolyScalarAniso::InelasticDefgradPolyScalarAniso(MAT::PAR::Parameter* params)
    : InelasticDefgradPolynomial(params), gp_(-1), concentrations_(NULL)
{
}


/*--------------------------------------------------------------------*
 | evaluate the inverse of the inelastic deformation                  |
 | gradient                                             civaner 03/19 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyScalarAniso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* defgrad, LINALG::Matrix<3, 3>& iFinM)
{
  // init and clear variable
  static LINALG::Matrix<3, 3> FinM(true);
  FinM.Clear();

  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const double detjacobian = defgrad->Determinant();
  const int Polyparanum = Parameter()->Polyparanum();
  const std::vector<double> Polyparams = Parameter()->Polyparams();
  const double Xmin = Parameter()->Xmin();
  const double Xmax = Parameter()->Xmax();
  const double Cmax = Parameter()->Cmax();
  const double MATPolynomReference = Parameter()->MATPolynomReference();

  // get polynomials
  const double MATPolynom = EvaluateMatPolynomial(Sc1, Polyparanum, Polyparams, Xmin, Xmax,
      concentrations_->at(gp_).at(Sc1 - 1), Cmax, detjacobian);
  // calculate growth
  const double deltagrowth = (MATPolynom - MATPolynomReference) / (MATPolynomReference + 1.0);
  // calculate inelastic deformation gradient and its inverse
  for (int i = 0; i < 3; ++i) FinM(i, i) = 1.0;

  // finalize inelastic deformation gradient matrix (FinM is calculated, such that the volume change
  // is a linear function of the scalar (mapped to reference frame) that causes it)
  FinM.Update(deltagrowth, Parameter()->Growthdirmat(), 1.0);
  // calculate inverse of inelastic deformation gradient matrix
  iFinM.Invert(FinM);

  return;
}


/*--------------------------------------------------------------------*
 | evaluate additional contribution to cmat             civaner 03/19 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyScalarAniso::EvaluateAdditionalCmat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 1>& iCV, const LINALG::Matrix<6, 9>& dSdiFinj,
    LINALG::Matrix<6, 6>& cmatadd)
{
  static LINALG::Matrix<3, 3> temp(true);
  static LINALG::Matrix<3, 3> iFinjGiFinj(true);
  static LINALG::Matrix<9, 1> diFinjddetJ9x1(true);
  static LINALG::Matrix<9, 6> diFinjdC(true);

  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const int Polyparanum = Parameter()->Polyparanum();
  const std::vector<double> Polyparams = Parameter()->Polyparams();
  const double Cmax = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();
  const double MATPolynomReference = Parameter()->MATPolynomReference();

  // get polynomial
  const double DerMATPolynom = EvaluateMatPolynomialDerivative(
      Sc1, Polyparanum, Polyparams, concentrations_->at(gp_).at(Sc1 - 1), Cmax, detjacobian);

  const double scalefac = -0.5 * detjacobian * concentrations_->at(gp_).at(Sc1 - 1) / Cmax *
                          DerMATPolynom / (MATPolynomReference + 1.0);

  // calculate F_{in,j}^{-1} . G . F_{in,j}^{-1} with F_{in,j}, the j-th factor of F_{in}
  temp.MultiplyNN(1.0, iFinjM, Parameter()->Growthdirmat(), 0.0);
  iFinjGiFinj.MultiplyNN(1.0, temp, iFinjM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(iFinjGiFinj, diFinjddetJ9x1);

  // diFinjdC = - growthfac*det(F)*(c-c_{ref})/2.0 (F_{in,j}^{-1} . G . F_{in,j}^{-1}) \otimes
  // C^{-1}
  diFinjdC.MultiplyNT(scalefac, diFinjddetJ9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.MultiplyNN(2.0, dSdiFinj, diFinjdC, 1.0);

  return;
}


/*--------------------------------------------------------------------*
 | evaluate the off-diagonal contribution to the                      |
 | stiffness matrix (for monolithic calculation)        civaner 03/19 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyScalarAniso::EvaluateODStiffMat(
    const LINALG::Matrix<3, 3>* const defgrad, const LINALG::Matrix<3, 3>& iFinjM,
    const LINALG::Matrix<6, 9>& dSdiFinj, LINALG::Matrix<6, 1>& dstressdc)
{
  // static variables
  static LINALG::Matrix<3, 3> tmp(true);
  static LINALG::Matrix<3, 3> diFinjdcM(true);
  static LINALG::Matrix<9, 1> diFinjdc9x1(true);

  // get parameters
  const int Polyparanum = Parameter()->Polyparanum();
  const double Sc1 = Parameter()->Scalar1();
  const std::vector<double> Polyparams = Parameter()->Polyparams();
  const double Cmax = Parameter()->Cmax();
  const double detjacobian = defgrad->Determinant();
  const double MATPolynomReference = Parameter()->MATPolynomReference();

  // get polynomials
  const double DerMATPolynom = EvaluateMatPolynomialDerivative(
      Sc1, Polyparanum, Polyparams, concentrations_->at(gp_).at(Sc1 - 1), Cmax, detjacobian);

  // scalefac
  const double scalefac = -detjacobian * DerMATPolynom / Cmax / (MATPolynomReference + 1.0);

  // diFinjdc = - growthfac*det(F) F_{in,j}^{-1} . G . F_{in,j}^{-1} with F_{in,j}, the j-th factor
  // of F_{in}
  tmp.MultiplyNN(1.0, iFinjM, Parameter()->Growthdirmat(), 0.0);
  diFinjdcM.MultiplyNN(scalefac, tmp, iFinjM, 0.0);
  UTILS::VOIGT::Matrix3x3to9x1(diFinjdcM, diFinjdc9x1);

  // dstressdc = dSdiFinj : diFinjdc
  dstressdc.MultiplyNN(1.0, dSdiFinj, diFinjdc9x1, 1.0);
  return;
}


/*--------------------------------------------------------------------*
 | pre evaluate                                         civaner 03/19 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradPolyScalarAniso::PreEvaluate(Teuchos::ParameterList& params)
{
  // get Gauss point number
  gp_ = params.get<int>("gp", -1);
  if (gp_ == -1) dserror("No Gauss point number provided in material.");

  // set pointer to vector of gp_conc, only if gp is 0, because this is the first gp
  if (gp_ == 0)
    concentrations_ = params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc",
        Teuchos::rcp(new std::vector<std::vector<double>>(30, std::vector<double>(20, 0.0))));

  return;
}
