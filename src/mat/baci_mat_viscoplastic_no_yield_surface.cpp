/*----------------------------------------------------------------------*/
/*! \file
\brief An elastic visco-plastic material law without yield surface

\level 2
 *----------------------------------------------------------------------*/

#include "baci_mat_viscoplastic_no_yield_surface.H"

#include "baci_lib_globalproblem.H"
#include "baci_lib_voigt_notation.H"
#include "baci_linalg_utils_densematrix_eigen.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_service.H"
#include "baci_so3_hex8.H"
#include "baci_utils_fad.H"

#include <vector>

// struct definition
struct MAT::PreCalculatedTerms
{
  double equ_tens_stress_const;
  double flow_res_const;
  double equ_tens_stress_flow_res_ratio;
  double flow_res_sat_deviation;
};

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::ViscoPlasticNoYieldSurface::ViscoPlasticNoYieldSurface(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      density_(matdata->GetDouble("DENS")),
      nue_(matdata->GetDouble("NUE")),
      young_(matdata->GetDouble("YOUNG")),
      temperature_(matdata->GetDouble("TEMPERATURE")),
      pre_exp_fac_(matdata->GetDouble("PRE_EXP_FAC")),
      activation_energy_(matdata->GetDouble("ACTIVATION_ENERGY")),
      gas_constant_(matdata->GetDouble("GAS_CONSTANT")),
      strain_rate_sensitivity_(matdata->GetDouble("STRAIN_RATE_SENS")),
      init_flow_res_(matdata->GetDouble("INIT_FLOW_RES")),
      flow_res_exp_(matdata->GetDouble("FLOW_RES_EXP")),
      flow_res_pre_fac_(matdata->GetDouble("FLOW_RES_PRE_FAC")),
      flow_res_sat_fac_(matdata->GetDouble("FLOW_RES_SAT_FAC")),
      flow_res_sat_exp_(matdata->GetDouble("FLOW_RES_SAT_EXP"))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ViscoPlasticNoYieldSurface::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ViscoPlasticNoYieldSurface(this));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ViscoPlasticNoYieldSurfaceType MAT::ViscoPlasticNoYieldSurfaceType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::ViscoPlasticNoYieldSurfaceType::Create(const std::vector<char>& data)
{
  auto* visco_plastic_no_yield_surface = new MAT::ViscoPlasticNoYieldSurface();
  visco_plastic_no_yield_surface->Unpack(data);
  return visco_plastic_no_yield_surface;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ViscoPlasticNoYieldSurface::ViscoPlasticNoYieldSurface()
    : params_(nullptr),
      last_plastic_defgrd_inverse_(Teuchos::null),
      last_flowres_isotropic_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ViscoPlasticNoYieldSurface::ViscoPlasticNoYieldSurface(
    MAT::PAR::ViscoPlasticNoYieldSurface* params)
    : params_(params)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoPlasticNoYieldSurface::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  // in case we are in post-process mode
  if (Parameter() != nullptr) matid = Parameter()->Id();
  AddtoPack(data, matid);

  // pack history data
  AddtoPack<3, 3>(data, last_plastic_defgrd_inverse_);
  AddtoPack(data, last_flowres_isotropic_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoPlasticNoYieldSurface::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<MAT::PAR::ViscoPlasticNoYieldSurface*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  // history data
  ExtractfromPack<3, 3>(position, data, last_plastic_defgrd_inverse_);
  ExtractfromPack(position, data, last_flowres_isotropic_);

  // no need to unpack this, just resize the data members
  current_flowres_isotropic_.resize(last_flowres_isotropic_.size(), 0.0);
  current_plastic_defgrd_inverse_.resize(last_plastic_defgrd_inverse_.size());

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
void MAT::ViscoPlasticNoYieldSurface::Setup(const int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // read initial flow resistance from line definition
  last_flowres_isotropic_.resize(numgp, params_->InitFlowRes());

  // initialize last inverse plastic deformation gradient as identity
  CORE::LINALG::Matrix<3, 3> id2(true);
  for (int i = 0; i < 3; ++i) id2(i, i) = 1.0;
  last_plastic_defgrd_inverse_.resize(numgp, id2);

  // initialize current history variables (values do not matter)
  current_flowres_isotropic_.resize(numgp, 0.0);
  current_plastic_defgrd_inverse_.resize(numgp, id2);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoPlasticNoYieldSurface::Update()
{
  // this update is done after the global newton loop is converged
  last_plastic_defgrd_inverse_ = current_plastic_defgrd_inverse_;
  last_flowres_isotropic_ = current_flowres_isotropic_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoPlasticNoYieldSurface::Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
    const CORE::LINALG::Matrix<6, 1>* strain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  // read input and history variables
  const double dt = params.get<double>("delta time");
  const CORE::LINALG::Matrix<3, 3>& last_iFv = last_plastic_defgrd_inverse_[gp];

  // trial (purely elastic) deformation gradient
  static CORE::LINALG::Matrix<3, 3> Fe_trial;
  Fe_trial.MultiplyNN(*defgrd, last_iFv);

  // define variables for eigenvalue analysis
  static CORE::LINALG::Matrix<3, 1> EigenvaluesFe_trial;
  static CORE::LINALG::Matrix<3, 3> EigenvectorsFe_trial;
  CalculateTrialElasticDefgradEigenvaluesAndEigenvectors(
      Fe_trial, EigenvaluesFe_trial, EigenvectorsFe_trial);

  // elastic trial rotation
  static CORE::LINALG::Matrix<3, 3> Re_trial;
  Re_trial = CalculateTrialElasticRotation(Fe_trial, EigenvectorsFe_trial, EigenvaluesFe_trial);

  // logarithmic strain
  static CORE::LINALG::Matrix<6, 1> Ee_trial_Vstrain;
  Ee_trial_Vstrain =
      CalculateLogElasticStrainInStrainLikeVoigtNotation(EigenvectorsFe_trial, EigenvaluesFe_trial);

  // setup elasticity tensor (stress-stress-like)
  SetupCmat(*cmat);

  // stress-like Voigt notation of stresses conjugated to logarithmic strains
  static CORE::LINALG::Matrix<6, 1> Me_trial_Vstress;
  // trial stress conjugate to logarithmic strain
  Me_trial_Vstress.MultiplyNN(*cmat, Ee_trial_Vstrain);

  // mean normal pressure (p = 1/3 * trace(Me_trial))
  const double p = 1.0 / 3.0 * (Me_trial_Vstress(0) + Me_trial_Vstress(1) + Me_trial_Vstress(2));

  static CORE::LINALG::Matrix<3, 3> Me_trial_dev;
  Me_trial_dev = CalculateDeviatoricTrialStresses(Me_trial_Vstress, p);
  const double Me_trial_eqv = CalculateTrialEquivalentStress(Me_trial_dev);

  // solution vector with variables x(0) = Me_eqv;  x(1) = S
  static CORE::LINALG::Matrix<2, 1> x;
  x(0) = Me_trial_eqv;
  x(1) = last_flowres_isotropic_[gp];

  // default value of eta = 1 (no plastic deformation) if trial equivalent stress is zero and no
  // mapping is needed
  double eta = 1.0;

  // do mapping if trial equivalent stress is not zero
  if (Me_trial_eqv != 0.0)
  {
    LocalNewtonLoop(x, dt);
    // calculate ratio of actual elastic stress (first component of solution vector) to trial stress
    eta = x(0) / Me_trial_eqv;
  }

  // conjugate stress (Me = eta * Me_trial_dev - p * id2)
  static CORE::LINALG::Matrix<3, 3> Me;
  // identity matrix
  static CORE::LINALG::Matrix<3, 3> id2(true);
  for (unsigned i = 0; i < 3; ++i) id2(i, i) = 1.0;
  Me.Update(eta, Me_trial_dev, p, id2);

  static CORE::LINALG::Matrix<3, 3> PK2;
  PK2 = CalculateSecondPiolaKirchhoffStresses(defgrd, Re_trial, Me);

  // current inverse plastic deformation gradient
  static CORE::LINALG::Matrix<3, 3> current_iFv;
  current_iFv = CalculateUpdatedInverseViscousDefgrad(
      last_iFv, EigenvectorsFe_trial, EigenvaluesFe_trial, eta);

  // update history variables of current Gauss-point
  current_plastic_defgrd_inverse_[gp] = current_iFv;
  current_flowres_isotropic_[gp] = x(1);

  // transform stresses to stress-like Voigt notation
  UTILS::VOIGT::Stresses::MatrixToVector(PK2, *stress);

  auto cmatel = CalculateElasticStiffness(EigenvectorsFe_trial, EigenvaluesFe_trial);
  *cmat = MAT::PullBackFourTensor<3>(current_iFv, cmatel);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<3, 3>& MAT::ViscoPlasticNoYieldSurface::CalculateDeviatoricTrialStresses(
    const CORE::LINALG::Matrix<6, 1>& Me_trial_Vstress, const double p) const
{
  // deviatoric trial stress in stress-like Voigt notation
  static CORE::LINALG::Matrix<6, 1> Me_trial_dev_Vstress;
  for (unsigned i = 0; i < 3; ++i) Me_trial_dev_Vstress(i) = Me_trial_Vstress(i) - p;
  for (unsigned i = 3; i < 6; ++i) Me_trial_dev_Vstress(i) = Me_trial_Vstress(i);

  static CORE::LINALG::Matrix<3, 3> Me_trial_dev;
  // transform deviatoric equivalent stress to matrix notation
  UTILS::VOIGT::Stresses::VectorToMatrix(Me_trial_dev_Vstress, Me_trial_dev);

  return Me_trial_dev;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<3, 3>& MAT::ViscoPlasticNoYieldSurface::CalculateUpdatedInverseViscousDefgrad(
    const CORE::LINALG::Matrix<3, 3>& last_iFv, const CORE::LINALG::Matrix<3, 3>& eigen_vectors,
    const CORE::LINALG::Matrix<3, 1>& eigen_values, const double eta) const
{
  static CORE::LINALG::Matrix<3, 3> scaledEigenvalues(true);

  // cubic root of determinant
  const double cubicRoot_stretchTensorTrial_determinant =
      std::cbrt(eigen_values(0) * eigen_values(1) * eigen_values(2));

  for (unsigned i = 0; i < 3; ++i)
  {
    scaledEigenvalues(i, i) =
        std::pow((eigen_values(i) / cubicRoot_stretchTensorTrial_determinant), (eta - 1.0));
  }

  static CORE::LINALG::Matrix<3, 3> tmp3x3;
  static CORE::LINALG::Matrix<3, 3> scaled_iUeTrialDeviatoric;
  static CORE::LINALG::Matrix<3, 3> current_iFv;
  // calculate inverse scaled trial stretch tensor
  tmp3x3.MultiplyNN(eigen_vectors, scaledEigenvalues);
  scaled_iUeTrialDeviatoric.MultiplyNT(tmp3x3, eigen_vectors);
  current_iFv.MultiplyNN(last_iFv, scaled_iUeTrialDeviatoric);

  return current_iFv;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::ViscoPlasticNoYieldSurface::CalculateTrialEquivalentStress(
    const CORE::LINALG::Matrix<3, 3>& Me_trial_dev) const
{
  static CORE::LINALG::Matrix<3, 3> Me_trial_dev_squared;
  Me_trial_dev_squared.Multiply(Me_trial_dev, Me_trial_dev);

  // calculate trace
  double trace = 0.0;
  for (unsigned i = 0; i < 3; ++i)
  {
    trace += Me_trial_dev_squared(i, i);
  }

  // trial equivalent stress = sqrt(3/2 * trace(Me_trial_dev . Me_trial_dev))
  return std::sqrt(3.0 / 2.0 * trace);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoPlasticNoYieldSurface::CalculateTrialElasticDefgradEigenvaluesAndEigenvectors(
    const CORE::LINALG::Matrix<3, 3>& Fe_trial, CORE::LINALG::Matrix<3, 1>& eigen_values,
    CORE::LINALG::Matrix<3, 3>& eigen_vectors) const
{
  // as we only have a method to calculate the eigenvalues and eigenvectors from symmetric matrices,
  // we calculate the right Cauchy-Green tensor first and then make use of the fact that we get the
  // squared eigen values

  // trial right Cauchy-Green tensor
  static CORE::LINALG::Matrix<3, 3> Ce_trial;
  Ce_trial.MultiplyTN(Fe_trial, Fe_trial);

  // squared principal stretches
  static CORE::LINALG::Matrix<3, 3> SquaredEigenvalues;
  // eigenvalue analysis of trial deformation
  CORE::LINALG::SYEV(Ce_trial, SquaredEigenvalues, eigen_vectors);

  // principal stretches
  for (unsigned i = 0; i < 3; ++i)
  {
    eigen_values(i) = std::sqrt(SquaredEigenvalues(i, i));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>&
MAT::ViscoPlasticNoYieldSurface::CalculateElasticStiffness(
    const CORE::LINALG::Matrix<3, 3>& eigen_vectors,
    const CORE::LINALG::Matrix<3, 1>& eigen_values) const
{
  // init and clear elastic stiffness matrix
  static CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> Ce;
  Ce.Clear();
  const double eps(1.0e-12);

  // elastic parameters
  const double E = params_->Young();
  const double nue = params_->Nue();
  const double G = E / (2.0 * (1.0 + nue));
  const double K = E / (3.0 * (1.0 - 2.0 * nue));

  // pre-calculation of terms
  static CORE::LINALG::Matrix<3, 1> LogEigenValues;
  static CORE::LINALG::Matrix<3, 1> SquaredEigenValues;
  static CORE::LINALG::Matrix<3, 1> ithEigenVector;
  std::vector<CORE::LINALG::Matrix<3, 1>> EigenVectorsVec;
  for (unsigned i = 0; i < 3; ++i)
  {
    LogEigenValues(i) = std::log(eigen_values(i));
    SquaredEigenValues(i) = std::pow(eigen_values(i), 2.0);

    // extract ith-eigenvector and add it to the vector
    for (unsigned j = 0; j < 3; ++j) ithEigenVector(j) = eigen_vectors(j, i);
    EigenVectorsVec.push_back(ithEigenVector);
  }

  const double traceLogEigenvalues = LogEigenValues(0) + LogEigenValues(1) + LogEigenValues(2);

  // define outer products of eigenvectors
  static CORE::LINALG::Matrix<3, 3> temp;
  static CORE::LINALG::Matrix<6, 1> NaNaV, NbNbV, NaNbV, NbNaV;
  for (unsigned a = 0; a < 3; ++a)
  {
    temp.MultiplyNT(EigenVectorsVec.at(a), EigenVectorsVec.at(a));
    UTILS::VOIGT::Stresses::MatrixToVector(temp, NaNaV);

    const double Se_a = (2.0 * G * (LogEigenValues(a) - 1.0 / 3.0 * traceLogEigenvalues) +
                            K * traceLogEigenvalues) /
                        SquaredEigenValues(a);

    for (unsigned b = 0; b < 3; ++b)
    {
      temp.MultiplyNT(EigenVectorsVec.at(a), EigenVectorsVec.at(b));
      UTILS::VOIGT::Stresses::MatrixToVector(temp, NaNbV);
      temp.MultiplyNT(EigenVectorsVec.at(b), EigenVectorsVec.at(a));
      UTILS::VOIGT::Stresses::MatrixToVector(temp, NbNaV);
      temp.MultiplyNT(EigenVectorsVec.at(b), EigenVectorsVec.at(b));
      UTILS::VOIGT::Stresses::MatrixToVector(temp, NbNbV);

      const double Se_b = (2.0 * G * (LogEigenValues(b) - 1.0 / 3.0 * traceLogEigenvalues) +
                              K * traceLogEigenvalues) /
                          SquaredEigenValues(b);

      double Ce_term1_fac = -(2.0 / 3.0 * G - K) / (SquaredEigenValues(a) * SquaredEigenValues(b));
      double Ce_term2_fac(0.0);

      if (a == b)
      {
        Ce_term1_fac -=
            (2.0 * G * (LogEigenValues(a) - 1.0 / 3.0 * traceLogEigenvalues - 1.0 / 2.0) +
                K * traceLogEigenvalues) *
            2.0 / (SquaredEigenValues(a) * SquaredEigenValues(b));
      }
      else
      {
        // principle stretches are almost equal -> l'Hospital is used
        if (std::abs(SquaredEigenValues(a) - SquaredEigenValues(b)) < eps)
        {
          // add derivative of S_e,b w.r.t. lambda_e,b^2
          Ce_term2_fac = -(2.0 * G * (LogEigenValues(b) - 1.0 / 3.0 * (traceLogEigenvalues + 1.0)) +
                             K * (traceLogEigenvalues - 1.0 / 2.0)) *
                         2.0 / (SquaredEigenValues(b) * SquaredEigenValues(b));
          // subtract derivative of S_e,a w.r.t. lambda_e,bn^2
          Ce_term2_fac -= Ce_term1_fac;
        }
        // principle stretches differ -> normal formula can be used
        else
        {
          Ce_term2_fac = (Se_b - Se_a) / (SquaredEigenValues(b) - SquaredEigenValues(a));
        }
      }

      Ce.MultiplyNT(Ce_term1_fac, NaNaV, NbNbV, 1.0);
      Ce.MultiplyNT(Ce_term2_fac, NaNbV, NaNbV, 1.0);
      Ce.MultiplyNT(Ce_term2_fac, NaNbV, NbNaV, 1.0);
    }
  }

  return Ce;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<6, 1>&
MAT::ViscoPlasticNoYieldSurface::CalculateLogElasticStrainInStrainLikeVoigtNotation(
    const CORE::LINALG::Matrix<3, 3>& eigen_vectors,
    const CORE::LINALG::Matrix<3, 1>& eigen_values) const
{
  // trial elastic material logarithmic strains
  static CORE::LINALG::Matrix<3, 3> Ee_trial;
  static CORE::LINALG::Matrix<6, 1> Ee_trial_Vstrain;
  static CORE::LINALG::Matrix<3, 3> tmp3x3;
  static CORE::LINALG::Matrix<3, 3> logEigenValues(true);

  for (unsigned i = 0; i < 3; ++i)
  {
    logEigenValues(i, i) = std::log(eigen_values(i));
  }

  tmp3x3.MultiplyNN(eigen_vectors, logEigenValues);
  Ee_trial.MultiplyNT(tmp3x3, eigen_vectors);

  // transform to strain-like Voigt notation of logarithmic elastic strain
  UTILS::VOIGT::Strains::MatrixToVector(Ee_trial, Ee_trial_Vstrain);

  return Ee_trial_Vstrain;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<3, 3>& MAT::ViscoPlasticNoYieldSurface::CalculateTrialElasticRotation(
    const CORE::LINALG::Matrix<3, 3>& Fe_trial, const CORE::LINALG::Matrix<3, 3>& eigen_vectors,
    const CORE::LINALG::Matrix<3, 1>& eigen_values) const
{
  static CORE::LINALG::Matrix<3, 3> Re_trial;
  static CORE::LINALG::Matrix<3, 3> tmp3x3;
  static CORE::LINALG::Matrix<3, 3> invEigenValues(true);

  for (unsigned i = 0; i < 3; ++i)
  {
    invEigenValues(i, i) = 1.0 / eigen_values(i);
  }

  // trial elastic rotation tensor (F = R * U --> R = F * U^(-1))
  tmp3x3.MultiplyNN(eigen_vectors, invEigenValues);
  // inverse trial stretch
  static CORE::LINALG::Matrix<3, 3> iUe_trial;
  iUe_trial.MultiplyNT(tmp3x3, eigen_vectors);
  Re_trial.MultiplyNN(Fe_trial, iUe_trial);

  return Re_trial;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<3, 3>& MAT::ViscoPlasticNoYieldSurface::CalculateSecondPiolaKirchhoffStresses(
    const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 3>& Re_trial,
    const CORE::LINALG::Matrix<3, 3>& Me) const
{
  static CORE::LINALG::Matrix<3, 3> iF;
  iF.Invert(*defgrd);

  // calculate 2nd PK stress (PK2 = F^(-1) * Re_trial * Me * Re_trial^(T) * F^(-T))
  // from left to right
  static CORE::LINALG::Matrix<3, 3> PK2;
  static CORE::LINALG::Matrix<3, 3> tmp3x3;
  tmp3x3.MultiplyNN(iF, Re_trial);
  PK2.MultiplyNN(tmp3x3, Me);
  tmp3x3.MultiplyNT(PK2, Re_trial);
  PK2.MultiplyNT(tmp3x3, iF);

  return PK2;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoPlasticNoYieldSurface::LocalNewtonLoop(
    CORE::LINALG::Matrix<2, 1>& x, const double dt)
{
  // predictor values
  const double equ_tens_trial_stress = x(0);
  const double flow_resistance_n = x(1);

  // Convergence tolerance for NR
  unsigned iter(0);
  const double tolNR(1.0e-12);
  const unsigned max_iter(100);
  // Jacobian matrix containing linearization
  static CORE::LINALG::Matrix<2, 2> J;
  // Inverse Jacobian
  static CORE::LINALG::Matrix<2, 2> iJ;
  // Increment of the solution variables
  static CORE::LINALG::Matrix<2, 1> dx;
  // Residual of both equations
  static CORE::LINALG::Matrix<2, 1> residual;

  // Newton Loop until convergence of the L2-norm, or maximum iteration reached
  while (true)
  {
    ++iter;

    // execute pre calculations
    const MAT::PreCalculatedTerms terms = PreCalculateTerms(x(0), x(1), dt);

    // Calculate residual and L2 Norm
    residual = CalculateResidual(x(0), equ_tens_trial_stress, x(1), flow_resistance_n, terms);
    double residual_norm2 = residual.Norm2();

    // skip if residual has converged
    if (residual_norm2 < tolNR) break;

    // Check if maximum iteration is reached
    if (iter > max_iter)
      dserror("Local Newton Raphson in ViscoPlasticNoYieldSurface material not converged");

    // calculate linearization
    J = CalculateLinearization(x(0), x(1), terms);

    // invert linearization and calculate increment
    iJ.Invert(J);
    dx.Multiply(-1.0, iJ, residual, 0.0);

    // update solution vector
    x.Update(1.0, dx, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<2, 1>& MAT::ViscoPlasticNoYieldSurface::CalculateResidual(
    const double equ_tens_stress_np, const double equ_tens_trial_stress,
    const double flow_resistance_np, const double flow_resistance_n,
    const PreCalculatedTerms& terms)
{
  // viscosity parameters
  const double m = params_->StrainRateSensitivity();
  const double a = params_->FlowResExp();

  static CORE::LINALG::Matrix<2, 1> residual;
  residual(0) =
      equ_tens_stress_np - equ_tens_trial_stress +
      terms.equ_tens_stress_const * std::pow(terms.equ_tens_stress_flow_res_ratio, (1.0 / m));

  const double temp_abs_pow = std::pow(std::abs(terms.flow_res_sat_deviation), a);
  const double flow_resistance_deriv_np = terms.flow_res_const * temp_abs_pow *
                                          CORE::FADUTILS::Signum(terms.flow_res_sat_deviation) *
                                          std::pow(terms.equ_tens_stress_flow_res_ratio, (1.0 / m));

  residual(1) = flow_resistance_np - flow_resistance_n - flow_resistance_deriv_np;

  return residual;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<2, 2>& MAT::ViscoPlasticNoYieldSurface::CalculateLinearization(
    const double equ_tens_stress_np, const double flow_resistance_np,
    const MAT::PreCalculatedTerms& terms)
{
  static CORE::LINALG::Matrix<2, 2> J;

  // viscosity parameters
  const double m = params_->StrainRateSensitivity();
  const double a = params_->FlowResExp();
  const double b = params_->FlowResSatExp();
  const double flow_res_sat_fac = params_->FlowResSatFac();

  J(0, 0) = 1.0 + terms.equ_tens_stress_const / (m * equ_tens_stress_np) *
                      std::pow(terms.equ_tens_stress_flow_res_ratio, 1.0 / m);

  J(0, 1) = -terms.equ_tens_stress_const / (m * flow_resistance_np) *
            std::pow(terms.equ_tens_stress_flow_res_ratio, 1.0 / m);

  J(1, 0) = -terms.flow_res_const / (m * equ_tens_stress_np) *
                std::pow(std::abs(terms.flow_res_sat_deviation), a) *
                CORE::FADUTILS::Signum(terms.flow_res_sat_deviation) *
                std::pow(terms.equ_tens_stress_flow_res_ratio, 1.0 / m) -
            terms.flow_res_const * b / (m * flow_res_sat_fac) * a *
                std::abs(terms.flow_res_sat_deviation) *
                std::pow(std::abs(terms.flow_res_sat_deviation), a - 2.0) *
                std::pow(terms.equ_tens_stress_flow_res_ratio, (1.0 - b - m) / m);

  J(1, 1) = 1.0 +
            terms.flow_res_const / (m * flow_resistance_np) *
                std::pow(std::abs(terms.flow_res_sat_deviation), a) *
                CORE::FADUTILS::Signum(terms.flow_res_sat_deviation) *
                std::pow(terms.equ_tens_stress_flow_res_ratio, 1.0 / m) +
            terms.flow_res_const / flow_res_sat_fac *
                std::pow(terms.equ_tens_stress_flow_res_ratio, (1.0 - b) / m) * a *
                std::abs(terms.flow_res_sat_deviation) *
                std::pow(std::abs(terms.flow_res_sat_deviation), a - 2.0) * (1.0 + b / m);

  return J;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PreCalculatedTerms MAT::ViscoPlasticNoYieldSurface::PreCalculateTerms(
    const double equ_tens_stress_np, const double flow_resistance_np, const double dt)
{
  MAT::PreCalculatedTerms terms;

  // Elasticity parameters
  const double E = params_->Young();
  const double nu = params_->Nue();
  const double G = E / (2.0 * (1.0 + nu));

  // viscosity parameters
  const double T = params_->Temperature();
  const double A = params_->PreExpFac();
  const double Q = params_->ActivationEnergy();
  const double R = params_->GasConstant();
  const double m = params_->StrainRateSensitivity();
  const double H_0 = params_->FlowResPreFac();
  const double flow_res_sat_fac = params_->FlowResSatFac();
  const double b = params_->FlowResSatExp();

  const double expterm = std::exp(-Q / (R * T));
  terms.equ_tens_stress_const = 3.0 * G * dt * A * expterm;
  terms.flow_res_const = dt * H_0 * A * expterm;
  terms.equ_tens_stress_flow_res_ratio = equ_tens_stress_np / flow_resistance_np;
  terms.flow_res_sat_deviation =
      1.0 - flow_resistance_np /
                (flow_res_sat_fac * std::pow(terms.equ_tens_stress_flow_res_ratio, b / m));

  return terms;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoPlasticNoYieldSurface::SetupCmat(
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& cmat)
{
  // get material parameters
  // Young's modulus
  const double youngs_mod = params_->Young();
  // Poisson's ratio
  const double nue = params_->Nue();

  // isotropic elasticity tensor C in Voigt matrix notation
  //                       [ 1-nu     nu     nu |          0    0    0 ]
  //                       [        1-nu     nu |          0    0    0 ]
  //           E           [               1-nu |          0    0    0 ]
  //   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~  ~~~ ]
  //       (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0    0 ]
  //                       [                    |      (1-2*nu)/2    0 ]
  //                       [ symmetric          |           (1-2*nu)/2 ]
  //
  const double mfac = youngs_mod / ((1.0 + nue) * (1.0 - 2.0 * nue));  // factor

  cmat.Clear();
  // write non-zero components --- axial
  cmat(0, 0) = mfac * (1.0 - nue);
  cmat(0, 1) = mfac * nue;
  cmat(0, 2) = mfac * nue;
  cmat(1, 0) = mfac * nue;
  cmat(1, 1) = mfac * (1.0 - nue);
  cmat(1, 2) = mfac * nue;
  cmat(2, 0) = mfac * nue;
  cmat(2, 1) = mfac * nue;
  cmat(2, 2) = mfac * (1.0 - nue);
  // write non-zero components --- shear
  cmat(3, 3) = mfac * 0.5 * (1.0 - 2.0 * nue);
  cmat(4, 4) = mfac * 0.5 * (1.0 - 2.0 * nue);
  cmat(5, 5) = mfac * 0.5 * (1.0 - 2.0 * nue);
}