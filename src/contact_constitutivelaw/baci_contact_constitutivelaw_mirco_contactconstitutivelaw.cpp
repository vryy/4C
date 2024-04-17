/*----------------------------------------------------------------------------*/
/*! \file
\brief implements a mirco contact constitutive law
\level 3

*----------------------------------------------------------------------*/


#include "baci_contact_constitutivelaw_mirco_contactconstitutivelaw.hpp"

#include "baci_global_data.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_mat_par_bundle.hpp"

#ifdef BACI_WITH_MIRCO

#include <mirco_evaluate.h>
#include <mirco_topology.h>
#include <mirco_topologyutilities.h>

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::MircoConstitutiveLawParams::MircoConstitutiveLawParams(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container)
    : CONTACT::CONSTITUTIVELAW::Parameter(container),
      firstmatid_(*container->Get<int>("FirstMatID")),
      secondmatid_(*container->Get<int>("SecondMatID")),
      lateralLength_(*container->Get<double>("LateralLength")),
      resolution_(*container->Get<int>("Resolution")),
      pressureGreenFunFlag_(*container->Get<bool>("PressureGreenFunFlag")),
      initialTopologyStdDeviation_(*container->Get<double>("InitialTopologyStdDeviation")),
      hurstExponent_(*container->Get<double>("HurstExponent")),
      randomTopologyFlag_(*container->Get<bool>("RandomTopologyFlag")),
      randomSeedFlag_(*container->Get<bool>("RandomSeedFlag")),
      randomGeneratorSeed_(*container->Get<int>("RandomGeneratorSeed")),
      tolerance_(*container->Get<double>("Tolerance")),
      maxIteration_(*container->Get<int>("MaxIteration")),
      warmStartingFlag_(*container->Get<bool>("WarmStartingFlag")),
      finiteDifferenceFraction_(*container->Get<double>("FiniteDifferenceFraction")),
      activeGapTolerance_(*container->Get<double>("ActiveGapTolerance")),
      topologyFilePath_(*(container->Get<std::string>("TopologyFilePath")))
{
  this->SetParameters();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw>
CONTACT::CONSTITUTIVELAW::MircoConstitutiveLawParams::CreateConstitutiveLaw()
{
  return Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::MircoConstitutiveLaw(this));
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::MircoConstitutiveLaw::MircoConstitutiveLaw(
    CONTACT::CONSTITUTIVELAW::MircoConstitutiveLawParams* params)
    : params_(params)
{
}
void CONTACT::CONSTITUTIVELAW::MircoConstitutiveLawParams::SetParameters()
{
  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> firstmat =
      GLOBAL::Problem::Instance(probinst)->Materials()->ById(GetFirstMatID());
  Teuchos::RCP<MAT::PAR::Material> secondmat =
      GLOBAL::Problem::Instance(probinst)->Materials()->ById(GetSecondMatID());

  const double E1 = *firstmat->Get<double>("YOUNG");
  const double E2 = *secondmat->Get<double>("YOUNG");
  const double nu1 = *firstmat->Get<double>("NUE");
  const double nu2 = *secondmat->Get<double>("NUE");

  // Composite Young's modulus
  compositeYoungs_ = pow(((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2)) / E2), -1);

  // Composite Shear modulus
  double G1 = E1 / (2 * (1 + nu1));
  double G2 = E2 / (2 * (1 + nu2));
  double CompositeShear = pow(((2 - nu1) / (4 * G1) + (2 - nu2) / (4 * G2)), -1);

  // Composite Poisson's ratio
  compositePoissonsRatio_ = compositeYoungs_ / (2 * CompositeShear) - 1;

  gridSize_ = lateralLength_ / (pow(2, resolution_) + 1);

  // Shape factors (See section 3.3 of https://doi.org/10.1007/s00466-019-01791-3)
  // These are the shape factors to calculate the elastic compliance correction of the micro-scale
  // contact constitutive law for various resolutions.
  // NOTE: Currently MIRCO works for resouluion of 1 to 8. The following map store the shape
  // factors for resolution of 1 to 8.

  // The following pressure based constants are calculated by solving a flat indentor problem in
  // MIRCO using the pressure based Green function described in Pohrt and Li (2014).
  // http://dx.doi.org/10.1134/s1029959914040109
  const std::map<int, double> shape_factors_pressure{{1, 0.961389237917602}, {2, 0.924715342432435},
      {3, 0.899837531880697}, {4, 0.884976751041942}, {5, 0.876753783192863},
      {6, 0.872397956576882}, {7, 0.871958228537090}, {8, 0.882669916668780}};

  // The following force based constants are taken from Table 1 of Bonari et al. (2020).
  // https://doi.org/10.1007/s00466-019-01791-3
  const std::map<int, double> shape_factors_force{{1, 0.778958541513360}, {2, 0.805513388666376},
      {3, 0.826126871395416}, {4, 0.841369158110513}, {5, 0.851733020725652},
      {6, 0.858342234203154}, {7, 0.862368243479785}, {8, 0.864741597831785}};

  const double ShapeFactor = pressureGreenFunFlag_ ? shape_factors_pressure.at(resolution_)
                                                   : shape_factors_force.at(resolution_);

  elasticComplianceCorrection_ = lateralLength_ * compositeYoungs_ / ShapeFactor;

  const int iter = int(ceil((lateralLength_ - (gridSize_ / 2)) / gridSize_));
  meshgrid_ = Teuchos::Ptr(new std::vector<double>(iter));
  MIRCO::CreateMeshgrid(*meshgrid_, iter, gridSize_);

  // Setup Topology
  const int N = pow(2, resolution_);
  topology_ = Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix(N + 1, N + 1));

  Teuchos::RCP<MIRCO::TopologyGeneration> surfacegenerator;
  // creating the correct surface object
  MIRCO::CreateSurfaceObject(resolution_, initialTopologyStdDeviation_, hurstExponent_,
      randomSeedFlag_, topologyFilePath_, randomTopologyFlag_, randomGeneratorSeed_,
      surfacegenerator);
  surfacegenerator->GetSurface(*topology_);

  auto max_and_mean = MIRCO::ComputeMaxAndMean(*topology_);
  maxTopologyHeight_ = max_and_mean.max_;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::MircoConstitutiveLaw::Evaluate(double gap)
{
  if (gap + params_->GetOffset() > 0.0)
  {
    dserror("You should not be here. The Evaluate function is only tested for active nodes. ");
  }
  if (-(gap + params_->GetOffset()) < params_->GetActiveGapTolerance())
  {
    return 0.0;
  }

  double pressure = 0.0;
  MIRCO::Evaluate(pressure, -(gap + params_->GetOffset()), params_->GetLateralLength(),
      params_->GetGridSize(), params_->GetTolerance(), params_->GetMaxIteration(),
      params_->GetCompositeYoungs(), params_->GetCompositePoissonsRatio(),
      params_->GetWarmStartingFlag(), params_->GetComplianceCorrection(), *params_->GetTopology(),
      params_->GetMaxTopologyHeight(), *params_->GetMeshGrid(), params_->GetPressureGreenFunFlag());

  return (-1 * pressure);
}  // end of mirco_coconstlaw evaluate
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double CONTACT::CONSTITUTIVELAW::MircoConstitutiveLaw::EvaluateDeriv(double gap)
{
  if (gap + params_->GetOffset() > 0.0)
  {
    dserror("You should not be here. The Evaluate function is only tested for active nodes.");
  }
  if (-(gap + params_->GetOffset()) < params_->GetActiveGapTolerance())
  {
    return 0.0;
  }

  double pressure1 = 0.0;
  double pressure2 = 0.0;
  // using backward difference approach
  MIRCO::Evaluate(pressure1, -1.0 * (gap + params_->GetOffset()), params_->GetLateralLength(),
      params_->GetGridSize(), params_->GetTolerance(), params_->GetMaxIteration(),
      params_->GetCompositeYoungs(), params_->GetCompositePoissonsRatio(),
      params_->GetWarmStartingFlag(), params_->GetComplianceCorrection(), *params_->GetTopology(),
      params_->GetMaxTopologyHeight(), *params_->GetMeshGrid(), params_->GetPressureGreenFunFlag());
  MIRCO::Evaluate(pressure2,
      -(1 - params_->GetFiniteDifferenceFraction()) * (gap + params_->GetOffset()),
      params_->GetLateralLength(), params_->GetGridSize(), params_->GetTolerance(),
      params_->GetMaxIteration(), params_->GetCompositeYoungs(),
      params_->GetCompositePoissonsRatio(), params_->GetWarmStartingFlag(),
      params_->GetComplianceCorrection(), *params_->GetTopology(), params_->GetMaxTopologyHeight(),
      *params_->GetMeshGrid(), params_->GetPressureGreenFunFlag());
  return ((pressure1 - pressure2) /
          (-(params_->GetFiniteDifferenceFraction()) * (gap + params_->GetOffset())));
}

FOUR_C_NAMESPACE_CLOSE

#endif
