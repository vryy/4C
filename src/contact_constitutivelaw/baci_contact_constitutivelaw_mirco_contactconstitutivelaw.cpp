/*----------------------------------------------------------------------------*/
/*! \file
\brief implements a mirco contact constitutive law
\level 3

*----------------------------------------------------------------------*/

#ifdef BACI_WITH_MIRCO

#include "baci_contact_constitutivelaw_mirco_contactconstitutivelaw.H"

#include "baci_lib_globalproblem.H"
#include "baci_mat_par_bundle.H"

#include <mirco_evaluate.h>
#include <mirco_topology.h>
#include <mirco_topologyutilities.h>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>

#include <vector>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::MircoConstitutiveLawParams::MircoConstitutiveLawParams(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container)
    : CONTACT::CONSTITUTIVELAW::Parameter(container),
      firstmatid_(container->GetDouble("FirstMatID")),
      secondmatid_(container->GetDouble("SecondMatID")),
      lateralLength_(container->GetDouble("LateralLength")),
      resolution_(container->GetDouble("Resolution")),
      initialTopologyStdDeviation_(container->GetDouble("InitialTopologyStdDeviation")),
      hurstExponent_(container->GetDouble("HurstExponent")),
      randomTopologyFlag_(container->GetDouble("RandomTopologyFlag")),
      randomSeedFlag_(container->GetDouble("RandomSeedFlag")),
      randomGeneratorSeed_(container->GetDouble("RandomGeneratorSeed")),
      tolerance_(container->GetDouble("Tolerance")),
      maxIteration_(container->GetDouble("MaxIteration")),
      warmStartingFlag_(container->GetDouble("WarmStartingFlag")),
      finiteDifferenceFraction_(container->GetDouble("FiniteDifferenceFraction")),
      activeGapTolerance_(container->GetDouble("ActiveGapTolerance")),
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
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> firstmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(GetFirstMatID());
  Teuchos::RCP<MAT::PAR::Material> secondmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(GetSecondMatID());

  const double E1 = firstmat->GetDouble("YOUNG");
  const double E2 = secondmat->GetDouble("YOUNG");
  const double nu1 = firstmat->GetDouble("NUE");
  const double nu2 = secondmat->GetDouble("NUE");

  compositeYoungs_ = pow(((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2)) / E2), -1);

  gridSize_ = lateralLength_ / (pow(2, resolution_) + 1);

  // Correction factor vector
  // These are the correction factors to calculate the elastic compliance of the micro-scale contact
  // constitutive law for various resolutions. These constants are taken from Table 1 of Bonari et
  // al. (2020). https://doi.org/10.1007/s00466-019-01791-3
  std::vector<double> alpha_con{0.778958541513360, 0.805513388666376, 0.826126871395416,
      0.841369158110513, 0.851733020725652, 0.858342234203154, 0.862368243479785,
      0.864741597831785};
  const double alpha = alpha_con[resolution_ - 1];
  elasticComplianceCorrection_ = lateralLength_ * compositeYoungs_ / alpha;

  const int iter = int(ceil((lateralLength_ - (gridSize_ / 2)) / gridSize_));
  meshgrid_ = Teuchos::Ptr(new std::vector<double>(iter));
  MIRCO::CreateMeshgrid(*meshgrid_, iter, gridSize_);

  // Setup Topology
  const int N = pow(2, resolution_);
  topology_ = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double>(N + 1, N + 1));

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
      params_->GetCompositeYoungs(), params_->GetWarmStartingFlag(),
      params_->GetComplianceCorrection(), *params_->GetTopology(), params_->GetMaxTopologyHeight(),
      *params_->GetMeshGrid());

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
      params_->GetCompositeYoungs(), params_->GetWarmStartingFlag(),
      params_->GetComplianceCorrection(), *params_->GetTopology(), params_->GetMaxTopologyHeight(),
      *params_->GetMeshGrid());
  MIRCO::Evaluate(pressure2,
      -(1 - params_->GetFiniteDifferenceFraction()) * (gap + params_->GetOffset()),
      params_->GetLateralLength(), params_->GetGridSize(), params_->GetTolerance(),
      params_->GetMaxIteration(), params_->GetCompositeYoungs(), params_->GetWarmStartingFlag(),
      params_->GetComplianceCorrection(), *params_->GetTopology(), params_->GetMaxTopologyHeight(),
      *params_->GetMeshGrid());
  return ((pressure1 - pressure2) /
          (-(params_->GetFiniteDifferenceFraction()) * (gap + params_->GetOffset())));
}
#endif
