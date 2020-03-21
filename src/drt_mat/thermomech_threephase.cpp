/*! \file
\brief A thermo-mechnical pseudo-solid material used for additive manufacturing applications

\level 3

\maintainer Sebastian Proell
*/

#include "thermomech_threephase.H"
#include "consolidation.H"
#include "fouriervar.H"

#include "../drt_lib/drt_globalproblem.H"
#include "matpar_bundle.H"
#include "stvenantkirchhoff.H"
#include "thermostvenantkirchhoff.H"

MAT::PAR::ThermoMechThreePhase::ThermoMechThreePhase(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      youngsfunct_(*(matdata->Get<std::vector<int>>("YOUNGFUNCT"))),
      poissonratio_(matdata->GetDouble("NUE")),
      density_(matdata->GetDouble("DENS")),
      thermexpansfunct_(*(matdata->Get<std::vector<int>>("THEXPANSFUNCT"))),
      thetaref_(matdata->GetDouble("THETAREF")),
      fouriervarmat_(matdata->GetInt("THERMOMAT")),
      consolmat_(matdata->GetInt("CONSOLMAT"))
{
  if (poissonratio_ >= 0.5 || poissonratio_ < -1.) dserror("Poisson's ratio must be in [-1;0.5)");
  if (youngsfunct_.size() != 3)
    dserror("Exactly three functions for Young's modulus in powder, melt and solid are required.");
  if (thermexpansfunct_.size() != 3)
    dserror(
        "Exactly three functions for thermal expansion in powder, melt and solid are required.");
}

Teuchos::RCP<MAT::Material> MAT::PAR::ThermoMechThreePhase::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ThermoMechThreePhase(this));
}

DRT::ParObject* MAT::ThermoMechThreePhaseType::Create(const std::vector<char>& data)
{
  auto* material = new MAT::ThermoMechThreePhase();
  material->Unpack(data);
  return material;
}

MAT::ThermoMechThreePhase::ThermoMechThreePhase()
    : params_(nullptr), consolidation_(Teuchos::null), fouriervar_(Teuchos::null)
{
}

MAT::ThermoMechThreePhase::ThermoMechThreePhase(MAT::PAR::ThermoMechThreePhase* params)
    : params_(params), consolidation_(Teuchos::null), fouriervar_(Teuchos::null)
{
  const int consolMatId = params_->consolmat_;
  auto mat = MAT::Material::Factory(consolMatId);
  if (mat == Teuchos::null) dserror("Failed to create consolidation material, id=%d", consolMatId);
  consolidation_ = Teuchos::rcp_dynamic_cast<MAT::Consolidation>(mat);

  CreateThermoMaterial();
}

void MAT::ThermoMechThreePhase::CreateThermoMaterial()
{
  const int fourierMatId = this->params_->fouriervarmat_;
  auto mat2 = MAT::Material::Factory(fourierMatId);
  if (mat2 == Teuchos::null) dserror("Failed to create FourierVar material, id=%d", fourierMatId);
  this->fouriervar_ = Teuchos::rcp_dynamic_cast<MAT::FourierVar>(mat2);

  // pass the same object for Consolidation to thermo material implementation
  fouriervar_->InjectConsolidation(consolidation_);
}

void MAT::ThermoMechThreePhase::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // packing matid will enable retrieving all parameters including the id of FourierVar
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // pack consolidation manager
  if (consolidation_ != Teuchos::null) consolidation_->Pack(data);
}

void MAT::ThermoMechThreePhase::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ThermoMechThreePhase*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  // only do this if not in post-processing
  if (params_ != NULL)
  {
    // unpack consolidation manager
    // get the data for the consolidation manager
    std::vector<char> consoldata;
    ExtractfromPack(position, data, consoldata);
    // construct it and unpack data
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(params_->consolmat_);
    consolidation_ = Teuchos::rcp_dynamic_cast<MAT::Consolidation>(mat);
    consolidation_->Unpack(consoldata);

    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);

    CreateThermoMaterial();
  }
}

void MAT::ThermoMechThreePhase::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  consolidation_->Setup(numgp);
}

/*----------------------------------------------------------------------*|
 *----------------------------------------------------------------------*/
void MAT::ThermoMechThreePhase::Reinit(double temperature, unsigned gp)
{
  currentTemperature_ = temperature;
  currentGp_ = gp;
  fouriervar_->Reinit(currentTemperature_, currentGp_);
}

/*----------------------------------------------------------------------*|
 *----------------------------------------------------------------------*/
void MAT::ThermoMechThreePhase::Reinit(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, double temperature, unsigned gp)
{
  Reinit(temperature, gp);
  currentGlstrain_ = glstrain;
}


/*----------------------------------------------------------------------*|
 * evaluate stress and constitutive tensor                               |
 *----------------------------------------------------------------------*/
void MAT::ThermoMechThreePhase::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int eleGID)
{
  // fixme this backwards compatibility modification should be moved outside
  double temperature = params.get<double>("scalartemp", currentTemperature_);
  unsigned gp = params.get<int>("gp", currentGp_);
  Reinit(defgrd, glstrain, temperature, gp);  // fixme call this before

  SetupCmat(*cmat);
  // purely mechanical part
  stress->MultiplyNN(*cmat, *currentGlstrain_);

  // additive thermal part
  double Tref = params_->thetaref_;
  double m = STModulus();

  // loop over the element nodes, non-zero entries only in main directions
  for (int i = 0; i < 3; ++i) (*stress)(i, 0) += m * (currentTemperature_ - Tref);
}

/*----------------------------------------------------------------------*|
 *----------------------------------------------------------------------*/
void MAT::ThermoMechThreePhase::GetdSdT(LINALG::Matrix<6, 1>* dS_dT)
{
  // calculate derivative of cmat w.r.t. T_{n+1}
  LINALG::Matrix<6, 6> cmat_T(false);
  SetupCmat_dT(cmat_T);

  // mechanical part
  dS_dT->MultiplyNN(cmat_T, *currentGlstrain_);

  // additive thermal part
  double Tref = params_->thetaref_;
  double m_dT = STModulus_dT();
  double m = STModulus();

  // compute derivative (product rule!)
  for (int i = 0; i < 3; ++i) (*dS_dT)(i, 0) += m_dT * (currentTemperature_ - Tref) + m;
}

void MAT::ThermoMechThreePhase::SetupCmat(LINALG::Matrix<6, 6>& cmat)
{
  const double nu = params_->poissonratio_;
  double Emod = GetMaterialParameter(params_->youngsfunct_);

  StVenantKirchhoff::FillCmat(cmat, Emod, nu);
}

void MAT::ThermoMechThreePhase::SetupCmat_dT(LINALG::Matrix<6, 6>& derivcmat)
{
  const double E_T = GetMaterialParameterThermalDerivative(params_->youngsfunct_);
  // Poisson's ratio (Querdehnzahl)
  const double nu = params_->poissonratio_;

  StVenantKirchhoff::FillCmat(derivcmat, E_T, nu);
}



double MAT::ThermoMechThreePhase::STModulus()
{
  const double Emod = GetMaterialParameter(params_->youngsfunct_);

  const double nu = params_->poissonratio_;

  const double thermexpans = GetMaterialParameter(params_->thermexpansfunct_);

  // plane strain, rotational symmetry
  // E / (1+nu)
  const double c1 = Emod / (1.0 + nu);
  // (E*nu) / ((1+nu)(1-2nu))
  const double b1 = c1 * nu / (1.0 - 2.0 * nu);

  // build the lame constants
  //            E
  //   mu = --------
  //        2*(1+nu)
  //                  E*nu
  //   lambda = ----------------
  //            (1+nu)*(1-2*nu)
  //
  //  \f \mu =  \frac{E}{2(1+\nu)} \f
  const double mu = 0.5 * c1;
  // lambda
  // \f \frac{E\,\nu}{(1-2\nu)(1+\nu)} \f
  const double lambda = b1;

  // stress-temperature modulus
  // \f m\, = \, -(2\,\cdot \mu \, +\, 3\cdot\lambda)\cdot\varalpha_T \f
  const double stmodulus = (-1.0) * (2.0 * mu + 3.0 * lambda) * thermexpans;

  return stmodulus;
}

double MAT::ThermoMechThreePhase::STModulus_dT()
{
  const double Emod = GetMaterialParameter(params_->youngsfunct_);
  const double Ederiv = GetMaterialParameterThermalDerivative(params_->youngsfunct_);

  const double thermexpans = GetMaterialParameter(params_->thermexpansfunct_);
  const double thermexpans_deriv =
      GetMaterialParameterThermalDerivative(params_->thermexpansfunct_);

  // initialize the parameters for the lame constants
  const double nu = params_->poissonratio_;

  // plane strain, rotational symmetry
  // E / (1+nu)
  const double c1_deriv = Ederiv / (1.0 + nu);
  const double c1 = Emod / (1.0 + nu);
  // (E . nu) / ((1+nu)(1-2nu))
  const double b1_deriv = c1_deriv * nu / (1.0 - 2.0 * nu);
  const double b1 = c1 * nu / (1.0 - 2.0 * nu);

  // build the lame constants
  const double mu = 0.5 * c1;
  const double mu_deriv = 0.5 * c1_deriv;
  const double lambda = b1;
  const double lambda_deriv = b1_deriv;

  // build the derivative of the stress-temperature modulus w.r.t. T_{n+1}
  // m =   -(2 . mu_deriv + 3 . lambda_deriv) . varalpha
  //       +(-(2 . mu + 3 . lambda)) . varalpha_deriv
  const double stmodulus_T = (-1.0) * (2.0 * mu_deriv + 3.0 * lambda_deriv) * thermexpans +
                             (-1.0) * (2.0 * mu + 3.0 * lambda) * thermexpans_deriv;
  return stmodulus_T;
}
double MAT::ThermoMechThreePhase::GetMaterialParameter(const std::vector<int>& functions) const
{
  return consolidation_->EvaluateTempDependentFunction(currentTemperature_, currentGp_, functions);
}

double MAT::ThermoMechThreePhase::GetMaterialParameterThermalDerivative(
    const std::vector<int>& functions) const
{
  return consolidation_->EvaluateTempDependentDerivative(
      currentTemperature_, currentGp_, functions);
}

void MAT::ThermoMechThreePhase::ResetCurrentState() { fouriervar_->ResetCurrentState(); }

void MAT::ThermoMechThreePhase::CommitCurrentState() { fouriervar_->CommitCurrentState(); }

// Forwarding pure thermo calls to FourierVar implementation

void MAT::ThermoMechThreePhase::Evaluate(const LINALG::Matrix<3, 1>& gradtemp,
    LINALG::Matrix<3, 3>& cmat, LINALG::Matrix<3, 1>& heatflux) const
{
  fouriervar_->Evaluate(gradtemp, cmat, heatflux);
}

void MAT::ThermoMechThreePhase::Evaluate(const LINALG::Matrix<2, 1>& gradtemp,
    LINALG::Matrix<2, 2>& cmat, LINALG::Matrix<2, 1>& heatflux) const
{
  fouriervar_->Evaluate(gradtemp, cmat, heatflux);
}

void MAT::ThermoMechThreePhase::Evaluate(const LINALG::Matrix<1, 1>& gradtemp,
    LINALG::Matrix<1, 1>& cmat, LINALG::Matrix<1, 1>& heatflux) const
{
  fouriervar_->Evaluate(gradtemp, cmat, heatflux);
}

void MAT::ThermoMechThreePhase::ConductivityDerivT(LINALG::Matrix<3, 3>& dCondDT) const
{
  fouriervar_->ConductivityDerivT(dCondDT);
}

void MAT::ThermoMechThreePhase::ConductivityDerivT(LINALG::Matrix<2, 2>& dCondDT) const
{
  fouriervar_->ConductivityDerivT(dCondDT);
}

void MAT::ThermoMechThreePhase::ConductivityDerivT(LINALG::Matrix<1, 1>& dCondDT) const
{
  fouriervar_->ConductivityDerivT(dCondDT);
}

double MAT::ThermoMechThreePhase::CapacityDerivT() const { return fouriervar_->CapacityDerivT(); }

double MAT::ThermoMechThreePhase::Capacity() const { return fouriervar_->Capacity(); }

double MAT::ThermoMechThreePhase::Density() const { return params_->density_; }


void MAT::ThermoMechThreePhase::StressTemperatureModulusAndDeriv(
    LINALG::Matrix<6, 1>& stm, LINALG::Matrix<6, 1>& stm_dT)
{
  const double stm_scalar = STModulus();
  const double stm_dT_scalar = STModulus_dT();

  ThermoStVenantKirchhoff::FillCthermo(stm, stm_scalar);
  ThermoStVenantKirchhoff::FillCthermo(stm_dT, stm_dT_scalar);
}

// static initializer
MAT::ThermoMechThreePhaseType MAT::ThermoMechThreePhaseType::instance_;
