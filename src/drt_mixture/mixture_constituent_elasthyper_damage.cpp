/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a hyperelastic constituent with a damage process

\level 3
*/
/*----------------------------------------------------------------------*/

#include "mixture_constituent_elasthyper_damage.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/multiplicative_split_defgrad_elasthyper_service.H"
#include "../drt_matelast/elast_aniso_structuraltensor_strategy.H"
#include "../drt_matelast/elast_isoneohooke.H"
#include "elastin_membrane_prestress_strategy.H"

// Constructor for the parameter class
MIXTURE::PAR::MixtureConstituent_ElastHyperDamage::MixtureConstituent_ElastHyperDamage(
    const Teuchos::RCP<MAT::PAR::Material>& matdata, const double ref_mass_fraction)
    : MixtureConstituent_ElastHyperBase(matdata, ref_mass_fraction),
      damage_function_id_(matdata->GetInt("DAMAGE_FUNCT"))
{
  // nothing to do here
}

// Create an instance of MIXTURE::MixtureConstituent_ElastHyper from the parameters
Teuchos::RCP<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituent_ElastHyperDamage::CreateConstituent()
{
  return Teuchos::rcp(new MIXTURE::MixtureConstituent_ElastHyperDamage(this));
}

// Constructor of the constituent holding the material parameters
MIXTURE::MixtureConstituent_ElastHyperDamage::MixtureConstituent_ElastHyperDamage(
    MIXTURE::PAR::MixtureConstituent_ElastHyperDamage* params)
    : MixtureConstituent_ElastHyperBase(params), params_(params)
{
  // nothing to do here
}

// Returns the material type
INPAR::MAT::MaterialType MIXTURE::MixtureConstituent_ElastHyperDamage::MaterialType() const
{
  return INPAR::MAT::mix_elasthyper_damage;
}

// Pack the constituent
void MIXTURE::MixtureConstituent_ElastHyperDamage::PackConstituent(DRT::PackBuffer& data) const
{
  MixtureConstituent_ElastHyperBase::PackConstituent(data);

  DRT::ParObject::AddtoPack(data, current_reference_density_);
}

// Unpack the constituent
void MIXTURE::MixtureConstituent_ElastHyperDamage::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureConstituent_ElastHyperBase::UnpackConstituent(position, data);

  DRT::ParObject::ExtractfromPack(position, data, current_reference_density_);
}

// Reads the element from the input file
void MIXTURE::MixtureConstituent_ElastHyperDamage::ReadElement(
    int numgp, DRT::INPUT::LineDefinition* linedef)
{
  MixtureConstituent_ElastHyperBase::ReadElement(numgp, linedef);

  current_reference_density_.resize(numgp, params_->RefMassFraction() * InitialRefDensity());
}

// Updates all summands
void MIXTURE::MixtureConstituent_ElastHyperDamage::Update(LINALG::Matrix<3, 3> const& defgrd,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  LINALG::Matrix<1, 3> gprefecoord(true);  // gp coordinates in reference configuration
  gprefecoord = params.get<LINALG::Matrix<1, 3>>("gprefecoord");

  double totaltime = params.get<double>("total time", -1);
  if (totaltime < 0.0)
  {
    dserror("Parameter 'total time' could not be read!");
  }

  current_reference_density_[gp] = params_->RefMassFraction() * InitialRefDensity() *
                                   DRT::Problem::Instance()
                                       ->Funct(params_->damage_function_id_ - 1)
                                       .Evaluate(0, gprefecoord.A(), totaltime);

  MixtureConstituent_ElastHyperBase::Update(defgrd, params, gp, eleGID);
}

// Returns the reference mass fraction of the constituent
double MIXTURE::MixtureConstituent_ElastHyperDamage::CurrentRefDensity(int gp) const
{
  return current_reference_density_[gp];
}

void MIXTURE::MixtureConstituent_ElastHyperDamage::Evaluate(const LINALG::Matrix<3, 3>& F,
    const LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  dserror("This constituent does not support Evaluation without an elastic part.");
}

void MIXTURE::MixtureConstituent_ElastHyperDamage::EvaluateElasticPart(
    const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<3, 3>& iFextin,
    Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
    int gp, int eleGID)
{
  // Compute total inelastic deformation gradient
  static LINALG::Matrix<3, 3> iFin(false);
  iFin.MultiplyNN(iFextin, PrestretchTensor(gp));

  // Evaluate 3D elastic part
  MAT::ElastHyperEvaluateElasticPart(
      F, iFin, S_stress, cmat, Summands(), SummandProperties(), gp, eleGID);

  S_stress.Scale(CurrentRefDensity(gp));
  cmat.Scale(CurrentRefDensity(gp));
}