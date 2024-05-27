/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a hyperelastic constituent with a damage process

\level 3
*/
/*----------------------------------------------------------------------*/

#include "4C_mixture_constituent_elasthyper_damage.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper_service.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_matelast_isoneohooke.hpp"
#include "4C_mixture_elastin_membrane_prestress_strategy.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

// Constructor for the parameter class
MIXTURE::PAR::MixtureConstituentElastHyperDamage::MixtureConstituentElastHyperDamage(
    const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata)
    : MixtureConstituentElastHyperBase(matdata),
      damage_function_id_(matdata->Get<int>("DAMAGE_FUNCT"))
{
  // nothing to do here
}

// Create an instance of MIXTURE::MixtureConstituentElastHyper from the parameters
std::unique_ptr<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituentElastHyperDamage::CreateConstituent(int id)
{
  return std::unique_ptr<MIXTURE::MixtureConstituentElastHyperDamage>(
      new MIXTURE::MixtureConstituentElastHyperDamage(this, id));
}

// Constructor of the constituent holding the material parameters
MIXTURE::MixtureConstituentElastHyperDamage::MixtureConstituentElastHyperDamage(
    MIXTURE::PAR::MixtureConstituentElastHyperDamage* params, int id)
    : MixtureConstituentElastHyperBase(params, id), params_(params)
{
  // nothing to do here
}

// Returns the material type
CORE::Materials::MaterialType MIXTURE::MixtureConstituentElastHyperDamage::MaterialType() const
{
  return CORE::Materials::mix_elasthyper_damage;
}

// Pack the constituent
void MIXTURE::MixtureConstituentElastHyperDamage::PackConstituent(
    CORE::COMM::PackBuffer& data) const
{
  MixtureConstituentElastHyperBase::PackConstituent(data);

  CORE::COMM::ParObject::AddtoPack(data, current_reference_growth_);
}

// Unpack the constituent
void MIXTURE::MixtureConstituentElastHyperDamage::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureConstituentElastHyperBase::UnpackConstituent(position, data);

  CORE::COMM::ParObject::ExtractfromPack(position, data, current_reference_growth_);
}

// Reads the element from the input file
void MIXTURE::MixtureConstituentElastHyperDamage::ReadElement(
    int numgp, INPUT::LineDefinition* linedef)
{
  MixtureConstituentElastHyperBase::ReadElement(numgp, linedef);

  current_reference_growth_.resize(numgp, 1.0);
}

// Updates all summands
void MIXTURE::MixtureConstituentElastHyperDamage::Update(CORE::LINALG::Matrix<3, 3> const& defgrd,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  const auto& reference_coordinates = params.get<CORE::LINALG::Matrix<3, 1>>("gp_coords_ref");

  double totaltime = params.get<double>("total time", -1);
  if (totaltime < 0.0)
  {
    FOUR_C_THROW("Parameter 'total time' could not be read!");
  }

  current_reference_growth_[gp] =
      GLOBAL::Problem::Instance()
          ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(params_->damage_function_id_ - 1)
          .Evaluate(reference_coordinates.A(), totaltime, 0);

  MixtureConstituentElastHyperBase::Update(defgrd, params, gp, eleGID);
}

double MIXTURE::MixtureConstituentElastHyperDamage::GetGrowthScalar(int gp) const
{
  return current_reference_growth_[gp];
}

void MIXTURE::MixtureConstituentElastHyperDamage::Evaluate(const CORE::LINALG::Matrix<3, 3>& F,
    const CORE::LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  FOUR_C_THROW("This constituent does not support Evaluation without an elastic part.");
}

void MIXTURE::MixtureConstituentElastHyperDamage::EvaluateElasticPart(
    const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<3, 3>& iFextin,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  // Compute total inelastic deformation gradient
  static CORE::LINALG::Matrix<3, 3> iFin(false);
  iFin.MultiplyNN(iFextin, prestretch_tensor(gp));

  // Evaluate 3D elastic part
  MAT::ElastHyperEvaluateElasticPart(
      F, iFin, S_stress, cmat, Summands(), SummandProperties(), gp, eleGID);
}
FOUR_C_NAMESPACE_CLOSE
