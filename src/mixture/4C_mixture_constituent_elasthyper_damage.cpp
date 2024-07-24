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
    const Core::Mat::PAR::Parameter::Data& matdata)
    : MixtureConstituentElastHyperBase(matdata),
      damage_function_id_(matdata.parameters.get<int>("DAMAGE_FUNCT"))
{
  // nothing to do here
}

// Create an instance of MIXTURE::MixtureConstituentElastHyper from the parameters
std::unique_ptr<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituentElastHyperDamage::create_constituent(int id)
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
Core::Materials::MaterialType MIXTURE::MixtureConstituentElastHyperDamage::material_type() const
{
  return Core::Materials::mix_elasthyper_damage;
}

// Pack the constituent
void MIXTURE::MixtureConstituentElastHyperDamage::pack_constituent(
    Core::Communication::PackBuffer& data) const
{
  MixtureConstituentElastHyperBase::pack_constituent(data);

  Core::Communication::ParObject::add_to_pack(data, current_reference_growth_);
}

// Unpack the constituent
void MIXTURE::MixtureConstituentElastHyperDamage::unpack_constituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureConstituentElastHyperBase::unpack_constituent(position, data);

  Core::Communication::ParObject::extract_from_pack(position, data, current_reference_growth_);
}

// Reads the element from the input file
void MIXTURE::MixtureConstituentElastHyperDamage::read_element(
    int numgp, const Core::IO::InputParameterContainer& container)
{
  MixtureConstituentElastHyperBase::read_element(numgp, container);

  current_reference_growth_.resize(numgp, 1.0);
}

// Updates all summands
void MIXTURE::MixtureConstituentElastHyperDamage::update(Core::LinAlg::Matrix<3, 3> const& defgrd,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  const auto& reference_coordinates = params.get<Core::LinAlg::Matrix<3, 1>>("gp_coords_ref");

  double totaltime = params.get<double>("total time", -1);
  if (totaltime < 0.0)
  {
    FOUR_C_THROW("Parameter 'total time' could not be read!");
  }

  current_reference_growth_[gp] =
      Global::Problem::instance()
          ->function_by_id<Core::UTILS::FunctionOfSpaceTime>(params_->damage_function_id_ - 1)
          .evaluate(reference_coordinates.data(), totaltime, 0);

  MixtureConstituentElastHyperBase::update(defgrd, params, gp, eleGID);
}

double MIXTURE::MixtureConstituentElastHyperDamage::get_growth_scalar(int gp) const
{
  return current_reference_growth_[gp];
}

void MIXTURE::MixtureConstituentElastHyperDamage::evaluate(const Core::LinAlg::Matrix<3, 3>& F,
    const Core::LinAlg::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>& S_stress, Core::LinAlg::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  FOUR_C_THROW("This constituent does not support Evaluation without an elastic part.");
}

void MIXTURE::MixtureConstituentElastHyperDamage::evaluate_elastic_part(
    const Core::LinAlg::Matrix<3, 3>& F, const Core::LinAlg::Matrix<3, 3>& iFextin,
    Teuchos::ParameterList& params, Core::LinAlg::Matrix<6, 1>& S_stress,
    Core::LinAlg::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  // Compute total inelastic deformation gradient
  static Core::LinAlg::Matrix<3, 3> iFin(false);
  iFin.multiply_nn(iFextin, prestretch_tensor(gp));

  // Evaluate 3D elastic part
  Mat::ElastHyperEvaluateElasticPart(
      F, iFin, S_stress, cmat, summands(), summand_properties(), gp, eleGID);
}
FOUR_C_NAMESPACE_CLOSE
