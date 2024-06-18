/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for homogenized constrained mixtures with mass fractions defined as discrete
values per element

\level 3


*/
/*----------------------------------------------------------------------*/
#include "4C_mixture_rule_map.hpp"

#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_constituent.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_ConfigDefs.h>
#include <Teuchos_RCP.hpp>

#include <algorithm>
#include <iosfwd>
#include <memory>
#include <string>
#include <unordered_map>


FOUR_C_NAMESPACE_OPEN

namespace
{
  std::vector<double> GetValidateMassFractions(
      const std::unordered_map<int, std::vector<double>>& mass_fractions_map, const int ele_id_key,
      const std::size_t num_constituents)
  {
    auto it = mass_fractions_map.find(ele_id_key);
    if (it == mass_fractions_map.end())
    {
      FOUR_C_THROW(
          "Element id %d not found in the mass fraction map supplied by csv file.", ele_id_key);
    }

    if (it->second.size() != num_constituents)
    {
      FOUR_C_THROW(
          "Number of mass fractions for element id %d does not match the number of constituents "
          "%d.",
          ele_id_key, num_constituents);
    }
    const std::vector<double> massfracs = it->second;

    // check, whether the mass frac sums up to 1
    const double sum = std::accumulate(massfracs.begin(), massfracs.end(), 0.0);
    if (std::abs(1.0 - sum) > 1e-8)
      FOUR_C_THROW(
          "Mass fractions for element id %d don't sum up to 1, which is unphysical.", ele_id_key);

    return massfracs;
  }
}  // namespace

MIXTURE::PAR::MapMixtureRule::MapMixtureRule(const Core::Mat::PAR::Parameter::Data& matdata)
    : MixtureRule(matdata),
      initial_reference_density_(matdata.parameters.get<double>("DENS")),
      num_constituents_(matdata.parameters.get<int>("NUMCONST")),
      mass_fractions_map_(matdata.parameters.get<std::unordered_map<int, std::vector<double>>>(
          "MASSFRACMAPFILE")){};

std::unique_ptr<MIXTURE::MixtureRule> MIXTURE::PAR::MapMixtureRule::create_rule()
{
  return std::make_unique<MIXTURE::MapMixtureRule>(this);
}

MIXTURE::MapMixtureRule::MapMixtureRule(MIXTURE::PAR::MapMixtureRule* params)
    : MixtureRule(params), params_(params)
{
}

void MIXTURE::MapMixtureRule::setup(Teuchos::ParameterList& params, const int eleGID)
{
  MixtureRule::setup(params, eleGID);
}

void MIXTURE::MapMixtureRule::unpack_mixture_rule(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MIXTURE::MixtureRule::unpack_mixture_rule(position, data);
}

void MIXTURE::MapMixtureRule::evaluate(const Core::LinAlg::Matrix<3, 3>& F,
    const Core::LinAlg::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>& S_stress, Core::LinAlg::Matrix<6, 6>& cmat, const int gp,
    const int eleGID)
{
  // define temporary matrices
  Core::LinAlg::Matrix<6, 1> cstress;
  Core::LinAlg::Matrix<6, 6> ccmat;

  // evaluate the mass fractions at the given element id (one based entires in the csv file)
  auto massfracs =
      GetValidateMassFractions(params_->mass_fractions_map_, eleGID + 1, constituents().size());

  // Iterate over all constituents and add all stress/cmat contributions
  for (std::size_t i = 0; i < constituents().size(); ++i)
  {
    MixtureConstituent& constituent = *constituents()[i];
    cstress.Clear();
    ccmat.Clear();
    constituent.evaluate(F, E_strain, params, cstress, ccmat, gp, eleGID);

    // add stress contribution to global stress
    double constituent_density = params_->initial_reference_density_ * massfracs[i];
    S_stress.Update(constituent_density, cstress, 1.0);
    cmat.Update(constituent_density, ccmat, 1.0);
  }
}


FOUR_C_NAMESPACE_CLOSE
