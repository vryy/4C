// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MIXTURE_RULE_MAP_HPP
#define FOUR_C_MIXTURE_RULE_MAP_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mixture_rule.hpp"

#include <memory>
#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Mixture
{
  // forward declaration
  class MapMixtureRule;

  namespace PAR
  {
    class MapMixtureRule : public Mixture::PAR::MixtureRule
    {
      friend class Mixture::MapMixtureRule;

     public:
      /// constructor
      explicit MapMixtureRule(const Core::Mat::PAR::Parameter::Data& matdata);

      /// Create mixturerule instance
      std::unique_ptr<Mixture::MixtureRule> create_rule() override;

      /// @name parameters of the mixture rule
      /// @{
      const double initial_reference_density_;
      const int num_constituents_;
      const std::unordered_map<int, std::vector<double>> mass_fractions_map_;
      /// @}
    };

  }  // namespace PAR

  /*!
   * \brief Mixture rule to be used in homogenized constrained mixture models. It scales the mass
   * fractions of the individual constituents by elementwise defined discrete values, that are
   * passed via an external '*.pattern' file.
   */
  class MapMixtureRule : public Mixture::MixtureRule
  {
   public:
    /// Constructor for mixture rule given the input parameters
    explicit MapMixtureRule(Mixture::PAR::MapMixtureRule* params);

    void evaluate(const Core::LinAlg::Matrix<3, 3>& F, const Core::LinAlg::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, Core::LinAlg::Matrix<6, 1>& S_stress,
        Core::LinAlg::Matrix<6, 6>& cmat, int gp, int eleGID) override;

    [[nodiscard]] double return_mass_density() const override
    {
      return params_->initial_reference_density_;
    };

    void setup(Teuchos::ParameterList& params, const int eleGID) override;

    void unpack_mixture_rule(Core::Communication::UnpackBuffer& buffer) override;

   private:
    ///! Rule parameters as defined in the input file
    PAR::MapMixtureRule* params_{};
  };
}  // namespace Mixture

FOUR_C_NAMESPACE_CLOSE

#endif