// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_VISCOELAST_CONTRIBUTION_HPP
#define FOUR_C_MAT_VISCOELAST_CONTRIBUTION_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_elasthyper_service.hpp"
#include "4C_mat_viscoelast_state.hpp"

#include <cstddef>
#include <memory>
#include <vector>

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

namespace Mat::ViscoElast
{
  class Summand;

  enum class ViscoModelKind
  {
    iso_rate,
    generalized_maxwell,
    fsls
  };

  struct ContributionPointContext
  {
    int visco_mat_id = -1;
    int gp = -1;
    int ele_gid = -1;
  };

  struct ContributionSetupContext
  {
    ContributionPointContext point;
    const std::vector<std::shared_ptr<Summand>>& visco_summands;
    const std::vector<int>& visco_summand_mat_ids;
    ViscoElastState::ActiveModels active_models;
  };

  struct ContributionUpdateContext
  {
    ContributionPointContext point;
    ViscoElastState& state;
  };

  struct ContributionEvaluateBase
  {
    ContributionPointContext point;
    double dt = 0.0;
    ViscoElastState& state;
    Core::LinAlg::Matrix<6, 1>& stress;
    Core::LinAlg::Matrix<6, 6>& cmat;
  };

  struct IsoRateEvaluateContext
  {
    ContributionEvaluateBase base;
    const Teuchos::ParameterList& params;
    const std::vector<std::shared_ptr<Summand>>& visco_summands;
    const SummandProperties& effective_properties;
    Core::LinAlg::Matrix<6, 1>& c_strain;
    Core::LinAlg::Matrix<6, 1>& c_stress;
    Core::LinAlg::Matrix<6, 1>& i_c_stress;
    Core::LinAlg::Matrix<3, 1>& prinv;
    Core::LinAlg::Matrix<3, 1>& modinv;
    Core::LinAlg::Matrix<7, 1>& rateinv;
    Core::LinAlg::Matrix<7, 1>& modrateinv;
    Core::LinAlg::Matrix<6, 1>& mod_c_strain;
    Core::LinAlg::Matrix<6, 1>& scgrate;
    Core::LinAlg::Matrix<6, 1>& modrcgrate;
    Core::LinAlg::Matrix<8, 1>& mu;
    Core::LinAlg::Matrix<8, 1>& modmu;
    Core::LinAlg::Matrix<33, 1>& xi;
    Core::LinAlg::Matrix<33, 1>& modxi;
    Core::LinAlg::Matrix<6, 1>& id2;
    Core::LinAlg::Matrix<6, 6>& id4;
    Core::LinAlg::Matrix<6, 6>& id4sharp;
  };

  struct GeneralizedMaxwellEvaluateContext
  {
    ContributionEvaluateBase base;
    const Core::LinAlg::Matrix<6, 1>& glstrain_mat;
  };

  struct FslsEvaluateContext
  {
    ContributionEvaluateBase base;
  };

  class Contribution
  {
   public:
    virtual ~Contribution() = default;
    [[nodiscard]] virtual ViscoModelKind kind() const = 0;
    virtual void setup(const ContributionSetupContext& context) = 0;
    virtual void update(const ContributionUpdateContext& context) = 0;
    [[nodiscard]] virtual std::size_t history_entry_count_for_setup() const { return 0; }
    [[nodiscard]] virtual unsigned int history_capacity_for_update() const { return 0; }
  };
}  // namespace Mat::ViscoElast

FOUR_C_NAMESPACE_CLOSE

#endif
