// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_CONTACT_BEAM_TO_BEAM_INPUT_HPP
#define FOUR_C_BEAMINTERACTION_CONTACT_BEAM_TO_BEAM_INPUT_HPP

#include "4C_config.hpp"

#include "4C_fem_condition.hpp"
#include "4C_io_input_spec.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Conditions
{
  class ConditionDefinition;
}  // namespace Core::Conditions

namespace BeamInteraction::Contact::BeamToBeam
{

  enum PenaltyLaw
  {
    pl_lp,     ///< linear penalty law
    pl_qp,     ///< quadratic penalty law
    pl_lnqp,   ///< linear penalty law with quadratic regularization for negative gaps
    pl_lpqp,   ///< linear penalty law with quadratic regularization for positive gaps
    pl_lpcp,   ///< linear penalty law with cubic regularization for positive gaps
    pl_lpdqp,  ///< linear penalty law with double quadratic regularization for positive gaps
    pl_lpep    ///< linear penalty law with exponential regularization for positive gaps
  };

  std::vector<Core::IO::InputSpec> valid_parameters();

  void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);

}  // namespace BeamInteraction::Contact::BeamToBeam

FOUR_C_NAMESPACE_CLOSE

#endif
