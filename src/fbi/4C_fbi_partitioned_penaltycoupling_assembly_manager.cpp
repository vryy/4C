// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fbi_partitioned_penaltycoupling_assembly_manager.hpp"

#include "4C_beaminteraction_contact_pair.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManager::
    PartitionedBeamInteractionAssemblyManager(
        std::vector<std::shared_ptr<BEAMINTERACTION::BeamContactPair>>& assembly_contact_elepairs)
    : assembly_contact_elepairs_(assembly_contact_elepairs)
{
}

FOUR_C_NAMESPACE_CLOSE
