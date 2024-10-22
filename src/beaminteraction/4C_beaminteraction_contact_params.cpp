// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_contact_params.hpp"

#include "4C_beaminteraction_beam_to_beam_contact_params.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_contact_params.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_meshtying_params.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_params.hpp"
#include "4C_beaminteraction_beam_to_sphere_contact_params.hpp"
#include "4C_beaminteraction_contact_runtime_visualization_output_params.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamContactParams::BeamContactParams()
    : beam_to_beam_contact_params_(Teuchos::null),
      beam_to_sphere_contact_params_(Teuchos::null),
      beam_to_solid_volume_meshtying_params_(Teuchos::null),
      beam_to_solid_surface_meshtying_params_(Teuchos::null),
      beam_to_solid_surface_contact_params_(Teuchos::null),
      beam_contact_runtime_output_params_(Teuchos::null)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactParams::build_beam_to_beam_contact_params()
{
  beam_to_beam_contact_params_ = Teuchos::make_rcp<BEAMINTERACTION::BeamToBeamContactParams>();
  beam_to_beam_contact_params_->init();
  beam_to_beam_contact_params_->setup();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactParams::build_beam_contact_runtime_output_params(
    const double restart_time)
{
  beam_contact_runtime_output_params_ =
      Teuchos::make_rcp<BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams>(restart_time);
  beam_contact_runtime_output_params_->init();
  beam_contact_runtime_output_params_->setup();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactParams::build_beam_to_sphere_contact_params()
{
  beam_to_sphere_contact_params_ = Teuchos::make_rcp<BEAMINTERACTION::BeamToSphereContactParams>();
  beam_to_sphere_contact_params_->init();
  beam_to_sphere_contact_params_->setup();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactParams::build_beam_to_solid_volume_meshtying_params()
{
  beam_to_solid_volume_meshtying_params_ =
      Teuchos::make_rcp<BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams>();
  beam_to_solid_volume_meshtying_params_->init();
  beam_to_solid_volume_meshtying_params_->setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactParams::build_beam_to_solid_surface_meshtying_params()
{
  beam_to_solid_surface_meshtying_params_ =
      Teuchos::make_rcp<BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams>();
  beam_to_solid_surface_meshtying_params_->init();
  beam_to_solid_surface_meshtying_params_->setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactParams::build_beam_to_solid_surface_contact_params()
{
  beam_to_solid_surface_contact_params_ =
      Teuchos::make_rcp<BEAMINTERACTION::BeamToSolidSurfaceContactParams>();
  beam_to_solid_surface_contact_params_->init();
  beam_to_solid_surface_contact_params_->setup();
}

FOUR_C_NAMESPACE_CLOSE
