// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_options.hpp"

#include "4C_cut_position.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/// Initializes Cut Parameters by Parameterlist (typically from *.dat-file section CUT GENERAL)
void Cut::Options::init_by_paramlist(const Teuchos::ParameterList& cutparams)
{
  geomintersect_floattype_ =
      Teuchos::getIntegralValue<CutFloatType>(cutparams, "KERNEL_INTERSECTION_FLOATTYPE");
  geomdistance_floattype_ =
      Teuchos::getIntegralValue<CutFloatType>(cutparams, "KERNEL_DISTANCE_FLOATTYPE");
  general_position_dist_floattype_ =
      Teuchos::getIntegralValue<CutFloatType>(cutparams, "GENERAL_POSITION_DISTANCE_FLOATTYPE");
  general_position_pos_floattype_ =
      Teuchos::getIntegralValue<CutFloatType>(cutparams, "GENERAL_POSITION_POSITION_FLOATTYPE");
  direct_divergence_refplane_ = Teuchos::getIntegralValue<CutDirectDivergenceRefplane>(
      cutparams, "DIRECT_DIVERGENCE_REFPLANE");
  Cut::PositionFactory::specify_general_dist_floattype(general_position_dist_floattype_);
  Cut::PositionFactory::specify_general_pos_floattype(general_position_pos_floattype_);
  split_cutsides_ = cutparams.get<bool>("SPLIT_CUTSIDES");
  do_selfcut_ = cutparams.get<bool>("DO_SELFCUT");
  selfcut_do_meshcorrection_ = cutparams.get<bool>("SELFCUT_DO_MESHCORRECTION");
  selfcut_island_geom_multiplicator_ = cutparams.get<int>("SELFCUT_MESHCORRECTION_MULTIPLICATOR");
  bc_cubaturedegree_ = cutparams.get<int>("BOUNDARYCELL_CUBATURDEGREE");
  integrate_inside_cells_ = cutparams.get<bool>("INTEGRATE_INSIDE_CELLS");
}

/// Initializes Cut Parameters for Cuttests (use full cln) -- slowest option
void Cut::Options::init_for_cuttests()
{
  geomintersect_floattype_ = floattype_cln;
  geomdistance_floattype_ = floattype_cln;
  general_position_dist_floattype_ = floattype_cln;
  general_position_pos_floattype_ = floattype_cln;
  direct_divergence_refplane_ = DirDiv_refplane_all;
  Cut::PositionFactory::specify_general_dist_floattype(general_position_dist_floattype_);
  Cut::PositionFactory::specify_general_pos_floattype(general_position_pos_floattype_);
  split_cutsides_ = true;
  selfcut_do_meshcorrection_ = false;
  bc_cubaturedegree_ = 20;
}

FOUR_C_NAMESPACE_CLOSE
