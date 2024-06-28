/*---------------------------------------------------------------------*/
/*! \file

\brief options to set up the cut

\level 2


*----------------------------------------------------------------------*/

#include "4C_cut_options.hpp"

#include "4C_cut_position.hpp"
#include "4C_global_data.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/// Initializes Cut Parameters by Parameterlist (typically from *.dat-file section CUT GENERAL)
void Core::Geo::Cut::Options::Init_by_Paramlist()
{
  Init_by_Paramlist(Global::Problem::Instance()->CutGeneralParams());
}

/// Initializes Cut Parameters by Parameterlist (typically from *.dat-file section CUT GENERAL)
void Core::Geo::Cut::Options::Init_by_Paramlist(const Teuchos::ParameterList& cutparams)
{
  geomintersect_floattype_ =
      Core::UTILS::IntegralValue<CutFloatType>(cutparams, "KERNEL_INTERSECTION_FLOATTYPE");
  geomdistance_floattype_ =
      Core::UTILS::IntegralValue<CutFloatType>(cutparams, "KERNEL_DISTANCE_FLOATTYPE");
  general_position_dist_floattype_ =
      Core::UTILS::IntegralValue<CutFloatType>(cutparams, "GENERAL_POSITON_DISTANCE_FLOATTYPE");
  general_position_pos_floattype_ =
      Core::UTILS::IntegralValue<CutFloatType>(cutparams, "GENERAL_POSITON_POSITION_FLOATTYPE");
  direct_divergence_refplane_ = Core::UTILS::IntegralValue<CutDirectDivergenceRefplane>(
      cutparams, "DIRECT_DIVERGENCE_REFPLANE");
  Core::Geo::Cut::PositionFactory::specify_general_dist_floattype(general_position_dist_floattype_);
  Core::Geo::Cut::PositionFactory::specify_general_pos_floattype(general_position_pos_floattype_);
  split_cutsides_ = Core::UTILS::IntegralValue<bool>(cutparams, "SPLIT_CUTSIDES");
  do_selfcut_ = Core::UTILS::IntegralValue<bool>(cutparams, "DO_SELFCUT");
  selfcut_do_meshcorrection_ =
      Core::UTILS::IntegralValue<bool>(cutparams, "SELFCUT_DO_MESHCORRECTION");
  selfcut_island_geom_multiplicator_ = cutparams.get<int>("SELFCUT_MESHCORRECTION_MULTIPLICATOR");
  bc_cubaturedegree_ = cutparams.get<int>("BOUNDARYCELL_CUBATURDEGREE");
}

/// Initializes Cut Parameters for Cuttests (use full cln) -- slowest option
void Core::Geo::Cut::Options::Init_for_Cuttests()
{
  geomintersect_floattype_ = floattype_cln;
  geomdistance_floattype_ = floattype_cln;
  general_position_dist_floattype_ = floattype_cln;
  general_position_pos_floattype_ = floattype_cln;
  direct_divergence_refplane_ = DirDiv_refplane_all;
  Core::Geo::Cut::PositionFactory::specify_general_dist_floattype(general_position_dist_floattype_);
  Core::Geo::Cut::PositionFactory::specify_general_pos_floattype(general_position_pos_floattype_);
  split_cutsides_ = true;
  selfcut_do_meshcorrection_ = false;
  bc_cubaturedegree_ = 20;
}

FOUR_C_NAMESPACE_CLOSE
