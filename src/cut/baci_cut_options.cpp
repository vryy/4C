/*---------------------------------------------------------------------*/
/*! \file

\brief options to set up the cut

\level 2


*----------------------------------------------------------------------*/

#include "baci_cut_options.hpp"

#include "baci_cut_position.hpp"
#include "baci_global_data.hpp"

FOUR_C_NAMESPACE_OPEN

/// Initializes Cut Parameters by Parameterlist (typically from *.dat-file section CUT GENERAL)
void CORE::GEO::CUT::Options::Init_by_Paramlist()
{
  Init_by_Paramlist(GLOBAL::Problem::Instance()->CutGeneralParams());
}

/// Initializes Cut Parameters by Parameterlist (typically from *.dat-file section CUT GENERAL)
void CORE::GEO::CUT::Options::Init_by_Paramlist(const Teuchos::ParameterList& cutparams)
{
  geomintersect_floattype_ = CORE::UTILS::IntegralValue<INPAR::CUT::CutFloattype>(
      cutparams, "KERNEL_INTERSECTION_FLOATTYPE");
  geomdistance_floattype_ =
      CORE::UTILS::IntegralValue<INPAR::CUT::CutFloattype>(cutparams, "KERNEL_DISTANCE_FLOATTYPE");
  general_position_dist_floattype_ = CORE::UTILS::IntegralValue<INPAR::CUT::CutFloattype>(
      cutparams, "GENERAL_POSITON_DISTANCE_FLOATTYPE");
  general_position_pos_floattype_ = CORE::UTILS::IntegralValue<INPAR::CUT::CutFloattype>(
      cutparams, "GENERAL_POSITON_POSITION_FLOATTYPE");
  direct_divergence_refplane_ = CORE::UTILS::IntegralValue<INPAR::CUT::CutDirectDivergenceRefplane>(
      cutparams, "DIRECT_DIVERGENCE_REFPLANE");
  CORE::GEO::CUT::PositionFactory::SpecifyGeneralDistFloattype(general_position_dist_floattype_);
  CORE::GEO::CUT::PositionFactory::SpecifyGeneralPosFloattype(general_position_pos_floattype_);
  split_cutsides_ = CORE::UTILS::IntegralValue<bool>(cutparams, "SPLIT_CUTSIDES");
  do_selfcut_ = CORE::UTILS::IntegralValue<bool>(cutparams, "DO_SELFCUT");
  selfcut_do_meshcorrection_ =
      CORE::UTILS::IntegralValue<bool>(cutparams, "SELFCUT_DO_MESHCORRECTION");
  selfcut_island_geom_multiplicator_ = cutparams.get<int>("SELFCUT_MESHCORRECTION_MULTIPLICATOR");
  bc_cubaturedegree_ = cutparams.get<int>("BOUNDARYCELL_CUBATURDEGREE");
}

/// Initializes Cut Parameters for Cuttests (use full cln) -- slowest option
void CORE::GEO::CUT::Options::Init_for_Cuttests()
{
  geomintersect_floattype_ = INPAR::CUT::floattype_cln;
  geomdistance_floattype_ = INPAR::CUT::floattype_cln;
  general_position_dist_floattype_ = INPAR::CUT::floattype_cln;
  general_position_pos_floattype_ = INPAR::CUT::floattype_cln;
  direct_divergence_refplane_ = INPAR::CUT::DirDiv_refplane_all;
  CORE::GEO::CUT::PositionFactory::SpecifyGeneralDistFloattype(general_position_dist_floattype_);
  CORE::GEO::CUT::PositionFactory::SpecifyGeneralPosFloattype(general_position_pos_floattype_);
  split_cutsides_ = true;
  selfcut_do_meshcorrection_ = false;
  bc_cubaturedegree_ = 20;
}

FOUR_C_NAMESPACE_CLOSE
