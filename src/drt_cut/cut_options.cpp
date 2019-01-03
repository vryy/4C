/*---------------------------------------------------------------------*/
/*!
\file cut_options.cpp

\brief options to set up the cut

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#include "cut_options.H"

#include "cut_position.H"
#include "../drt_lib/drt_globalproblem.H"

/// Initializes Cut Parameters by Parameterlist (typically from *.dat-file section CUT GENERAL)
void GEO::CUT::Options::Init_by_Paramlist()
{
  Init_by_Paramlist(DRT::Problem::Instance()->CutGeneralParams());
}

/// Initializes Cut Parameters by Parameterlist (typically from *.dat-file section CUT GENERAL)
void GEO::CUT::Options::Init_by_Paramlist(const Teuchos::ParameterList& cutparams)
{
  geomintersect_floattype_ = DRT::INPUT::IntegralValue<INPAR::CUT::CUT_Floattype>(
      cutparams, "KERNEL_INTERSECTION_FLOATTYPE");
  geomdistance_floattype_ =
      DRT::INPUT::IntegralValue<INPAR::CUT::CUT_Floattype>(cutparams, "KERNEL_DISTANCE_FLOATTYPE");
  general_position_dist_floattype_ = DRT::INPUT::IntegralValue<INPAR::CUT::CUT_Floattype>(
      cutparams, "GENERAL_POSITON_DISTANCE_FLOATTYPE");
  general_position_pos_floattype_ = DRT::INPUT::IntegralValue<INPAR::CUT::CUT_Floattype>(
      cutparams, "GENERAL_POSITON_POSITION_FLOATTYPE");
  GEO::CUT::PositionFactory::SpecifyGeneralDistFloattype(general_position_dist_floattype_);
  GEO::CUT::PositionFactory::SpecifyGeneralPosFloattype(general_position_pos_floattype_);
}

/// Initializes Cut Parameters for Cuttests (use full cln) -- slowest option
void GEO::CUT::Options::Init_for_Cuttests()
{
  geomintersect_floattype_ = INPAR::CUT::floattype_cln;
  geomdistance_floattype_ = INPAR::CUT::floattype_cln;
  general_position_dist_floattype_ = INPAR::CUT::floattype_cln;
  general_position_pos_floattype_ = INPAR::CUT::floattype_double;
  GEO::CUT::PositionFactory::SpecifyGeneralDistFloattype(general_position_dist_floattype_);
  GEO::CUT::PositionFactory::SpecifyGeneralPosFloattype(general_position_pos_floattype_);
}
