/*----------------------------------------------------------------------*/
/*!
\file inpar_cavitation.cpp

\brief Input parameters for cavitation problems

\level 1

\maintainer Georg Hammerl
*-----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_cavitation.H"
#include "inpar_parameterlist_utils.H"
#include "../drt_lib/drt_conditiondefinition.H"


void INPAR::CAVITATION::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& cavitationdyn = list->sublist(
      "CAVITATION DYNAMIC",
      false,
      "control parameters for cavitation problems\n");

  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&cavitationdyn);
  IntParameter("NUMSTEP",20,"Total number of time steps",&cavitationdyn);
  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&cavitationdyn);
  IntParameter("RESULTSEVRY",1,"Increment for writing solution",&cavitationdyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&cavitationdyn);
  IntParameter("RESTARTSTEP_PARTICLES",-1,"Step for restarting particles in cavitation problem",&cavitationdyn);
  IntParameter("ITEMAX",500,"maximum number of iterations over fields",&cavitationdyn);
  DoubleParameter("CONVTOL",1e-5,"tolerance for convergence check of iteratively coupled cavitation problems",&cavitationdyn);

  // Coupling strategy
  setStringToIntegralParameter<int>(
                              "COUPALGO","cavitation_twowaymomentum",
                              "Coupling strategies for cavitation problems",
                              tuple<std::string>(
                                "cavitation_oneway",
                                "cavitation_twowaymomentum",
                                "cavitation_voidfrac_only",
                                "cavitation_twowayfull_weak",
                                "cavitation_twowayfull_strong"
                                ),
                              tuple<int>(
                                OneWay,
                                TwoWayMomentum,
                                VoidFracOnly,
                                TwoWayFull_weak,
                                TwoWayFull_strong
                                ),
                              &cavitationdyn);
  setStringToIntegralParameter<int>(
                              "VOID_FRACTION_CALC","analytical_quadraticpoly",
                              "Void fraction calculation strategy",
                              tuple<std::string>(
                                "analytical_constpoly",
                                "analytical_quadraticpoly",
                                "analytical_quarticpoly",
                                "gaussian_integration"
                                ),
                              tuple<int>(
                                analytical_constpoly,
                                analytical_quadraticpoly,
                                analytical_quarticpoly,
                                gaussian_integration
                                ),
                              &cavitationdyn);
  IntParameter("NUM_GP_VOID_FRACTION",4,"Number of gauss points in each direction for void fraction calculation",&cavitationdyn);

  BoolParameter("APPROX_ELECOORDS_INIT","no","switch on/off approximate initial guess for computing element coordinates",&cavitationdyn);

  BoolParameter("SIMPLIFIED_BUBBLE_FORCES","no","switch on/off simplified bubble force computation",&cavitationdyn);

  setStringToIntegralParameter<int>("FLUIDFRAC_PROJ_METHOD", "superconvergent_patch_recovery",
                               "Flag to decide which reconstruction type for the fluid fraction is chosen.",
                               tuple<std::string>(
                                 "superconvergent_patch_recovery",
                                 "L2_projection"),
                               tuple<std::string>(
                                 "fluid fraction reconstruction via superconvergent patch recovery",
                                 "fluid fraction reconstruction via L2-projection"),
                               tuple<int>(
                                 fluidfracreco_spr,
                                 fluidfracreco_l2
                               ),
                               &cavitationdyn);

  IntParameter("FLUIDFRAC_PROJ_SOLVER",-1,"Number of linear solver used for L2 projection",&cavitationdyn);

  IntParameter("TIME_STEP_SIZE_RATIO",1,"Ration between fluid and particle time step size in cavitation problems",&cavitationdyn);

  BoolParameter("INFLOW_RADIUS_BLENDING","yes","switch on/off blending of radius of inflowing particles",&cavitationdyn);

  BoolParameter("COMPUTE_RADIUS_RP_BASED","no","switch on/off radius calculation based on Ralyeigh-Plesset equation",&cavitationdyn);

  BoolParameter("INIT_BUBBLEVEL_FROM_FLUID","no","interpolate initial velocity for particles from fluid field",&cavitationdyn);

  DoubleParameter("INFLUENCE_SCALING",1.2,"Scale for bubble radius to obtain influence radius for void frac computation",&cavitationdyn);
}
