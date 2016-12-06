/*----------------------------------------------------------------------*/
/*!
\file inpar_combust.cpp

\brief Input parameters for combustion

\level 3

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_combust.H"



void INPAR::COMBUST::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& combustcontrol = list->sublist("COMBUSTION CONTROL",false,
      "control parameters for a combustion problem");

  DoubleParameter("MAXTIME",10.0,"Total simulation time",&combustcontrol);
  IntParameter("NUMSTEP",100,"Total number of timesteps",&combustcontrol);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&combustcontrol);
  IntParameter("ITEMAX",1,"Total number of FG iterations",&combustcontrol);
  DoubleParameter("CONVTOL",1e-6,"Tolerance for iteration over fields",&combustcontrol);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&combustcontrol);
  IntParameter("RESULTSEVRY",1,"Increment for writing solution",&combustcontrol);

  //Redistribution parameter, redistributes according to "expensive" XFEM elements to processors.
  // Currently out of use!
  DoubleParameter("PARALLEL_REDIST_RATIO_FAC",0.0,"Factor by which the max to min element evaluation time has to increase to trigger a redistribution",&combustcontrol);

  //Parameter to switch off G-function field.
  // Used to precompute a single phase flow solution, e.g. for turbulent problems.
  BoolParameter("GENERATE_FLOW_FIELD","No","Do not solve g-function, since we merely want to set up a fluid field",&combustcontrol);

  //Read no ScaTra restart, as none is existent.
  BoolParameter("RESTART_FROM_FLUID","No","Restart from a standard fluid problem (no scalar transport field). No XFEM dofs allowed!",&combustcontrol);
  //Overwrites ScaTra restart by initial field given in input file.
  BoolParameter("RESTART_SCATRA_INPUT","No","Use ScaTra field from .dat-file instead",&combustcontrol);

  BoolParameter("WRITE_CENTER_OF_MASS","No","write center of mass to file",&combustcontrol);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrolfluid = combustcontrol.sublist("COMBUSTION FLUID",false,
      "control parameters for the fluid field of a combustion problem");

  setStringToIntegralParameter<int>("COMBUSTTYPE","Premixed_Combustion",
      "Type of combustion problem",
      tuple<std::string>(
          "Premixed_Combustion",         //jump-enr vel and press,         with combustion interface cond
          "Two_Phase_Flow",              //kink-enr vel and press
          "Two_Phase_Flow_Surf",         //kink-enr vel and jump-enr press
          "Two_Phase_Flow_Jumps"),       //jump-enr vel and press
          tuple<int>(
              combusttype_premixedcombustion,
              combusttype_twophaseflow,
              combusttype_twophaseflow_surf,
              combusttype_twophaseflowjump),
              &combustcontrolfluid);

  //Face-oriented GP and fluid stab terms for Cut elements.
  // See paper Rasthofer, Schott....
  BoolParameter("XFEMSTABILIZATION","Yes","Switch on/off face integrals based on edge-based stabilization, i.e., ghost penalty terms",&combustcontrolfluid);

  setStringToIntegralParameter<int>("XFEMINTEGRATION","Cut",
      "Type of integration strategy for intersected elements",
      tuple<std::string>(
          "Cut",                 //standard Cut
          "Tetrahedra",          //construct of integration cells on predefined of a hexahedra element (impl by a Student not carefully tested)
          "Hexahedra"),          //subdivides the Cut element (by refinement) into hexahedras to approximate the volumecells (Turned out to be insufficient.)
          tuple<int>(
              xfemintegration_cut,
              xfemintegration_tetrahedra,
              xfemintegration_hexahedra),
              &combustcontrolfluid);

  setStringToIntegralParameter<int>("XFEMTIMEINT","DoNothing",
      "Type of time integration strategy for standard degrees of freedom",
      tuple<std::string>(
          "DoNothing",                                //no standard DoF change
          "SemiLagrange",
          "ExtrapolationOld",                         //@Martin W. ?
          "ExtrapolationNew",
          "MixedSemiLagrangeExtrapolation",
          "MixedSemiLagrangeExtrapolationNew",
          "MixedGhostvalSemiLagrange",
          "MixedGhostvalExtrapolation",
          "MixedGhostvalSemiLagrangeExtrapolation"),
          tuple<int>(
              xfemtimeint_donothing,
              xfemtimeint_semilagrange,
              xfemtimeint_extrapolationold,
              xfemtimeint_extrapolationnew,
              xfemtimeint_mixedSLExtrapol,
              xfemtimeint_mixedSLExtrapolNew,
              xfemtimeint_mixedghostSL,
              xfemtimeint_mixedghostExtrapol,
              xfemtimeint_mixedghostSLExtrapol),
              &combustcontrolfluid);

  setStringToIntegralParameter<int>("XFEMTIMEINT_ENR","DoNothing",
      "Type of time integration strategy for enrichment degrees of freedom",
      tuple<std::string>(
          "DoNothing",
          "QuasiStatic",                          //Keep only the std DoFs
          "Projection",
          "ProjectionScalar",
          "Extrapolation",
          "ExtrapolationScalar"),
          tuple<int>(
              xfemtimeintenr_donothing,
              xfemtimeintenr_quasistatic,
              xfemtimeintenr_project,
              xfemtimeintenr_project_scalar,
              xfemtimeintenr_extrapolate,
              xfemtimeintenr_extrapolate_scalar),
              &combustcontrolfluid);

  setStringToIntegralParameter<int>("XFEMTIMEINT_ENR_COMP","Standard",
      "Type of time integration strategy for enrichment computation",
      tuple<std::string>(
          "Standard",   //@Martin W.
          "Full",
          "Minimal"),
          tuple<int>(
              xfemtimeintenr_standard,
              xfemtimeintenr_full,
              xfemtimeintenr_minimal),
              &combustcontrolfluid);

  setStringToIntegralParameter<int>("INITIALFIELD","zero_field","Initial field for fluid problem",
      tuple<std::string>(
          "zero_field",
          "field_by_function",
          "disturbed_field_by_function",
          "flame_vortex_interaction",
          "darrieus_landau_instability",
          "beltrami_flow"),
          tuple<int>(
              initfield_zero_field,
              initfield_field_by_function,
              initfield_disturbed_field_by_function,
              initfield_flame_vortex_interaction,
              initfield_darrieus_landau_instability,
              initfield_beltrami_flow),
              &combustcontrolfluid);
  setStringToIntegralParameter<int>("NITSCHE_ERROR","nitsche_error_none","To which analyt. solution do we compare?",
      tuple<std::string>(
          "nitsche_error_none",
          "nitsche_error_static_bubble_nxnx1",
          "nitsche_error_static_bubble_nxnxn",
          "nitsche_error_shear",
          "nitsche_error_couette_20x20x1",
          "nitsche_error_straight_bodyforce",
          "nitsche_error_ellipsoid_bubble_2D",
          "nitsche_error_ellipsoid_bubble_3D",
          "nitsche_error_beltrami"),
          tuple<int>(
              nitsche_error_none,
              nitsche_error_static_bubble_nxnx1,
              nitsche_error_static_bubble_nxnxn,
              nitsche_error_shear,
              nitsche_error_couette_20x20x1,
              nitsche_error_straight_bodyforce,
              nitsche_error_ellipsoid_bubble_2D,
              nitsche_error_ellipsoid_bubble_3D,
              nitsche_error_beltrami),
              &combustcontrolfluid);

  setStringToIntegralParameter<int>("SURFTENSAPPROX","surface_tension_approx_none","Type of surface tension approximation",
      tuple<std::string>(
          "surface_tension_approx_none",                         //none
          "surface_tension_approx_fixed_curvature",              //prescribed curvature
          "surface_tension_approx_divgrad",                      //calcs curvature at GP using the smoothed grad_phi and smoothed grad_phi for normal  //Do Not Migrate.
          "surface_tension_approx_divgrad_normal",               //calcs curvature at GP using the smoothed grad_phi and normal on Boundary Cell for normal
          "surface_tension_approx_nodal_curvature",              //calcs curvature at nodes and normal on Boundary Cell for normal
          "surface_tension_approx_laplacebeltrami",              //standard Laplace-Beltrami (see e.g. Fries 2009)
          "surface_tension_approx_laplacebeltrami_smoothed"),    //         Laplace-Beltrami, includes additional projection based on the smoothed normal (see e.g. Gross, Reusken 2009)
          tuple<int>(
              surface_tension_approx_none,
              surface_tension_approx_fixed_curvature,
              surface_tension_approx_divgrad,
              surface_tension_approx_divgrad_normal,
              surface_tension_approx_nodal_curvature,
              surface_tension_approx_laplacebeltrami,
              surface_tension_approx_laplacebeltrami_smoothed),
              &combustcontrolfluid);

  DoubleParameter("VARIABLESURFTENS",0.0,"Variable surface tension coefficient",&combustcontrolfluid);
  setStringToIntegralParameter<int>("SMOOTHGRADPHI","smooth_grad_phi_none","Type of smoothing for grad(phi)",
      tuple<std::string>(
          "smooth_grad_phi_none",
          "smooth_grad_phi_meanvalue",
          "smooth_grad_phi_leastsquares_3D",
          "smooth_grad_phi_leastsquares_2Dx",
          "smooth_grad_phi_leastsquares_2Dy",
          "smooth_grad_phi_leastsquares_2Dz",
          "smooth_grad_phi_l2_projection"),
          tuple<int>(
              smooth_grad_phi_none,
              smooth_grad_phi_meanvalue,
              smooth_grad_phi_leastsquares_3D,
              smooth_grad_phi_leastsquares_2Dx,
              smooth_grad_phi_leastsquares_2Dy,
              smooth_grad_phi_leastsquares_2Dz,
              smooth_grad_phi_l2_projection),
              &combustcontrolfluid);
  // set parameters VELOCITY_JUMP_TYPE and FLUX_JUMP_TYPE in case of CombustType Premixed_Combustion
  // Two_Phase_Flow_Jumps is equal to Premixed_Combustion & vel_jump_none & flux_jump_surface_tension
  setStringToIntegralParameter<int>("VELOCITY_JUMP_TYPE","vel_jump_none","Type of velocity jump",
      tuple<std::string>(
          "vel_jump_none",
          "vel_jump_const",
          "vel_jump_premixed_combustion"),
          tuple<int>(
              vel_jump_none,
              vel_jump_const,
              vel_jump_premixed_combustion),
              &combustcontrolfluid);
  setStringToIntegralParameter<int>("FLUX_JUMP_TYPE","flux_jump_none","Type of flux jump",
      tuple<std::string>(
          "flux_jump_none",
          "flux_jump_const",
          "flux_jump_premixed_combustion",
          "flux_jump_surface_tension"),
          tuple<int>(
              flux_jump_none,
              flux_jump_const,
              flux_jump_premixed_combustion,
              flux_jump_surface_tension),
              &combustcontrolfluid);
  IntParameter("INITFUNCNO",-1,"Function for initial field",&combustcontrolfluid);

  //params @Martin W.
  IntParameter("ITE_MAX_FRS",1,"The maximal number of iterations between fluid and recomputation of reference solution",&combustcontrolfluid);
  DoubleParameter("LAMINAR_FLAMESPEED",1.0,"The laminar flamespeed incorporates all chemical kinetics into the problem for now",&combustcontrolfluid);
  DoubleParameter("MOL_DIFFUSIVITY",0.0,"Molecular diffusivity",&combustcontrolfluid);
  DoubleParameter("MARKSTEIN_LENGTH",0.0,"The Markstein length takes flame curvature into account",&combustcontrolfluid);
  DoubleParameter("NITSCHE_VELOCITY",100.0,"Nitsche parameter to stabilize/penalize the velocity jump",&combustcontrolfluid);
  //Not used anymore...
  DoubleParameter("NITSCHE_PRESSURE",0.0,"Nitsche parameter to stabilize/penalize the pressure jump",&combustcontrolfluid);
  //@ Ursula when back
  setStringToIntegralParameter<int>("NITSCHE_CONVFLUX","Yes","(De)activate Nitsche convective flux term",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("NITSCHE_CONVSTAB","Yes","(De)activate Nitsche convective stabilization term",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("NITSCHE_CONVPENALTY","No","(De)activate Nitsche convective penalty term",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("NITSCHE_MASS","No","(De)activate Nitsche mass conservation term",
                                     yesnotuple,yesnovalue,&combustcontrolfluid); //
  setStringToIntegralParameter<int>("NITSCHE_WEIGHT","intersection_visc_based_harmonic","Definition of weighting",
                                    tuple<std::string>(
                                    "visc_based_harmonic",
                                    "intersection_visc_based_harmonic"),
                                    tuple<int>(
                                    weight_visc_based_harmonic,
                                    weight_intersection_visc_based_harmonic),
                                    &combustcontrolfluid);
  setStringToIntegralParameter<int>("CONNECTED_INTERFACE","No","laplace-beltrami surface-tension evaluation only: consider boundary integrals if interface is not closed",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("SMOOTHED_BOUNDARY_INTEGRATION","No","Turn on/off usage of smoothed normal vectors at interface",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("INITSTATSOL","No","Compute stationary solution as initial solution",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);

  BoolParameter("L2_PROJECTION_SECOND_DERIVATIVES","No","L2 Projection Second Derivatives of Level Set",&combustcontrolfluid);
  setStringToIntegralParameter<int>("NODAL_CURVATURE","l2_projected","Type of calculation of nodal curvature value",
      tuple<std::string>(
          "l2_projected",
          "averaged"),
          tuple<int>(
              l2_projected,
              averaged),
              &combustcontrolfluid);
  DoubleParameter("SMOOTHING_PARAMETER",0.0,"Diffusion Coefficient for Smoothing",&combustcontrolfluid); //Added diffusion for L2_projection

  // for selection of enriched fields (velocity, pressure, velocity+pressure)
  setStringToIntegralParameter<int>("SELECTED_ENRICHMENT","both",
       "select fields which get enriched dofs",
       tuple<std::string>(
           "both",
           "velocity",
           "pressure",
           "none"),
           tuple<int>(
               selectedenrichment_both,
               selectedenrichment_velocity,
               selectedenrichment_pressure,
               selectedenrichment_none),
               &combustcontrolfluid);

  BoolParameter("REPELLANT_FORCE","No","Activate repellant force for turbulent bubbly channel flow",&combustcontrolfluid);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrolgfunc = combustcontrol.sublist("COMBUSTION GFUNCTION",false,
      "control parameters for the G-function (level set) field of a combustion problem");

  setStringToIntegralParameter<int>("REFINEMENT","No","Turn refinement strategy for level set function on/off",
                                     yesnotuple,yesnovalue,&combustcontrolgfunc);
  IntParameter("REFINEMENTLEVEL",-1,"number of refinement level for refinement strategy",&combustcontrolgfunc);
  /*----------------------------------------------------------------------*/
}
