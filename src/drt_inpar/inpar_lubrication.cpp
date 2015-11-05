/*--------------------------------------------------------------------------*/
/*!
\file inpar_lubrication.cpp

\brief Lubrication dynamic parameters

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/


#include "drt_validparameters.H"

#include "inpar_lubrication.H"

void INPAR::LUBRICATION::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
//  using Teuchos::tuple;
//  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& reynolsdyn = list->sublist(
      "LUBRICATION DYNAMIC",
      false,
      "control parameters for Lubrication problems\n");

//  setStringToIntegralParameter<int>("SOLVERTYPE","linear_full",
//                               "type of scalar transport solver",
//                               tuple<std::string>(
//                                 "linear_full",
//                                 "linear_incremental",
//                                 "nonlinear"
//                                 ),
//                               tuple<int>(
//                                   solvertype_linear_full,
//                                   solvertype_linear_incremental,
//                                   solvertype_nonlinear),
//                               &reynolsdyn);
//
//  setStringToIntegralParameter<int>("TIMEINTEGR","One_Step_Theta",
//                               "Time Integration Scheme",
//                               tuple<std::string>(
//                                 "Stationary",
//                                 "One_Step_Theta",
//                                 "BDF2",
//                                 "Gen_Alpha"
//                                 ),
//                               tuple<int>(
//                                   timeint_stationary,
//                                   timeint_one_step_theta,
//                                   timeint_bdf2,
//                                   timeint_gen_alpha
//                                 ),
//                               &reynolsdyn);

  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&reynolsdyn);
  IntParameter("NUMSTEP",20,"Total number of time steps",&reynolsdyn);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&reynolsdyn);
//  DoubleParameter("THETA",0.5,"One-step-theta time integration factor",&reynolsdyn);
//  DoubleParameter("ALPHA_M",0.5,"Generalized-alpha time integration factor",&reynolsdyn);
//  DoubleParameter("ALPHA_F",0.5,"Generalized-alpha time integration factor",&reynolsdyn);
//  DoubleParameter("GAMMA",0.5,"Generalized-alpha time integration factor",&reynolsdyn);
  //IntParameter("WRITESOLEVRY",1,"Increment for writing solution",&reynolsdyn);
  IntParameter("UPRES",1,"Increment for writing solution",&reynolsdyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&reynolsdyn);
//  IntParameter("MATID",-1,"Material ID for automatic mesh generation",&reynolsdyn);
//
//  setStringToIntegralParameter<int>("VELOCITYFIELD","zero",
//                               "type of velocity field used for scalar transport problems",
//                               tuple<std::string>(
//                                 "zero",
//                                 "function",
//                                 "function_and_curve",
//                                 "Navier_Stokes"
//                                 ),
//                               tuple<int>(
//                                   velocity_zero,
//                                   velocity_function,
//                                   velocity_function_and_curve,
//                                   velocity_Navier_Stokes),
//                               &reynolsdyn);
//
//  IntParameter("VELFUNCNO",-1,"function number for scalar transport velocity field",&reynolsdyn);
//
//  IntParameter("VELCURVENO",-1,"curve number for time-dependent scalar transport velocity field",&reynolsdyn);
//
//  {
//    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
//    Teuchos::Tuple<std::string,12> name;
//    Teuchos::Tuple<int,12> label;
//    name[ 0] = "zero_field";                   label[ 0] = initfield_zero_field;
//    name[ 1] = "field_by_function";            label[ 1] = initfield_field_by_function;
//    name[ 2] = "field_by_condition";           label[ 2] = initfield_field_by_condition;
//    name[ 3] = "disturbed_field_by_function";  label[ 3] = initfield_disturbed_field_by_function;
//    name[ 4] = "1D_DISCONTPV";                 label[ 4] = initfield_discontprogvar_1D;
//    name[ 5] = "FLAME_VORTEX_INTERACTION";     label[ 5] = initfield_flame_vortex_interaction;
//    name[ 6] = "RAYTAYMIXFRAC";                label[ 6] = initfield_raytaymixfrac;
//    name[ 7] = "L_shaped_domain";              label[ 7] = initfield_Lshapeddomain;
//    name[ 8] = "facing_flame_fronts";          label[ 8] = initfield_facing_flame_fronts;
//    name[ 9] = "oracles_flame";                label[ 9] = initfield_oracles_flame;
//    name[10] = "high_forced_hit";              label[10] = initialfield_forced_hit_high_Sc;
//    name[11] = "low_forced_hit";               label[11] = initialfield_forced_hit_low_Sc;
//
//    setStringToIntegralParameter<int>(
//        "INITIALFIELD",
//        "zero_field",
//        "Initial Field for scalar transport problem",
//        name,
//        label,
//        &reynolsdyn);
//  }
//
//  IntParameter("INITFUNCNO",-1,"function number for scalar transport initial field",&reynolsdyn);
//
//  setStringToIntegralParameter<int>("CALCERROR","No",
//                               "compute error compared to analytical solution",
//                               tuple<std::string>(
//                                 "No",
//                                 "Kwok_Wu",
//                                 "ConcentricCylinders",
//                                 "Electroneutrality",
//                                 "error_by_function",
//                                 "SphereDiffusion"
//                                 ),
//                               tuple<int>(
//                                   calcerror_no,
//                                   calcerror_Kwok_Wu,
//                                   calcerror_cylinder,
//                                   calcerror_electroneutrality,
//                                   calcerror_byfunction,
//                                   calcerror_spherediffusion
//                                   ),
//                               &reynolsdyn);
//
//  IntParameter("CALCERRORNO",-1,"function number for scalar transport error computation",&reynolsdyn);
//
//  setStringToIntegralParameter<int>("WRITEFLUX","No","output of diffusive/total flux vectors",
//                               tuple<std::string>(
//                                 "No",
//                                 "totalflux_domain",
//                                 "diffusiveflux_domain",
//                                 "totalflux_boundary",
//                                 "diffusiveflux_boundary",
//                                 "convectiveflux_boundary"
//                                 ),
//                               tuple<int>(
//                                   flux_no,
//                                   flux_total_domain,
//                                   flux_diffusive_domain,
//                                   flux_total_boundary,
//                                   flux_diffusive_boundary,
//                                   flux_convective_boundary),
//                               &reynolsdyn);
//
//  // Parameters for reaction-diffusion systems (for example cardiac electrophysiology)
//  IntParameter("WRITEMAXINTSTATE",0,"number of maximal internal state variables to be postprocessed",&reynolsdyn);
//  IntParameter("WRITEMAXIONICCURRENTS",0,"number of maximal ionic currents to be postprocessed",&reynolsdyn);
//  DoubleParameter("ACTTHRES",1.0,"threshold for the potential for computing and postprocessing activation time ",&reynolsdyn);
//
//  setNumericStringParameter("WRITEFLUX_IDS","-1",
//      "Write diffusive/total flux vector fields for these scalar fields only (starting with 1)",
//      &reynolsdyn);
//
//  BoolParameter("OUTMEAN","No","Output of mean values for scalars and density",&reynolsdyn);
//  BoolParameter("OUTINTEGRREAC","No","Output of integral reaction values",&reynolsdyn);
//  BoolParameter("OUTPUT_GMSH","No","Do you want to write Gmsh postprocessing files?",&reynolsdyn);
//
//  BoolParameter("MATLAB_STATE_OUTPUT","No","Do you want to write the state solution to Matlab file?",&reynolsdyn);
//
//  setStringToIntegralParameter<int>("CONVFORM","convective","form of convective term",
//                               tuple<std::string>(
//                                 "convective",
//                                 "conservative"
//                                 ),
//                               tuple<int>(
//                                 convform_convective,
//                                 convform_conservative),
//                               &reynolsdyn);
//
//  BoolParameter("NEUMANNINFLOW",
//      "no","Flag to (de)activate potential Neumann inflow term(s)",&reynolsdyn);
//
//  BoolParameter("CONV_HEAT_TRANS",
//      "no","Flag to (de)activate potential convective heat transfer boundary conditions",&reynolsdyn);
//
//  BoolParameter("SKIPINITDER",
//      "no","Flag to skip computation of initial time derivative",&reynolsdyn);
//
//  setStringToIntegralParameter<int>("FSSUGRDIFF",
//                               "No",
//                               "fine-scale subgrid diffusivity",
//                               tuple<std::string>(
//                                 "No",
//                                 "artificial",
//                                 "Smagorinsky_all",
//                                 "Smagorinsky_small"
//                                 ),
//                               tuple<int>(
//                                   fssugrdiff_no,
//                                   fssugrdiff_artificial,
//                                   fssugrdiff_smagorinsky_all,
//                                   fssugrdiff_smagorinsky_small),
//                               &reynolsdyn);
//
//  setStringToIntegralParameter<int>("MESHTYING", "no", "Flag to (de)activate mesh tying algorithm",
//                                  tuple<std::string>(
//                                      "no",
//                                      "Condensed_Smat",
//                                      "Condensed_Bmat",
//                                      "Condensed_Bmat_merged"), //use the condensed_bmat_merged strategy
//                                    tuple<int>(
//                                        INPAR::FLUID::no_meshtying,
//                                        INPAR::FLUID::condensed_smat,
//                                        INPAR::FLUID::condensed_bmat,
//                                        INPAR::FLUID::condensed_bmat_merged),   //use the condensed_bmat_merged strategy
//                                    &reynolsdyn);

  // linear solver id used for scalar transport/elch problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for the Lubrication problem",&reynolsdyn);
  //IntParameter("SIMPLER_SOLVER",-1,"number of linear solver used for ELCH (solved with SIMPLER)...",&reynolsdyn);

  // parameters for natural convection effects
//  BoolParameter("NATURAL_CONVECTION","No","Include natural convection effects",&reynolsdyn);
//  IntParameter("NATCONVITEMAX",10,"Maximum number of outer iterations for natural convection",&reynolsdyn);
//  DoubleParameter("NATCONVCONVTOL",1e-6,"Convergence check tolerance for outer loop for natural convection",&reynolsdyn);
//
//  // parameters for finite difference check
//  setStringToIntegralParameter<int>("FDCHECK", "none", "flag for finite difference check: none, local, or global",
//                                    tuple<std::string>(
//                                      "none",
//                                      "local",    // perform finite difference check on element level
//                                      "global"),  // perform finite difference check on time integrator level
//                                    tuple<int>(
//                                        fdcheck_none,
//                                        fdcheck_local,
//                                        fdcheck_global),
//                                    &reynolsdyn);
//  DoubleParameter("FDCHECKEPS",1.e-6,"dof perturbation magnitude for finite difference check (1.e-6 seems to work very well, whereas smaller values don't)",&reynolsdyn);
//  DoubleParameter("FDCHECKTOL",1.e-6,"relative tolerance for finite difference check",&reynolsdyn);
//
//  // parameter for optional computation of domain and boundary integrals, i.e., of surface areas and volumes associated with specified nodesets
//  setStringToIntegralParameter<int>(
//      "COMPUTEINTEGRALS",
//      "none",
//      "flag for optional computation of domain integrals",
//      tuple<std::string>(
//          "none",
//          "initial",
//          "repeated"
//          ),
//      tuple<int>(
//          computeintegrals_none,
//          computeintegrals_initial,
//          computeintegrals_repeated
//          ),
//      &reynolsdyn
//      );

}
