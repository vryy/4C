/*----------------------------------------------------------------------*/
/*!
\file inpar_cell.cpp

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_cell.H"



void INPAR::CELL::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& celldyn = list->sublist("CELL DYNAMIC",false,
                                                  "Cell Migration Simulation");

  Teuchos::Tuple<std::string,3> coupname;
  Teuchos::Tuple<int,3> couplabel;

  coupname[ 0] = "basic_sequ_stagg";                              couplabel[ 0] = cell_basic_sequ_stagg;
  coupname[ 1] = "iter_stagg_fixed_rel_param";                    couplabel[ 1] = cell_iter_stagg_fixed_rel_param;
  coupname[ 2] = "iter_stagg_AITKEN_rel_param";                   couplabel[ 2] = cell_iter_stagg_AITKEN_rel_param;

  setStringToIntegralParameter<int>("COUPALGO","iter_stagg_AITKEN_rel_param",
                                    "Iteration Scheme over the fields",
                                    coupname,
                                    couplabel,
                                    &celldyn);


  setStringToIntegralParameter<int>(
                                    "PSEUDO2D","no",
                                    "True if a thin quasi 2-dimensional problem is modeled",
                                    tuple<std::string>(
                                                       "no",
                                                       "yes"),
                                    tuple<int>(
                                               0,
                                               1),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "SIMTYPE","pureFSI",
                                    "Simulation Type",
                                    tuple<std::string>(
                                                       "pureFSI",
                                                       "pureAdhesion",
                                                       "pureCompression",
                                                       "pureGrowth",
                                                       "Multiphysics"),
                                    tuple<int>(
                                               sim_type_pureFSI,
                                               sim_type_pureAdhesion,
                                               sim_type_pureCompression,
                                               sim_type_pureGrowth,
                                               sim_type_multiphysics),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "MIGRATIONTYPE","undefined",
                                    "Migration with or without ScaTra.",
                                    tuple<std::string>(
                                                       "undefined",
                                                       "ameboid",
                                                       "proteolytic"),
                                    tuple<int>(
                                               cell_migration_undefined,
                                               cell_migration_ameboid,
                                               cell_migration_proteolytic),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "SEGREGATION","undefined",
                                    "Segregation of chemical species by cell in volume or on surface",
                                    tuple<std::string>(
                                                       "undefined",
                                                       "volume",
                                                       "surface"),
                                    tuple<int>(
                                               segregation_undefined,
                                               segregation_volumetric,
                                               segregation_surface),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "SEGREGATION_BY","dirichlet",
                                    "Segregation via dirichlet or neumann condition",
                                    tuple<std::string>(
                                                       "undefined",
                                                       "dirichlet",
                                                       "neumann"),
                                    tuple<int>(
                                               segregation_by_undefined,
                                               segregation_by_dirichlet,
                                               segregation_by_neumann),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "SEGREGATION_LAW","undefined",
                                    "Segregation of chemical species by cell in volume or on surface",
                                    tuple<std::string>(
                                                       "undefined",
                                                       "constant",
                                                       "forcedependent"),
                                    tuple<int>(
                                               segregation_law_undefined,
                                               segregation_law_constant,
                                               segregation_law_forcedependent),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "INITIALIZE_CELL","no",
                                    "Use first time step as pre-simulation to initialize Cell",
                                    tuple<std::string>(
                                                       "no",
                                                       "yes"),
                                    tuple<int>(
                                               0,
                                               1),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "ECM_INTERACTION","no",
                                    "Switch on/off cell-ecm interaction model",
                                    tuple<std::string>(
                                                       "no",
                                                       "yes"),
                                    tuple<int>(
                                               0,
                                               1),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "FLUID_INTERACTION","yes",
                                    "Switch on/off cell-fluid interaction model",
                                    tuple<std::string>(
                                                       "no",
                                                       "yes"),
                                    tuple<int>(
                                               0,
                                               1),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "ADHESION_DYNAMICS","no",
                                    "include adhesion dynamics in simulation",
                                    tuple<std::string>(
                                                       "no",
                                                       "yes"),
                                    tuple<int>(
                                               0,
                                               1),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "SSI_CELL","no",
                                    "build cell as ssi problem to enable intracellular and membrane transport",
                                    tuple<std::string>(
                                                       "no",
                                                       "yes"),
                                    tuple<int>(
                                               0,
                                               1),
                                    &celldyn);


  DoubleParameter("SEGREGATION_CONST",0.0,"basic constant for segregation equation",&celldyn);

  IntParameter("NUMSTEP",200,"Total number of Timesteps",&celldyn);
  IntParameter("UPRES",1,"Increment for writing solution",&celldyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&celldyn);

  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&celldyn);
  DoubleParameter("INITIAL_TIMESTEP",0.1,"Time increment dt for first time step (pre-simulation)",&celldyn);
  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&celldyn);

  DoubleParameter("PENALTY_INIT",100.0,"Penalty parameter for the celll initialization",&celldyn);
  DoubleParameter("PENALTY_START",100.0,"Penalty parameter at the beginning of each time step",&celldyn);
  DoubleParameter("ECM_FIBER_RADIUS",0.3,"Average radius of ECM fibers",&celldyn);

  IntParameter("INITIALIZATION_STEPS",1,"Num of time steps for initialization (pre-simulation)",&celldyn);
}
