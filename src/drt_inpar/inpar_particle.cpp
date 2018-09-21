/*---------------------------------------------------------------------------*/
/*!
\file inpar_particle.cpp

\brief input parameters for particle problems

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "inpar_particle.H"
#include "drt_validparameters.H"
#include "inpar_parameterlist_utils.H"

/*---------------------------------------------------------------------------*
 | set the particle parameters                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void INPAR::PARTICLE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /*-------------------------------------------------------------------------*
   | general control parameters for particle simulations                     |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledyn =
      list->sublist("PARTICLE DYNAMIC", false, "control parameters for particle simulations\n");

  // type of particle time integration
  setStringToIntegralParameter<int>("DYNAMICTYP", "VelocityVerlet",
      "type of particle time integration",
      tuple<std::string>("SemiImplicitEuler", "VelocityVerlet"),
      tuple<int>(INPAR::PARTICLE::dyna_semiimpliciteuler, INPAR::PARTICLE::dyna_velocityverlet),
      &particledyn);

  // type of particle interaction
  setStringToIntegralParameter<int>("INTERACTION", "None", "type of particle interaction",
      tuple<std::string>("None"), tuple<int>(INPAR::PARTICLE::interaction_none), &particledyn);

  // output type
  IntParameter(
      "RESULTSEVRY", 1, "write particle runtime output every RESULTSEVRY steps", &particledyn);
  IntParameter("RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &particledyn);

  // data format for written numeric data via vtp
  setStringToIntegralParameter<int>("OUTPUT_DATA_FORMAT", "Binary",
      "data format for written numeric data", tuple<std::string>("Binary", "ASCII"),
      tuple<int>(INPAR::PARTICLE::binary, INPAR::PARTICLE::ascii), &particledyn);

  // write ghosted particles
  BoolParameter(
      "WRITE_GHOSTED_PARTICLES", "no", "write ghosted particles (debug feature)", &particledyn);

  // time loop control
  DoubleParameter("TIMESTEP", 0.01, "time step size", &particledyn);
  IntParameter("NUMSTEP", 100, "maximum number of steps", &particledyn);
  DoubleParameter("MAXTIME", 1.0, "maximum time", &particledyn);

  // gravity acceleration control
  setNumericStringParameter(
      "GRAVITY_ACCELERATION", "0.0 0.0 0.0", "acceleration due to gravity", &particledyn);
  IntParameter("GRAVITY_RAMP_FUNCT", -1, "number of function governing gravity ramp", &particledyn);

  // transfer particles to new bins every time step
  BoolParameter(
      "TRANSFER_EVERY", "no", "transfer particles to new bins every time step", &particledyn);

  // relate particle phase to material id
  StringParameter("PHASE_TO_MATERIAL_ID", "", "relate particle phase to material id", &particledyn);

  /*-------------------------------------------------------------------------*
   | control parameters for initial/boundary conditions                      |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledynconditions =
      particledyn.sublist("INITIAL AND BOUNDARY CONDITIONS", false,
          "control parameters for initial/boundary conditions in particle simulations\n");

  // initial velocity field of particle phase given by function
  StringParameter("INITIAL_VELOCITY_FIELD", "",
      "initial velocity field of particle phase given by function", &particledynconditions);

  // initial acceleration field of particle phase given by function
  StringParameter("INITIAL_ACCELERATION_FIELD", "",
      "initial acceleration field of particle phase given by function", &particledynconditions);

  // dirichlet boundary condition of particle phase given by function
  StringParameter("DIRICHLET_BOUNDARY_CONDITION", "",
      "dirichlet boundary condition of particle phase given by function", &particledynconditions);
}
