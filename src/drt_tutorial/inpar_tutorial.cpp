/*----------------------------------------------------------------------*/
/*!
\file inpar_tutorial.cpp

\brief Input parameters for tutorial

\maintainer Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de

*/

/*----------------------------------------------------------------------*/
#include "inpar_tutorial.H"

#include "../drt_inpar/drt_validparameters.H"


void INPAR::TUTORIAL::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& tutorialparams =
      list->sublist("TUTORIAL DYNAMIC", false, "General parameters for any tutorial");

  setStringToIntegralParameter<int>("TYPE", "NonlinearTruss", "Choose type of tutorial.",
      tuple<std::string>("NonlinearTruss", "FixedPointPartitioned", "StructAleCoupling"),
      tuple<int>(nonlinear_truss, partitioned_fixed_point_scheme, struct_ale_coupling),
      &tutorialparams);

  DoubleParameter("TIMESTEP", 0.1, "Size of time step (n -> n+1)", &tutorialparams);
  DoubleParameter("MAXTIME", 1.0, "Maximum Time of Simulation", &tutorialparams);


  /* Nonlinear Truss Tutorial Sublist*/

  Teuchos::ParameterList& nlntrussdyn =
      tutorialparams.sublist("NONLINEAR TRUSS", false, "Control the nonlinear Truss Tutorial");

  DoubleParameter("CONVTOL_RES", 1e-07, "Convergence Criterion for Residual", &nlntrussdyn);
  DoubleParameter("CONVTOL_INC", 1e-07, "Convergence Criterion for Increment", &nlntrussdyn);


  /* Fixed Point Scheme Tutorial Sublist*/

  Teuchos::ParameterList& fixedpointdyn =
      tutorialparams.sublist("FIXED POINT SCHEME", false, "Control the Fixed-Point Tutorial");

  DoubleParameter("RELAX_PARAMETER", 1.0, "Relaxation Parameter omega", &fixedpointdyn);
  DoubleParameter("CONVTOL", 1e-07, "Convergence Criterion for Increment of Partitioned Scheme",
      &fixedpointdyn);

  return;
}
