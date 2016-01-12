/*----------------------------------------------------------------------*/
/*!
\file inpar_immersed.cpp

\brief Input parameters for immersed

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_immersed.H"
#include "inpar_structure.H"
#include "inpar_fsi.H"
#include "../drt_lib/drt_conditiondefinition.H"


void INPAR::IMMERSED::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;
  Teuchos::ParameterList& immersedmethod = list->sublist(
    "IMMERSED METHOD",false,
    "General parameters for any immersed problem"
    );

  setStringToIntegralParameter<int>(
                               "COUPALGO","partitioned",
                               "Coupling strategies for immersed method.",
                               tuple<std::string>(
                                 "partitioned",
                                 "monolithic"),
                                 tuple<int>(
                                 partitioned,
                                 monolithic),
                                 &immersedmethod);

  setStringToIntegralParameter<int>(
                               "SCHEME","dirichletneumann",
                               "Coupling schemes for partitioned immersed method.",
                               tuple<std::string>(
                                 "neumannneumann",
                                 "dirichletneumann"),
                                 tuple<int>(
                                 neumannneumann,
                                 dirichletneumann),
                                 &immersedmethod);


  setStringToIntegralParameter<int>(
                               "PROJECTION","shapefunctions",
                               "Projection of nodal values between the coupled fields.",
                               tuple<std::string>(
                                 "shapefunctions",
                                 "mortar"),
                                 tuple<int>(
                                 shapefunctions,
                                 mortar),
                                 &immersedmethod);

  setStringToIntegralParameter<int>(
                               "APPLY_FORCE_RELAX","globally",
                               "Relax whole force vector or not.",
                               tuple<std::string>(
                                 "globally",
                                 "selectively"),
                                 tuple<int>(
                                 globally,
                                 selectively),
                                 &immersedmethod);

  setStringToIntegralParameter<int>(
                               "APPLY_VEL_RELAX","globally",
                               "Relax whole velocity vector or not.",
                               tuple<std::string>(
                                 "globally",
                                 "selectively"),
                                 tuple<int>(
                                 globally,
                                 selectively),
                                 &immersedmethod);

  setStringToIntegralParameter<int>(
                               "DIVERCONT","stop",
                               "What to do after maxiter is reached.",
                               tuple<std::string>(
                                 "stop",
                                 "continue"),
                                 tuple<int>(
                                 nlnsolver_stop,
                                 nlnsolver_continue),
                                 &immersedmethod);

  setStringToIntegralParameter<int>(
                               "OUTPUT_EVRY_NLNITER","no",
                               "write output after every solution step of the nonlin. part. iter. scheme",
                               tuple<std::string>(
                                 "yes",
                                 "no"),
                                 tuple<int>(
                                 1,
                                 0),
                                 &immersedmethod);

  setStringToIntegralParameter<int>(
                               "CORRECT_BOUNDARY_VELOCITIES","no",
                               "correct velocities in fluid elements cut by surface of immersed structure",
                               tuple<std::string>(
                                 "yes",
                                 "no"),
                                 tuple<int>(
                                 1,
                                 0),
                                 &immersedmethod);

  setStringToIntegralParameter<int>(
                               "TIMESTATS","everyiter",
                               "summarize time monitor every nln iteration",
                               tuple<std::string>(
                                 "everyiter",
                                 "endofsim"),
                                 tuple<int>(
                                 1,
                                 0),
                                 &immersedmethod);

  DoubleParameter("FORCE_RELAX",1.0,"Force Relaxaton Parameter"    ,&immersedmethod);
  DoubleParameter("VEL_RELAX"  ,1.0,"Velocity Relaxation Parameter",&immersedmethod);
  DoubleParameter("FLD_SRCHRADIUS_FAC",1.0,"fac times fluid ele. diag. length",&immersedmethod);
  DoubleParameter("STRCT_SRCHRADIUS_FAC",0.5,"fac times structure bounding box diagonal",&immersedmethod);
  IntParameter("NUM_GP_FLUID_BOUND",8,"number of gp in fluid elements cut by surface of immersed structure (higher number yields better mass conservation)",&immersedmethod);

  /*----------------------------------------------------------------------*/
  /* parameters for paritioned immersed solvers */
  Teuchos::ParameterList& immersedpart = immersedmethod.sublist("PARTITIONED SOLVER",false,"");

  setStringToIntegralParameter<int>(
                                 "PARTITIONED","DirichletNeumann",
                                 "Coupling strategies for partitioned FSI solvers.",
                                 tuple<std::string>(
                                   "DirichletNeumann",
                                   "DirichletNeumannSlideALE"
                                   ),
                                 tuple<int>(
                                   INPAR::FSI::DirichletNeumann,
                                   INPAR::FSI::DirichletNeumannSlideale
                                   ),
                                 &immersedpart);

  setStringToIntegralParameter<int>("PREDICTOR","d(n)",
                                 "Predictor for interface displacements",
                                 tuple<std::string>(
                                   "d(n)",
                                   "d(n)+dt*(1.5*v(n)-0.5*v(n-1))",
                                   "d(n)+dt*v(n)",
                                   "d(n)+dt*v(n)+0.5*dt^2*a(n)"
                                   ),
                                 tuple<int>(1,2,3,4),
                                 &immersedpart);

    setStringToIntegralParameter<int>("COUPVARIABLE_FSI","Displacement",
                                 "Coupling variable at the fsi interface",
                                 tuple<std::string>("Displacement","Force"),
                                 tuple<int>(0,1),
                                 &immersedpart);

    setStringToIntegralParameter<int>("COUPVARIABLE_ADHESION","Displacement",
                                 "Coupling variable at the adhesion interface",
                                 tuple<std::string>("Displacement","Force"),
                                 tuple<int>(0,1),
                                 &immersedpart);

    setStringToIntegralParameter<int>("COUPVARIABLE_PROTRUSION","Displacement",
                                 "Coupling variable at the protruding interface",
                                 tuple<std::string>("Displacement","Force"),
                                 tuple<int>(0,1),
                                 &immersedpart);

    setStringToIntegralParameter<int>("COUPMETHOD","immersed",
                                 "Coupling Method Mortar (mtr) or conforming nodes at interface",
                                 tuple<std::string>(
                                   "MTR",
                                   "Mtr",
                                   "mtr",
                                   "conforming",
                                   "immersed"
                                   ),
                                 tuple<int>(0,0,0,1,2),
                                 &immersedpart);

    DoubleParameter("BASETOL",1e-3,
                    "Basic tolerance for adaptive convergence check in monolithic FSI.\n"
                    "This tolerance will be used for the linear solve of the FSI block system.\n"
                    "The linear convergence test will always use the relative residual norm (AZ_r0).\n"
                    "Not to be confused with the Newton tolerance (CONVTOL) that applies\n"
                    "to the nonlinear convergence test using a absolute residual norm.",
                    &immersedpart);

    DoubleParameter("CONVTOL",1e-6,"Tolerance for iteration over fields in case of partitioned scheme",&immersedpart);
    DoubleParameter("RELAX",1.0,"fixed relaxation parameter for partitioned FSI solvers",&immersedpart);
    DoubleParameter("MAXOMEGA",0.0,"largest omega allowed for Aitken relaxation (0.0 means no constraint)",&immersedpart);
    IntParameter("ITEMAX",100,"Maximum number of iterations over fields",&immersedpart);

}



void INPAR::IMMERSED::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
    // IMMERSED FSI

    Teuchos::RCP<ConditionDefinition> immersedsearchbox =
      Teuchos::rcp(new ConditionDefinition("DESIGN VOLUME IMMERSED SEARCHBOX",
                                           "ImmersedSearchbox",
                                           "Immersed Searchbox",
                                           DRT::Condition::ImmersedSearchbox,
                                           true,
                                           DRT::Condition::Volume));

    condlist.push_back(immersedsearchbox);

    /*--------------------------------------------------------------------*/
      // IMMERSED COUPLING

    std::vector<Teuchos::RCP<ConditionComponent> > immersedcomponents;

    immersedcomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

    Teuchos::RCP<ConditionDefinition> lineimmersed =
      Teuchos::rcp(new ConditionDefinition("DESIGN IMMERSED COUPLING LINE CONDITIONS",
                                           "IMMERSEDCoupling",
                                           "IMMERSED Coupling",
                                           DRT::Condition::IMMERSEDCoupling,
                                           true,
                                           DRT::Condition::Line));
    Teuchos::RCP<ConditionDefinition> surfimmersed =
      Teuchos::rcp(new ConditionDefinition("DESIGN IMMERSED COUPLING SURF CONDITIONS",
                                           "IMMERSEDCoupling",
                                           "IMMERSED Coupling",
                                           DRT::Condition::IMMERSEDCoupling,
                                           true,
                                           DRT::Condition::Surface));

    for (unsigned i=0; i<immersedcomponents.size(); ++i)
    {
      lineimmersed->AddComponent(immersedcomponents[i]);
      surfimmersed->AddComponent(immersedcomponents[i]);
    }

    condlist.push_back(lineimmersed);
    condlist.push_back(surfimmersed);


}

