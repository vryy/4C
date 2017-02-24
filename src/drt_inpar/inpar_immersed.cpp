/*----------------------------------------------------------------------*/
/*!
\file inpar_immersed.cpp

\brief Input parameters for immersed

\level 1

<pre>
\maintainer Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/


#include "inpar_immersed.H"

#include "inpar_fsi.H"
#include "inpar_structure.H"

#include "drt_validparameters.H"
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

  Teuchos::Tuple<std::string,3> coupname;
  Teuchos::Tuple<int,3> couplabel;

  coupname[ 0] = "basic_sequ_stagg";                              couplabel[ 0] = cell_basic_sequ_stagg;
  coupname[ 1] = "iter_stagg_fixed_rel_param";                    couplabel[ 1] = cell_iter_stagg_fixed_rel_param;
  coupname[ 2] = "iter_stagg_AITKEN_rel_param";                   couplabel[ 2] = cell_iter_stagg_AITKEN_rel_param;



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
                               "DEFORM_BACKGROUND_MESH","no",
                               "switch between immersed with fixed or deformable background mesh",
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

  DoubleParameter("FLD_SRCHRADIUS_FAC",1.0,"fac times fluid ele. diag. length",&immersedmethod);
  DoubleParameter("STRCT_SRCHRADIUS_FAC",0.5,"fac times structure bounding box diagonal",&immersedmethod);
  IntParameter("NUM_GP_FLUID_BOUND",8,"number of gp in fluid elements cut by surface of immersed structure (higher number yields better mass conservation)",&immersedmethod);

  /*----------------------------------------------------------------------*/
  /* parameters for paritioned immersed solvers */
  Teuchos::ParameterList& immersedpart = immersedmethod.sublist("PARTITIONED SOLVER",false,"");

  setStringToIntegralParameter<int>("COUPALGO","iter_stagg_fixed_rel_param",
                                    "Iteration Scheme over the fields",
                                    coupname,
                                    couplabel,
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

    setStringToIntegralParameter<int>("COUPVARIABLE","Displacement",
                                 "Coupling variable at the fsi interface",
                                 tuple<std::string>("Displacement","Force"),
                                 tuple<int>(0,1),
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

