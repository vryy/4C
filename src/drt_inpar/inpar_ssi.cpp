/*----------------------------------------------------------------------*/
/*!
\file inpar_ssi.cpp

<pre>
Maintainer: Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_ssi.H"
#include "inpar_scatra.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::SSI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& ssidyn = list->sublist(
   "SSI CONTROL",false,
   "Control paramters for scatra structure interaction"
   );

  // Output type
  DoubleParameter("RESTARTEVRYTIME",0,"write restart possibility every RESTARTEVRY steps",&ssidyn);
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&ssidyn);
  // Time loop control
  IntParameter("NUMSTEP",200,"maximum number of Timesteps",&ssidyn);
  DoubleParameter("MAXTIME",1000.0,"total simulation time",&ssidyn);
  DoubleParameter("TIMESTEP",-1,"time step size dt",&ssidyn);
  BoolParameter("DIFFTIMESTEPSIZE","No","use different step size for scatra and solid",&ssidyn);
  DoubleParameter("UPRESTIME",0,"increment for writing solution",&ssidyn);
  IntParameter("UPRES",1,"increment for writing solution",&ssidyn);
  IntParameter("ITEMAX",10,"maximum number of iterations over fields",&ssidyn);
  BoolParameter("SCATRA_FROM_RESTART_FILE","No","read scatra result from restart files (use option 'restartfromfile' during execution of baci)",&ssidyn);
  StringParameter("SCATRA_FILENAME","nil","Control-file name for reading scatra results in SSI",&ssidyn);

  // Type of coupling strategy between the two fields
  setStringToIntegralParameter<int>(
                              "FIELDCOUPLING","matching",
                              "Type of coupling strategy between fields",
                              tuple<std::string>(
                                "matching",
                                "meshtying",
                                "volmortar"
                                ),
                              tuple<int>(
                                  coupling_match,
                                  coupling_meshtying,
                                  coupling_volmortar
                                ),
                              &ssidyn);

  // Coupling strategy for SSI solvers
  setStringToIntegralParameter<int>(
                              "COUPALGO","ssi_IterStagg",
                              "Coupling strategies for SSI solvers",
                              tuple<std::string>(
                                "ssi_OneWay_ScatraToSolid",
                                "ssi_OneWay_SolidToScatra",
//                                "ssi_SequStagg_ScatraToSolid",
//                                "ssi_SequStagg_SolidToScatra",
                                "ssi_IterStagg",
                                "ssi_IterStaggFixedRel_ScatraToSolid",
                                "ssi_IterStaggFixedRel_SolidToScatra",
                                "ssi_IterStaggAitken_ScatraToSolid",
                                "ssi_IterStaggAitken_SolidToScatra"
                                ),
                              tuple<int>(
                                ssi_OneWay_ScatraToSolid,
                                ssi_OneWay_SolidToScatra,
//                                ssi_SequStagg_ScatraToSolid,
//                                ssi_SequStagg_SolidToScatra,
                                ssi_IterStagg,
                                ssi_IterStaggFixedRel_ScatraToSolid,
                                ssi_IterStaggFixedRel_SolidToScatra,
                                ssi_IterStaggAitken_ScatraToSolid,
                                ssi_IterStaggAitken_SolidToScatra
                                ),
                              &ssidyn);

  // type of scalar transport
  setStringToIntegralParameter<int>("SCATRATYPE","Undefined",
                               "Type of scalar transport problem",
                               tuple<std::string>(
                                 "Undefined",
                                 "ConvectionDiffusion",
                                 "Advanced_Reaction",
                                 "Cardiac_Monodomain",
                                 "Chemotaxis",
                                 "Chemo_Reac"),
                               tuple<int>(
                                 INPAR::SCATRA::impltype_undefined,
                                 INPAR::SCATRA::impltype_std,
                                 INPAR::SCATRA::impltype_advreac,
                                 INPAR::SCATRA::impltype_cardiac_monodomain,
                                 INPAR::SCATRA::impltype_chemo,
                                 INPAR::SCATRA::impltype_chemoreac),
                                 &ssidyn);

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned SSI */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ssidynpart = ssidyn.sublist(
      "PARTITIONED",false,
      "Partitioned Structure Scalar Interaction\n"
      "Control section for partitioned SSI"
       );

  // Solver parameter for relaxation of iterative staggered partitioned SSI
  DoubleParameter("MAXOMEGA",10.0,"largest omega allowed for Aitken relaxation",&ssidynpart);
  DoubleParameter("MINOMEGA",0.1,"smallest omega allowed for Aitken relaxation",&ssidynpart);
  DoubleParameter("STARTOMEGA",1.0,"fixed relaxation parameter",&ssidynpart);

  // convergence tolerance of outer iteration loop
  DoubleParameter("CONVTOL",1e-6,"tolerance for convergence check of outer iteration within partitioned SSI",&ssidynpart);
}



void INPAR::SSI::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  Teuchos::RCP<ConditionDefinition> linessi =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING LINE CONDITIONS",
          "SSICoupling",
          "SSI Coupling",
          DRT::Condition::SSICoupling,
          true,
          DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfssi =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SURF CONDITIONS",
          "SSICoupling",
          "SSI Coupling",
          DRT::Condition::SSICoupling,
          true,
          DRT::Condition::Surface));

  condlist.push_back(linessi);
  condlist.push_back(surfssi);

}

