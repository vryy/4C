/*----------------------------------------------------------------------*/
/*!
\file inpar_ssi.cpp

\brief input parameters for solid-scatra-interaction

<pre>
\level 2

\maintainer Moritz Thon
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
  DoubleParameter("RESULTSEVRYTIME",0,"increment for writing solution",&ssidyn);
  IntParameter("RESULTSEVRY",1,"increment for writing solution",&ssidyn);
  IntParameter("ITEMAX",10,"maximum number of iterations over fields",&ssidyn);
  BoolParameter("SCATRA_FROM_RESTART_FILE","No","read scatra result from restart files (use option 'restartfromfile' during execution of baci)",&ssidyn);
  StringParameter("SCATRA_FILENAME","nil","Control-file name for reading scatra results in SSI",&ssidyn);

  // Type of coupling strategy between the two fields
  setStringToIntegralParameter<int>(
                              "FIELDCOUPLING","volume_matching",
                              "Type of coupling strategy between fields",
                              tuple<std::string>(
                                "volume_matching",
                                "volume_nonmatching",
                                "boundary_nonmatching",
                                "volumeboundary_matching"
                                ),
                              tuple<int>(
                                  coupling_volume_match,
                                  coupling_volume_nonmatch,
                                  coupling_boundary_nonmatch,
                                  coupling_volumeboundary_match
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
                                 "RefConc_Reac",
                                 "Cardiac_Monodomain",
                                 "Chemotaxis",
                                 "Chemo_Reac",
                                 "Bond_Reac"),
                               tuple<int>(
                                 INPAR::SCATRA::impltype_undefined,
                                 INPAR::SCATRA::impltype_std,
                                 INPAR::SCATRA::impltype_advreac,
                                 INPAR::SCATRA::impltype_refconcreac,
                                 INPAR::SCATRA::impltype_cardiac_monodomain,
                                 INPAR::SCATRA::impltype_chemo,
                                 INPAR::SCATRA::impltype_chemoreac,
                                 INPAR::SCATRA::impltype_bondreac),
                                 &ssidyn);

  //Restart from Structure problem instead of SSI
  BoolParameter("RESTART_FROM_STRUCTURE","no","restart from structure problem (e.g. from prestress calculations) instead of ssi",&ssidyn);

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
  Teuchos::RCP<ConditionDefinition> linessiplain =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING LINE CONDITIONS",
          "SSICoupling",
          "SSI Coupling",
          DRT::Condition::SSICoupling,
          true,
          DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfssiplain =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SURF CONDITIONS",
          "SSICoupling",
          "SSI Coupling",
          DRT::Condition::SSICoupling,
          true,
          DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volssiplain =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING VOL CONDITIONS",
          "SSICoupling",
          "SSI Coupling",
          DRT::Condition::SSICoupling,
          true,
          DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<ConditionComponent> > ssicoupcomponentsplain;
  ssicoupcomponentsplain.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  // insert input file line components into condition definitions
  for (unsigned i=0; i<ssicoupcomponentsplain.size(); ++i)
  {
    linessiplain->AddComponent(ssicoupcomponentsplain[i]);
    surfssiplain->AddComponent(ssicoupcomponentsplain[i]);
    volssiplain ->AddComponent(ssicoupcomponentsplain[i]);
  }

  condlist.push_back(linessiplain);
  condlist.push_back(surfssiplain);
  condlist.push_back(volssiplain);

  /*--------------------------------------------------------------------*/
  //! set solid dofset on scatra discretization
  Teuchos::RCP<ConditionDefinition> linessi =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SOLIDTOSCATRA LINE CONDITIONS",
          "SSICouplingSolidToScatra",
          "SSI Coupling SolidToScatra",
          DRT::Condition::SSICouplingSolidToScatra,
          true,
          DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfssi =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SOLIDTOSCATRA SURF CONDITIONS",
          "SSICouplingSolidToScatra",
          "SSI Coupling SolidToScatra",
          DRT::Condition::SSICouplingSolidToScatra,
          true,
          DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volssi =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SOLIDTOSCATRA VOL CONDITIONS",
          "SSICouplingSolidToScatra",
          "SSI Coupling SolidToScatra",
          DRT::Condition::SSICouplingSolidToScatra,
          true,
          DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<ConditionComponent> > ssicoupcomponents;
  ssicoupcomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  // insert input file line components into condition definitions
  for (unsigned i=0; i<ssicoupcomponents.size(); ++i)
  {
    linessi->AddComponent(ssicoupcomponents[i]);
    surfssi->AddComponent(ssicoupcomponents[i]);
    volssi->AddComponent(ssicoupcomponents[i]);
  }

  condlist.push_back(linessi);
  condlist.push_back(surfssi);
  condlist.push_back(volssi);

  /*--------------------------------------------------------------------*/
  //! set scatra dofset on solid discretization
  Teuchos::RCP<ConditionDefinition> linessi2 =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SCATRATOSOLID LINE CONDITIONS",
          "SSICouplingScatraToSolid",
          "SSI Coupling ScatraToSolid",
          DRT::Condition::SSICouplingScatraToSolid,
          true,
          DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfssi2 =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SCATRATOSOLID SURF CONDITIONS",
          "SSICouplingScatraToSolid",
          "SSI Coupling ScatraToSolid",
          DRT::Condition::SSICouplingScatraToSolid,
          true,
          DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volssi2 =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SCATRATOSOLID VOL CONDITIONS",
          "SSICouplingScatraToSolid",
          "SSI Coupling ScatraToSolid",
          DRT::Condition::SSICouplingScatraToSolid,
          true,
          DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<ConditionComponent> > ssicoupcomponents2;
  ssicoupcomponents2.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  // insert input file line components into condition definitions
  for (unsigned i=0; i<ssicoupcomponents2.size(); ++i)
  {
    linessi2->AddComponent(ssicoupcomponents2[i]);
    surfssi2->AddComponent(ssicoupcomponents2[i]);
    volssi2 ->AddComponent(ssicoupcomponents2[i]);
  }

  condlist.push_back(linessi2);
  condlist.push_back(surfssi2);
  condlist.push_back(volssi2);

}

