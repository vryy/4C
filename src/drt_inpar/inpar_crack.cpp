/*----------------------------------------------------------------------*/
/*!
\file inpar_crack.cpp

\brief Input parameters for crack propagation

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_crack.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::CRACK::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& crackdyn = list->sublist("CRACK",false,"");

  // type of crack propagation model -- either using linear elastic fracture mechanics concepts, or cohesive crack models
  setStringToIntegralParameter<int>("CRACK_MODEL","none",
                                      "type of crack propagation modeling",
                                      tuple<std::string>(
                                        "none",
                                        "lefm",
                                        "cohesive"),
                                      tuple<int>(
                                        crack_none,
                                        crack_lefm,
                                        crack_cohesive),
                                      &crackdyn);

  DoubleParameter("CRITICAL_K1",50000.0,"Critical stress intensity factor in normal mode",&crackdyn);

  DoubleParameter("CRITICAL_K2",50000.0,"Critical stress intensity factor in shear mode",&crackdyn);

  StringParameter("THICKNESS_ASSUMPTION","plane_strain",
                  "Whether is this plane strain or plane stress problem?",
                   &crackdyn);

  // type of cohesive crack propagation modeling
  setStringToIntegralParameter<int>("COHESIVE_CRACK_MODEL","none",
                                  "type of cohesive crack propagation modeling",
                                  tuple<std::string>(
                                    "none",
                                    "dczm",
                                    "ddzm"),
                                  tuple<int>(
                                    cohesive_none,
                                    cohesive_dczm,
                                    cohesive_ddzm),
                                  &crackdyn);

  setStringToIntegralParameter<int>("TRACTION_SEPARATION_LAW","exponential",
                                    "type of traction-separation law used for cohesive elements",
                                    tuple<std::string>(
                                      "linear",
                                      "trapezoidal",
                                      "exponential",
                                      "sinusoidal",
                                      "ppr"),
                                    tuple<int>(
                                      linear,
                                      trapezoidal,
                                      exponential,
                                      sinusoidal,
                                      ppr),
                                    &crackdyn);

  // are we modeling cracks with known propagation direction?
  setStringToIntegralParameter<int>("IS_CRACK_PREDEFINED","No","Have you already predefined the crack path?",
                               yesnotuple,yesnovalue,&crackdyn);

  DoubleParameter("NORMAL_COHESIVE_STRENGTH",50000.0,"Cohesive strength for normal separation",&crackdyn);
  DoubleParameter("SHEAR_COHESIVE_STRENGTH",500000000.0,"Cohesive strength for shear separation",&crackdyn);
  DoubleParameter("G_I",1.0,"Model I fracture energy (normal separation)",&crackdyn);
  DoubleParameter("G_II",1.0,"Model II fracture energy (shear separation)",&crackdyn);

  DoubleParameter("ALFA_PPR",3.0,"Constant alpha in PPR model",&crackdyn);
  DoubleParameter("BETA_PPR",3.0,"Constant beta in PPR model",&crackdyn);

  DoubleParameter("SLOPE_NORMAL",0.02,"Initial slope indicator in normal direction for PPR model",&crackdyn);
  DoubleParameter("SLOPE_SHEAR",0.02,"Initial slope indicator in normal direction for PPR model",&crackdyn);

  setStringToIntegralParameter<int>("GMSH_OUT","No","Do you want to write Gmsh output of displacement each timestep?",
                                 yesnotuple,yesnovalue,&crackdyn);

  IntParameter("START_NEW_NODE_ID",0,"Id of first node that will be introduced into discretization while propagating crack. This "
      "should be set greater than total no of nodes in the initial discretization",&crackdyn);

  IntParameter("START_NEW_ELE_ID",0,"Id of first wedge element that will be introduced into discretization while propagating crack. This "
        "should be set greater than total no of elements in the initial discretization",&crackdyn);

  // type of crack propagation model -- either using linear elastic fracture mechanics concepts, or cohesive crack models
  setStringToIntegralParameter<int>("CRACK_PROPAGATION_CRITERION","displacement_correlation",
                                      "Crack propagation criterion used for LEFM",
                                      tuple<std::string>(
                                        "displacement_correlation",
                                        "J_Integral"),
                                      tuple<int>(
                                        displacementCorrelation,
                                        J_Integral),
                                      &crackdyn);

  IntParameter("NO_LAYERS_J_INT",4,"No of element layers used for J-integral calculation",&crackdyn);

  DoubleParameter("CRITICAL_J",1.0,"Critical energy release rate",&crackdyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fsi_crackdyn = list->sublist("FSI CRACK",false,"");

  setStringToIntegralParameter<int>("CHECK_CONDITION","No","Do you want to check crack mouth opening condition?",
                                     yesnotuple,yesnovalue,&fsi_crackdyn);

  DoubleParameter("CRACK_OPENING_DIST",0.01,"Critial crack mouth opening distance",&fsi_crackdyn);

}



void INPAR::CRACK::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // crack

  std::vector<Teuchos::RCP<ConditionComponent> > crackdef;
  Teuchos::RCP<ConditionDefinition> crmaster =
      Teuchos::rcp(new ConditionDefinition("DESIGN CRACK MASTER SURFACE",
                                           "masterCrackSurface",
                                           "master crack surface",
                                           DRT::Condition::CrackMastersurface,
                                           true,
                                           DRT::Condition::Surface));

  Teuchos::RCP<ConditionDefinition> crslave =
        Teuchos::rcp(new ConditionDefinition("DESIGN CRACK SLAVE SURFACE",
                                             "slaveCrackSurface",
                                             "slave crack surface",
                                             DRT::Condition::CrackSlavesurface,
                                             true,
                                             DRT::Condition::Surface));

  Teuchos::RCP<ConditionDefinition> initialCrack =
          Teuchos::rcp(new ConditionDefinition("DESIGN CRACK INITIATION POINTS",
                                               "CrackInitiationPoints",
                                               "Points at which crack initiates",
                                               DRT::Condition::CrackInitPoints,
                                               true,
                                               DRT::Condition::Point));

  Teuchos::RCP<ConditionDefinition> crackBoundary =
          Teuchos::rcp(new ConditionDefinition("DESIGN CRACK BOUNDARY POINTS",
                                               "CrackBoundaryPoints",
                                               "All nodes on the boundary of the domain",
                                               DRT::Condition::CrackBoundaryPoints,
                                               true,
                                               DRT::Condition::Point));

  Teuchos::RCP<ConditionDefinition> crackInnerLayer =
        Teuchos::rcp(new ConditionDefinition("DESIGN CRACK INNER LAYER POINTS",
                                             "CrackInnerLayerPoints",
                                             "Points on the inner layer while performing crack initiation",
                                             DRT::Condition::CrackInnerLayerPoints,
                                             true,
                                             DRT::Condition::Point));

  Teuchos::RCP<ConditionDefinition> crackOuterLayer =
        Teuchos::rcp(new ConditionDefinition("DESIGN CRACK OUTER LAYER POINTS",
                                             "CrackOuterLayerPoints",
                                             "Points on the outer layer while performing crack initiation",
                                             DRT::Condition::CrackOuterLayerPoints,
                                             true,
                                             DRT::Condition::Point));

  for (unsigned i=0; i<crackdef.size(); ++i)
  {
    crmaster->AddComponent(crackdef[i]);
    crslave->AddComponent(crackdef[i]);
    initialCrack->AddComponent(crackdef[i]);
    crackBoundary->AddComponent(crackdef[i]);
    crackInnerLayer->AddComponent(crackdef[i]);
    crackOuterLayer->AddComponent(crackdef[i]);
  }

  condlist.push_back(crmaster);
  condlist.push_back(crslave);
  condlist.push_back(initialCrack);
  condlist.push_back(crackBoundary);
  condlist.push_back(crackInnerLayer);
  condlist.push_back(crackOuterLayer);

}
