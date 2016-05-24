/*----------------------------------------------------------------------*/
/*!
\file inpar_invanalysis.cpp

\brief Input parameters for inverse analysis

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@lnm.mw.tum.de
            089-289-10361
</pre>
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_invanalysis.H"
#include "inpar_statinvanalysis.H"
#include "inpar_structure.H"
#include "../drt_lib/drt_conditiondefinition.H"


void INPAR::INVANA::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& iap = list->sublist("INVERSE ANALYSIS",false,"");

  // Inverse Analysis
  setStringToIntegralParameter<int>("INV_ANALYSIS","none",
                               "types of inverse analysis and on/off switch",
                               tuple<std::string>(
                                 "none",
                                 "lung",
                                 "gen"),
                               tuple<int>(
                                 INPAR::STR::inv_none,
                                 INPAR::STR::inv_lung,
                                 INPAR::STR::inv_generalized),
                               &iap);

  DoubleParameter("MC_X_0",0.0,"measured displacment of the tension testing in x dir",&iap);
  DoubleParameter("MC_X_1",0.0,"measured displacment of the tension testing in x dir",&iap);
  DoubleParameter("MC_X_2",0.0,"measured displacment of the tension testing in x dir",&iap);
  DoubleParameter("MC_Y_0",0.0,"measured displacment of the tension testing in y dir",&iap);
  DoubleParameter("MC_Y_1",0.0,"measured displacment of the tension testing in y dir",&iap);
  DoubleParameter("MC_Y_2",0.0,"measured displacment of the tension testing in y dir",&iap);

  // tolerance for inv_analysis
  DoubleParameter("INV_ANA_TOL",1.0,"tolerance for inverse analysis",&iap);
  IntParameter("INV_ANA_MAX_RUN",100,"max iterations for inverse analysis",&iap);

  // perturbation parameters
  DoubleParameter("INV_ALPHA",1.0e-3,"perturbation parameters",&iap);
  DoubleParameter("INV_BETA",1.0e-3,"perturbation parameters",&iap);

  // initial regularization parameter
  DoubleParameter("INV_INITREG",1.0,"initial regularization parameter",&iap);

  // strategy to update regularization parameter
  setStringToIntegralParameter<int>("UPDATE_REG","RES","Update strategy for regularization parameter ",
                                    tuple<std::string>("RES","res",
                                                        "GRAD","grad"),
                                    tuple<int>(
                                        INPAR::STR::reg_update_res,INPAR::STR::reg_update_res,
                                        INPAR::STR::reg_update_grad,INPAR::STR::reg_update_grad
                                  ),
                                 &iap);


  StringParameter("MONITORFILE","none.monitor",
                  "filename of file containing measured displacements",
                  &iap);

  setNumericStringParameter("INV_LIST","-1",
                            "IDs of materials that have to be fitted",
                            &iap);

  setNumericStringParameter("INV_EH_LIST","-1",
                            "IDs of materials that have to be fitted",
                            &iap);

  setStringToIntegralParameter<int>("NEW_FILES","yes",
                                    "new result files for each run",
                                    yesnotuple,yesnovalue,&iap);
  setStringToIntegralParameter<int>("PARAM_BOUNDS","no",
                                      "Reset parameters if optstep predicts negative values",
                                      yesnotuple,yesnovalue,&iap);

  BoolParameter("PATCHES","No","Do you want to use smoothed patches?",&iap);
  StringParameter("DEFINEPATCHES","MaterialNumber",
      "define how the patches are defined: MaterialNumber or Uniform",
      &iap);
  IntParameter("NUMPATCHES",0,"number of patches",&iap);
  setNumericStringParameter("INV_LIST_PATCHES","-1",
                      "IDs of materials that are included in the patches",
                      &iap);
  IntParameter("SMOOTHINGSTEPSPATCHES",1,"number of smoothing steps that are performed",&iap);
  setNumericStringParameter("STARTVALUESFORPATCHES","1.0",
                  "startvalues for the patches, only needed for Uniform Patches",
                  &iap);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& statinvp = list->sublist("STAT INVERSE ANALYSIS",false,"");

  // Statistical Inverse Analysis switch
  setStringToIntegralParameter<int>("STAT_INV_ANALYSIS","none",
                                    "types of statistical inverse analysis and on/off switch",
                                    tuple<std::string>(
                                      "none",
                                      "GradientDescent",
                                      "MonteCarlo",
                                      "LBFGS"),
                                    tuple<int>(
                                      stat_inv_none,
                                      stat_inv_graddesc,
                                      stat_inv_mc,
                                      stat_inv_lbfgs),
                                    &statinvp);

  // initial scaling for the LBFGS algorithm
  BoolParameter("LBFGSINITSCAL","yes","want initial scaling for the LBFGS?", &statinvp);

  // step to restart from
  IntParameter("FPRESTART",0,"forward problem restart",&statinvp);

  // write restart info every so often
  IntParameter("RESTARTEVRY",1,"write restart information every x-th step",&statinvp);

  // decide which parametrization of material parameters to use
  setStringToIntegralParameter<int>("PARAMETRIZATION","none",
                                      "how to parametrize the parameter field",
                                    tuple<std::string>(
                                      "none",
                                      "kernelsmoothing",
                                      "elementwise",
                                      "uniform"),
                                    tuple<int>(
                                      stat_inv_mp_none,
                                      stat_inv_mp_smoothkernel,
                                      stat_inv_mp_elementwise,
                                      stat_inv_mp_uniform),
                                    &statinvp);

  // want some regularization
  setStringToIntegralParameter<int>("REGULARIZATION","none",
                                    "want regularization? ('tikhonov', 'totalvariation', 'none')",
                                    tuple<std::string>(
                                      "none",
                                      "tikhonov",
                                      "totalvariation"),
                                    tuple<int>(
                                      stat_inv_reg_none,
                                      stat_inv_reg_tikhonov,
                                      stat_inv_reg_totalvariation),
                                    &statinvp);

  // objective function
  setStringToIntegralParameter<int>("OBJECTIVEFUNCT","none",
                                    "choose type of objective function ('displacements', 'surfcurr')",
                                    tuple<std::string>(
                                      "none",
                                      "displacements",
                                      "surfcurr"),
                                    tuple<int>(
                                      stat_inv_obj_none,
                                      stat_inv_obj_disp,
                                      stat_inv_obj_surfcurr),
                                    &statinvp);

  // scaling of objective function
  BoolParameter("OBJECTIVEFUNCTSCAL","No","want scaling of objective function?", &statinvp);

  // monitorfile to provide measurements
  StringParameter("MONITORFILE","none.monitor",
                  "filename of file containing measured displacements",
                  &statinvp);

  // target discretization for surface currents
  StringParameter("TARGETDISCRETIZATION", "none.dat",
                  "datfile containing target discretization",
                  &statinvp);

  // list of parameters for the respective material
  StringParameter("PARAMLIST","none",
                  "list of std::string of parameters to be optimized, order as in INV_LIST e.g. 1 YOUNG BETA",
                  &statinvp);

  // number of optimization steps
  IntParameter("MAXITER",100,"max iterations for inverse analysis",&statinvp);

  // number of optimization steps before using
  // parameter continuation in the forward problem
  IntParameter("ITERTOPC",10,"iterations before parameter continuation in the forward problem",&statinvp);

  // give prestressing method to be used for the adjoint formulation
  setStringToIntegralParameter<int>("PRESTRESS","none","prestressing takes values none mulf id",
                               tuple<std::string>("none","None","NONE",
                                                  "mulf","Mulf","MULF",
                                                  "id","Id","ID"),
                               tuple<int>(INPAR::STR::prestress_none,INPAR::STR::prestress_none,INPAR::STR::prestress_none,
                                          INPAR::STR::prestress_mulf,INPAR::STR::prestress_mulf,INPAR::STR::prestress_mulf,
                                          INPAR::STR::prestress_id,INPAR::STR::prestress_id,INPAR::STR::prestress_id),
                               &statinvp);

  // stepsize for deterministic gradient based schemes
  DoubleParameter("STEPSIZE",1.0,"stepsize for the gradient descent scheme",&statinvp);

  // convergence criterion tolerance
  DoubleParameter("CONVTOL",1.0e-06,"stop optimizaiton iterations for convergence criterion below this value",&statinvp);

  // weight of the Tikhonov regularization
  DoubleParameter("REG_WEIGHT",0.1,"weight of the regularization functional",&statinvp);

  // mean value of the Tikhonov regularization
  DoubleParameter("REG_MEAN",0.0,"mean value for tikhonov regularization",&statinvp);

  // regularization of the totalvariation functional
  DoubleParameter("TVD_EPS",0.001,"differentiation epsilon for total variation",&statinvp);

  // number of optimization steps
  IntParameter("SIZESTORAGE",20,"number of vectors to keep in storage; defaults to 20 (lbfgs usage only)",&statinvp);

  // meta parametrization of material parameters
  setStringToIntegralParameter<int>("METAPARAMS","none",
                                    "choose type of metaparametrization (none/quad/arctan)",
                                    tuple<std::string>(
                                      "none",
                                      "quad",
                                      "arctan"),
                                    tuple<int>(
                                      stat_inv_meta_none,
                                      stat_inv_meta_quad,
                                      stat_inv_meta_arctan),
                                    &statinvp);


  // scale of the kernel functions used in surface currents
  DoubleParameter("KERNELSCALE", -1.0, "scale of the kernel function", &statinvp);

}



void INPAR::INVANA::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // inverse analysis fitted surface

  std::vector<Teuchos::RCP<ConditionComponent> > invanacomponents;
  invanacomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  Teuchos::RCP<ConditionDefinition> surfinvana =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE INV ANALYSIS",
                                         "SurfInvAna",
                                         "Inverse Analysis Surface",
                                         DRT::Condition::InvAnaSurface,
                                         true,
                                         DRT::Condition::Surface));

  surfinvana->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  condlist.push_back(surfinvana);

  /*--------------------------------------------------------------------*/
  // Surface current evaluation condition

  std::vector<Teuchos::RCP<ConditionComponent> > surfcurrcomponents;
  surfcurrcomponents.push_back(Teuchos::rcp(new IntConditionComponent("matching id")));

  Teuchos::RCP<ConditionDefinition> surfcurrcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE CURRENT EVALUATION CONDITION",
                                         "SurfaceCurrent",
                                         "Surface current",
                                         DRT::Condition::SurfaceCurrent,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<surfcurrcomponents.size(); ++i)
  {
    surfcurrcond->AddComponent(surfcurrcomponents[i]);
  }

  condlist.push_back(surfcurrcond);

  /*--------------------------------------------------------------------*/
  // Uncertain surface condition

  std::vector<Teuchos::RCP<ConditionComponent> > uncertsurfcomponents;
  uncertsurfcomponents.push_back(Teuchos::rcp(new IntConditionComponent("matching id")));

  Teuchos::RCP<ConditionDefinition> uncertsurfcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN UNCERTAIN SURFACE CONDITION",
                                         "UncertainSurface",
                                         "Uncertain surface",
                                         DRT::Condition::UncertainSurface,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<uncertsurfcomponents.size(); ++i)
  {
    uncertsurfcond->AddComponent(uncertsurfcomponents[i]);
  }

  condlist.push_back(uncertsurfcond);

}

