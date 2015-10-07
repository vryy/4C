/*----------------------------------------------------------------------*/
/*!
\file inpar_acou.cpp

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_acou.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::ACOU::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& acousticdyn = list->sublist("ACOUSTIC DYNAMIC",false,"control parameters for acoustic or photoacoustic problems\n");

  DoubleParameter("TIMESTEP",0.01,"Time increment dt",&acousticdyn);
  IntParameter("NUMSTEP",100,"Total number of time steps",&acousticdyn);
  DoubleParameter("MAXTIME",1.0,"Total simulation time",&acousticdyn);
  IntParameter("CALCERRORFUNCNO",-1,"Function for Error Calculation",&acousticdyn);

  IntParameter("UPRES",1,"Increment for writing solution",&acousticdyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&acousticdyn);
  IntParameter("LINEAR_SOLVER",-1,"Number of linear solver used for acoustical problem",&acousticdyn);
  IntParameter("STARTFUNCNO",-1,"Function for Initial Starting Field",&acousticdyn);
  IntParameter("SOURCETERMFUNCNO",-1,"Function for source term in volume",&acousticdyn);

  BoolParameter("ALLELESEQUAL","No","Yes, if all elements have same shape and material",&acousticdyn);

  // distinguish viscous and lossless flows
  setStringToIntegralParameter<int>("PHYSICAL_TYPE","lossless",
                    "fluid properties",
                    tuple<std::string>(
                    "lossless",
                    "solid"),
                    tuple<int>(
                    acou_lossless,
                    acou_solid),
                    &acousticdyn);
  // for viscous flows, one can specify if the displacement gradient or the stresses are outputted
  BoolParameter("WRITESTRESS","Yes","Output of stresses instead of displacement gradient",&acousticdyn);

  // photoacoustics
  BoolParameter("PHOTOACOU","No","Coupling with Scatra for Diffusive Light Transport",&acousticdyn);
  BoolParameter("MESHCONFORM","No","Conformity of scatra and acoustical mesh",&acousticdyn);

  // local postprocessing and p-adaptivity
  BoolParameter("ERRORMAPS","No","Output of error maps obtained by local postprocessing",&acousticdyn);
  BoolParameter("P_ADAPTIVITY","No","p-adaptivity in time integration",&acousticdyn);
  DoubleParameter("P_ADAPT_TOL",1.0e-15,"Error tolerance for p-adaptivity",&acousticdyn);

  // time integration
  setStringToIntegralParameter<int>("TIMEINT","impleuler",
                    "Type of time integration scheme",
                    tuple<std::string>(
                    "impleuler",
                    "expleuler",
                    "classrk4",
                    "lsrk45reg2",
                    "lsrk33reg2",
                    "lsrk45reg3",
                    "ssprk"),
                    tuple<int>(
                    acou_impleuler,
                    acou_expleuler,
                    acou_classrk4,
                    acou_lsrk45reg2,
                    acou_lsrk33reg2,
                    acou_lsrk45reg3,
                    acou_ssprk),
                    &acousticdyn);

  BoolParameter("WRITEMONITOR","No","Write a monitor file for Pressure Monitor Condition",&acousticdyn);
  setStringToIntegralParameter<int>("INV_ANALYSIS","none",
                 "Types of inverse analysis and on/off switch",
                 tuple<std::string>(
                   "none",
                   "pat"), // here, backprojection could be added
                 tuple<int>(
                   inv_none,
                   inv_pat),
                 &acousticdyn);


  Teuchos::ParameterList& acou_inv = acousticdyn.sublist("PA IMAGE RECONSTRUCTION",false,"");

  setStringToIntegralParameter<int>("OPTIMIZATION","LBFGS",
                                    "types of optimization algorithm",
                                    tuple<std::string>(
                                      "GradientDescent",
                                      "LBFGS"),
                                    tuple<int>(
                                      inv_gd,
                                      inv_lbfgs),
                                    &acou_inv);

  StringParameter("MONITORFILE","none.monitor","Filename of file containing measured pressure values",&acou_inv);
  BoolParameter("FDCHECK","No","Finite difference check",&acou_inv);
  DoubleParameter("INV_TOL",1e-16,"Tolerance for objective function of inverse pat analysis",&acou_inv);
  BoolParameter("INV_TOL_GRAD_YN","No","Flag to indicate check of the norm of the gradient",&acou_inv);
  DoubleParameter("INV_TOL_GRAD",0.0,"Tolerance for norm of gradient of inverse pat analysis",&acou_inv);
  BoolParameter("ELE_SCALING","No","Should gradient be scaled by element size?",&acou_inv);
  IntParameter("INV_MAX_RUN",10,"Maximal run number for inverse pat analysis",&acou_inv);
  IntParameter("INV_LS_MAX_RUN",10,"Maximal run number for line search in inverse pat analysis",&acou_inv);
  DoubleParameter("LS_DECREASECOND",0.0,"coefficient for calculation of sufficient decrease condition",&acou_inv);
  DoubleParameter("LS_STEPLENGTHRED",0.5,"step length is multiplied by this value if line search not yet sufficient",&acou_inv);
  DoubleParameter("EPSILON",-1.0,"tolerated distance in which measured curve=nod curve",&acou_inv);
  BoolParameter("BACKPROJECTION","No","PreEstimation of material parameters with back projection",&acou_inv);
  StringParameter("ACOUPARAMLIST","none","list of std::string of acoustical parameters to be optimized",&acou_inv);
  BoolParameter("SUCCESSIVEUPDATE","No","Optical and acoustical parameters update simultaneous or successive",&acou_inv);
  BoolParameter("SEGMENTATION","No","Update of material parameters in dependence, material identification",&acou_inv);
  StringParameter("SEGMENTATIONMATS","none.material","Filename of file containing table of materials",&acou_inv);
  DoubleParameter("ADAPTREGUPARAMS",0.0,"Adapt regularization paramaters for each iteration such that regularization contribution takes this relative value",&acou_inv);
  StringParameter("TIKH_ACOU_PARAMS_BORDERS","0.0 0.0 0.0 0.0","lower and upper allowed value for rho and c (rho_min rho_max c_min c_max)",&acou_inv);
  StringParameter("TIKH_OPTI_PARAMS_BORDERS","0.0 0.0 0.0 0.0","lower and upper allowed value for reac and D (reac_min reac_max D_min D_max)",&acou_inv);
  DoubleParameter("EQUALITYPENALTY",1.0,"penalty coefficient for the equality constraint for segmentation",&acou_inv);

}


/// set specific acoustic conditions
void INPAR::ACOU::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // absorbing boundary condition for acoustic problems
  // line
  Teuchos::RCP<ConditionDefinition> absorbing_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN ABSORBING LINE CONDITIONS",
                                         "Absorbing",
                                         "Absorbing line for acoustics",
                                         DRT::Condition::Absorb,
                                         true,
                                         DRT::Condition::Line));
  condlist.push_back(absorbing_line);

  // surface
  Teuchos::RCP<ConditionDefinition> absorbing_surface =
    Teuchos::rcp(new ConditionDefinition("DESIGN ABSORBING SURF CONDITIONS",
                                         "Absorbing",
                                         "Absorbing surface for acoustics",
                                         DRT::Condition::Absorb,
                                         true,
                                         DRT::Condition::Surface));
  condlist.push_back(absorbing_surface);

  /*--------------------------------------------------------------------*/
  // monitor condition for acoustic problems
  // line
  Teuchos::RCP<ConditionDefinition> pressmon_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN PRESSURE MONITOR LINE CONDITIONS",
                                         "PressureMonitor",
                                         "Pressure monitor line for acoustics",
                                         DRT::Condition::PressureMonitor,
                                         true,
                                         DRT::Condition::Line));
  condlist.push_back(pressmon_line);

  // surface
  Teuchos::RCP<ConditionDefinition> pressmon_surface =
    Teuchos::rcp(new ConditionDefinition("DESIGN PRESSURE MONITOR SURF CONDITIONS",
                                         "PressureMonitor",
                                         "Pressure monitor surface for acoustics",
                                         DRT::Condition::PressureMonitor,
                                         true,
                                         DRT::Condition::Surface));
  condlist.push_back(pressmon_surface);

}

