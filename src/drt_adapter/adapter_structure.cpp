/*----------------------------------------------------------------------*/
/*!
\file adapter_structure.cpp

\brief Structure field adapter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_structure.H"
#include "adapter_structure_strugenalpha.H"
#include "adapter_structure_timint.H"
#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::Structure::~Structure()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureBaseAlgorithm::StructureBaseAlgorithm(const Teuchos::ParameterList& prbdyn)
{
  SetupStructure(prbdyn);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureBaseAlgorithm::~StructureBaseAlgorithm()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithm::SetupStructure(const Teuchos::ParameterList& prbdyn)
{
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // major switch to different time integrators
  switch (Teuchos::getIntegralValue<int>(sdyn,"DYNAMICTYP"))
  {
  case STRUCT_DYNAMIC::centr_diff :
    dserror("no central differences in DRT");
    break;
  case STRUCT_DYNAMIC::gen_alfa :
  case STRUCT_DYNAMIC::gen_alfa_statics :
    SetupStruGenAlpha(prbdyn);  // <-- here is the show
    break;
  case STRUCT_DYNAMIC::Gen_EMM :
    dserror("Gen_EMM not supported");
    break;
  case STRUCT_DYNAMIC::genalpha :
  case STRUCT_DYNAMIC::onesteptheta :
  case STRUCT_DYNAMIC::gemm :
    SetupTimIntImpl(prbdyn);  // <-- here is the show
    break;
  case STRUCT_DYNAMIC::ab2 :
    dserror("explicitly no");
    break;
  default :
    dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
    break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithm::SetupStruGenAlpha(const Teuchos::ParameterList& prbdyn)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("ADAPTER::StructureBaseAlgorithm::SetupStructure");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf,0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter> output =
    rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR*         actsolv  = &solv[genprob.numsf];

  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& sdyn     = DRT::Problem::Instance()->StructuralDynamicParams();

  //const Teuchos::ParameterList& size     = DRT::Problem::Instance()->ProblemSizeParams();

  if ((actdis->Comm()).MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, sdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<ParameterList> solveparams = rcp(new ParameterList());
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // create a generalized alpha time integrator
  // -------------------------------------------------------------------
  RCP<ParameterList> genalphaparams = rcp(new ParameterList());
  StruGenAlpha::SetDefaults(*genalphaparams);

  genalphaparams->set<string>("DYNAMICTYP",sdyn.get<string>("DYNAMICTYP"));
  genalphaparams->set<bool>  ("damping",Teuchos::getIntegralValue<int>(sdyn,"DAMPING"));
  genalphaparams->set<double>("damping factor K",sdyn.get<double>("K_DAMP"));
  genalphaparams->set<double>("damping factor M",sdyn.get<double>("M_DAMP"));

  genalphaparams->set<double>("beta",sdyn.get<double>("BETA"));
  genalphaparams->set<double>("gamma",sdyn.get<double>("GAMMA"));
  genalphaparams->set<double>("alpha m",sdyn.get<double>("ALPHA_M"));
  genalphaparams->set<double>("alpha f",sdyn.get<double>("ALPHA_F"));

  genalphaparams->set<double>("total time",0.0);
  genalphaparams->set<double>("delta time",prbdyn.get<double>("TIMESTEP"));
  genalphaparams->set<int>   ("step",0);
  genalphaparams->set<int>   ("nstep",prbdyn.get<int>("NUMSTEP"));
  genalphaparams->set<int>   ("max iterations",sdyn.get<int>("MAXITER"));
  genalphaparams->set<int>   ("num iterations",-1);
  genalphaparams->set<double>("tolerance displacements",sdyn.get<double>("TOLDISP"));
  genalphaparams->set<string>("convcheck",sdyn.get<string>("CONV_CHECK"));

  genalphaparams->set<bool>  ("io structural disp",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_DISP"));
  genalphaparams->set<int>   ("io disp every nstep",prbdyn.get<int>("UPRES"));

  switch (Teuchos::getIntegralValue<STRUCT_STRESS_TYP>(ioflags,"STRUCT_STRESS"))
  {
  case struct_stress_none:
    genalphaparams->set<string>("io structural stress", "none");
    break;
  case struct_stress_cauchy:
    genalphaparams->set<string>("io structural stress", "cauchy");
    break;
  case struct_stress_pk:
    genalphaparams->set<string>("io structural stress", "2PK");
    break;
  default:
    genalphaparams->set<string>("io structural stress", "none");
    break;
  }

  genalphaparams->set<int>   ("io stress every nstep",sdyn.get<int>("RESEVRYSTRS"));

  genalphaparams->set<int>   ("restart",probtype.get<int>("RESTART"));
  genalphaparams->set<int>   ("write restart every",prbdyn.get<int>("RESTARTEVRY"));

  genalphaparams->set<bool>  ("print to screen",true);
  genalphaparams->set<bool>  ("print to err",true);
  genalphaparams->set<FILE*> ("err file",allfiles.out_err);

  switch (Teuchos::getIntegralValue<int>(sdyn,"NLNSOL"))
  {
  case STRUCT_DYNAMIC::fullnewton:
    genalphaparams->set<string>("equilibrium iteration","full newton");
    break;
  case STRUCT_DYNAMIC::lsnewton:
    genalphaparams->set<string>("equilibrium iteration","line search newton");
    break;
  case STRUCT_DYNAMIC::modnewton:
    genalphaparams->set<string>("equilibrium iteration","modified newton");
    break;
  case STRUCT_DYNAMIC::nlncg:
    genalphaparams->set<string>("equilibrium iteration","nonlinear cg");
    break;
  case STRUCT_DYNAMIC::ptc:
    genalphaparams->set<string>("equilibrium iteration","ptc");
    break;
  default:
    genalphaparams->set<string>("equilibrium iteration","full newton");
    break;
  }

  switch (Teuchos::getIntegralValue<STRUCT_STRAIN_TYP>(ioflags,"STRUCT_STRAIN"))
  {
  case struct_strain_none:
    genalphaparams->set<string>("io structural strain", "none");
  break;
  case struct_strain_ea:
    genalphaparams->set<string>("io structural strain", "euler_almansi");
  break;
  case struct_strain_gl:
    genalphaparams->set<string>("io structural strain", "green_lagrange");
  break;
  default:
    genalphaparams->set<string>("io structural strain", "none");
  break;
  }


  // set predictor (takes values "constant" "consistent")
  switch (Teuchos::getIntegralValue<int>(sdyn,"PREDICT"))
  {
  case STRUCT_DYNAMIC::pred_vague:
    dserror("You have to define the predictor");
    break;
  case STRUCT_DYNAMIC::pred_constdis:
    genalphaparams->set<string>("predictor","consistent");
    break;
  case STRUCT_DYNAMIC::pred_constdisvelacc:
    genalphaparams->set<string>("predictor","constant");
    break;
  default:
    dserror("Cannot cope with choice of predictor");
    break;
  }

  // sanity checks and default flags
  if (genprob.probtyp == prb_fsi)
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

    // robin flags
    INPUTPARAMS::FSIPartitionedCouplingMethod method =
      Teuchos::getIntegralValue<INPUTPARAMS::FSIPartitionedCouplingMethod>(fsidyn,"PARTITIONED");
    genalphaparams->set<bool>  ("structrobin",
                                method==INPUTPARAMS::fsi_DirichletRobin or method==INPUTPARAMS::fsi_RobinRobin);

    genalphaparams->set<double>("alpha s",fsidyn.get<double>("ALPHA_S"));

    int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling == fsi_iter_monolithic or
        coupling == fsi_iter_monolithiclagrange or
        coupling == fsi_iter_monolithicstructuresplit)
    {
      if (Teuchos::getIntegralValue<int>(sdyn,"PREDICT")!=STRUCT_DYNAMIC::pred_constdisvelacc)
        dserror("only constant structure predictor with monolithic FSI possible");

#if 0
      // overwrite time integration flags
      genalphaparams->set<double>("gamma",fsidyn.get<double>("GAMMA"));
      genalphaparams->set<double>("alpha m",fsidyn.get<double>("ALPHA_M"));
      genalphaparams->set<double>("alpha f",fsidyn.get<double>("ALPHA_F"));
#endif
    }
  }

  structure_ = rcp(new StructureGenAlpha(genalphaparams,actdis,solver,output));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureBaseAlgorithm::SetupTimIntImpl(const Teuchos::ParameterList& prbdyn)
{
  // this is not exactly a one hundred meter race, but we need timing
  Teuchos::RCP<Teuchos::Time> t 
    = Teuchos::TimeMonitor::getNewTimer("ADAPTER::StructureTimIntBaseAlgorithm::SetupStructure");
  Teuchos::TimeMonitor monitor(*t);

  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf, 0);

  // set degrees of freedom in the discretization
  if (not actdis->Filled()) actdis->FillComplete();

  // context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output
    = Teuchos::rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0, 0.0);

  // set some pointers and variables
  SOLVAR* actsolv = &solv[genprob.numsf];

  // get input parameter lists and copy them, because a few parameters are overwritten
  //const Teuchos::ParameterList& probtype 
  //  = DRT::Problem::Instance()->ProblemTypeParams();
  Teuchos::RCP<Teuchos::ParameterList> ioflags 
    = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->IOParams()));
  Teuchos::RCP<Teuchos::ParameterList> sdyn 
    = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->StructuralDynamicParams()));
  //const Teuchos::ParameterList& size
  //  = DRT::Problem::Instance()->ProblemSizeParams();

  // show default parameters
  if ((actdis->Comm()).MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, *sdyn);

  // add extra parameters (a kind of work-around)
  Teuchos::RCP<Teuchos::ParameterList> xparams 
    = Teuchos::rcp(new Teuchos::ParameterList());
  xparams->set<FILE*>("err file", allfiles.out_err);

  // overrule certain parameters
  sdyn->set<double>("TIMESTEP", prbdyn.get<double>("TIMESTEP"));
  sdyn->set<int>("NUMSTEP", prbdyn.get<int>("NUMSTEP"));
  sdyn->set<int>("RESEVRYDISP", prbdyn.get<int>("UPRES"));
  sdyn->set<int>("RESTARTEVRY", prbdyn.get<int>("RESTARTEVRY"));

  // sanity checks and default flags
  if (genprob.probtyp == prb_fsi)
  {
    // FSI input parameters
    const Teuchos::ParameterList& fsidyn 
      = DRT::Problem::Instance()->FSIDynamicParams();

    // Robin flags
    INPUTPARAMS::FSIPartitionedCouplingMethod method
      = Teuchos::getIntegralValue<INPUTPARAMS::FSIPartitionedCouplingMethod>(fsidyn,"PARTITIONED");
    xparams->set<bool>("structrobin",
                       ( (method==INPUTPARAMS::fsi_DirichletRobin) 
                         or (method==INPUTPARAMS::fsi_RobinRobin) ));

    // THIS SHOULD GO, OR SHOULDN'T IT?
    xparams->set<double>("alpha s", fsidyn.get<double>("ALPHA_S"));

    // check if predictor fits to FSI algo
    int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
    if ( (coupling == fsi_iter_monolithic)
         or (coupling == fsi_iter_monolithiclagrange)
         or (coupling == fsi_iter_monolithicstructuresplit) )
    {
      if (Teuchos::getIntegralValue<int>(*sdyn,"PREDICT")
          != STRUCT_DYNAMIC::pred_constdisvelacc)
      {
        dserror("only constant structure predictor with monolithic FSI possible");
      }
    }
  }

  // create a solver
  Teuchos::RCP<ParameterList> solveparams 
    = Teuchos::rcp(new ParameterList());
  Teuchos::RCP<LINALG::Solver> solver
    = Teuchos::rcp(new LINALG::Solver(solveparams,
                                      actdis->Comm(),
                                      allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams, actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // create marching time integrator
  structure_ = Teuchos::rcp(new StructureTimInt(ioflags, sdyn, xparams,
                                                actdis, solver, output));

  // see you
  return;
}

/*----------------------------------------------------------------------*/
#endif
