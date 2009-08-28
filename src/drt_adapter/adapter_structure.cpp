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
#include "adapter_structure_contactstrugenalpha.H"
#include "adapter_structure_timint.H"
#include "adapter_structure_constrained.H"
#include "adapter_structure_wrapper.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_contact.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

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
  switch (Teuchos::getIntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
  {
  case INPAR::STR::dyna_centr_diff :
    dserror("no central differences in DRT");
    break;
  case INPAR::STR::dyna_gen_alfa :
  case INPAR::STR::dyna_gen_alfa_statics :
    SetupStruGenAlpha(prbdyn);  // <-- here is the show
    break;
  case INPAR::STR::dyna_Gen_EMM :
    dserror("Gen_EMM not supported");
    break;
  case INPAR::STR::dyna_statics :
  case INPAR::STR::dyna_genalpha :
  case INPAR::STR::dyna_onesteptheta :
  case INPAR::STR::dyna_gemm :
    SetupTimIntImpl(prbdyn);  // <-- here is the show
    break;
  case INPAR::STR::dyna_ab2 :
//  case STRUCT_DYNAMIC::euma :
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
  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& sdyn     = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& scontact = DRT::Problem::Instance()->StructuralContactParams();

  //const Teuchos::ParameterList& size     = DRT::Problem::Instance()->ProblemSizeParams();

  if ((actdis->Comm()).MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, sdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(DRT::Problem::Instance()->StructSolverParams(),
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // -------------------------------------------------------------------
  // create a generalized alpha time integrator
  // -------------------------------------------------------------------
  RCP<ParameterList> genalphaparams = rcp(new ParameterList());
  StruGenAlpha::SetDefaults(*genalphaparams);

  genalphaparams->set<string>("DYNAMICTYP",sdyn.get<string>("DYNAMICTYP"));
  genalphaparams->set<bool>  ("damping",(not (sdyn.get<std::string>("DAMPING") == "no"
                                              or sdyn.get<std::string>("DAMPING") == "No"
                                              or sdyn.get<std::string>("DAMPING") == "NO")));
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

  INPAR::STR::StressType iostress = Teuchos::getIntegralValue<INPAR::STR::StressType>(ioflags,"STRUCT_STRESS");
  genalphaparams->set<INPAR::STR::StressType>("io structural stress", iostress);
  genalphaparams->set<int>   ("io stress every nstep",sdyn.get<int>("RESEVRYSTRS"));
  INPAR::STR::StrainType iostrain = Teuchos::getIntegralValue<INPAR::STR::StrainType>(ioflags,"STRUCT_STRAIN");
  genalphaparams->set<INPAR::STR::StrainType>("io structural strain", iostrain);

  genalphaparams->set<int>   ("restart",probtype.get<int>("RESTART"));
  genalphaparams->set<int>   ("write restart every",prbdyn.get<int>("RESTARTEVRY"));

  genalphaparams->set<bool>  ("print to screen",true);
  genalphaparams->set<bool>  ("print to err",true);
  genalphaparams->set<FILE*> ("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

  switch (Teuchos::getIntegralValue<INPAR::STR::NonlinSolTech>(sdyn,"NLNSOL"))
  {
  case INPAR::STR::soltech_newtonfull:
    genalphaparams->set<string>("equilibrium iteration","full newton");
    break;
  case INPAR::STR::soltech_newtonls:
    genalphaparams->set<string>("equilibrium iteration","line search newton");
    break;
  case INPAR::STR::soltech_newtonmod:
    genalphaparams->set<string>("equilibrium iteration","modified newton");
    break;
  case INPAR::STR::soltech_nlncg:
    genalphaparams->set<string>("equilibrium iteration","nonlinear cg");
    break;
  case INPAR::STR::soltech_ptc:
    genalphaparams->set<string>("equilibrium iteration","ptc");
    break;
  case INPAR::STR::soltech_newtonuzawalin:
    genalphaparams->set<string>("equilibrium iteration","newtonlinuzawa");
  break;
  case INPAR::STR::soltech_newtonuzawanonlin:
    genalphaparams->set<string>("equilibrium iteration","augmentedlagrange");
    break;
  default:
    genalphaparams->set<string>("equilibrium iteration","full newton");
    break;
  }

  // set predictor (takes values "constant" "consistent")
  switch (Teuchos::getIntegralValue<INPAR::STR::PredEnum>(sdyn,"PREDICT"))
  {
  case INPAR::STR::pred_vague:
    dserror("You have to define the predictor");
    break;
  case INPAR::STR::pred_constdis:
    genalphaparams->set<string>("predictor","consistent");
    break;
  case INPAR::STR::pred_constdisvelacc:
    genalphaparams->set<string>("predictor","constant");
    break;
  case INPAR::STR::pred_tangdis:
    genalphaparams->set<string>("predictor","tangdis");
    break;
  default:
    dserror("Cannot cope with choice of predictor");
    break;
  }

  genalphaparams->set<double>("UZAWAPARAM",sdyn.get<double>("UZAWAPARAM"));
  genalphaparams->set<double>("UZAWATOL",sdyn.get<double>("UZAWATOL"));
  genalphaparams->set<int>   ("UZAWAMAXITER",sdyn.get<int>("UZAWAMAXITER"));
  genalphaparams->set<INPAR::STR::ConSolveAlgo>("UZAWAALGO",getIntegralValue<INPAR::STR::ConSolveAlgo>(sdyn,"UZAWAALGO"));

  // sanity checks and default flags
  if (genprob.probtyp == prb_fsi)
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

    // robin flags
    INPAR::FSI::PartitionedCouplingMethod method =
      Teuchos::getIntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn,"PARTITIONED");
    genalphaparams->set<bool>  ("structrobin",
                                method==INPAR::FSI::DirichletRobin or method==INPAR::FSI::RobinRobin);

    genalphaparams->set<double>("alpha s",fsidyn.get<double>("ALPHA_S"));

    int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling == fsi_iter_monolithicfluidsplit or
        coupling == fsi_iter_monolithiclagrange or
        coupling == fsi_iter_monolithicstructuresplit)
    {
      if (Teuchos::getIntegralValue<INPAR::STR::PredEnum>(sdyn,"PREDICT")!=INPAR::STR::pred_constdisvelacc)
        dserror("only constant structure predictor with monolithic FSI possible");

#if 0
      // overwrite time integration flags
      genalphaparams->set<double>("gamma",fsidyn.get<double>("GAMMA"));
      genalphaparams->set<double>("alpha m",fsidyn.get<double>("ALPHA_M"));
      genalphaparams->set<double>("alpha f",fsidyn.get<double>("ALPHA_F"));
#endif
    }
  }

  RCP<Structure> tmpstr;
    tmpstr = Teuchos::rcp(new StructureGenAlpha(genalphaparams,actdis,solver,output));

  if (tmpstr->HaveConstraint())
  {
    structure_ = rcp(new StructureNOXCorrectionWrapper(
                       rcp(new StructureConstrained(tmpstr))));
  }
  else
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    if (Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO")==fsi_iter_monolithicxfem)
    {
      // no NOX correction required
      structure_ = tmpstr;
    }
    else
      structure_ = rcp(new StructureNOXCorrectionWrapper(tmpstr));
  }

  // invoke contact strugen alpha
  bool contact = false;
  switch (Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(scontact,"CONTACT"))
  {
    case INPAR::CONTACT::contact_none:
      contact = false;
      break;
    case INPAR::CONTACT::contact_normal:
      contact = true;
      break;
    case INPAR::CONTACT::contact_frictional:
      contact = true;
      break;
    case INPAR::CONTACT::contact_meshtying:
      contact = true;
      break;
    default:
      dserror("Cannot cope with choice of contact type");
      break;
  }

  if(contact)
  {
    structure_ = rcp(new StructureNOXCorrectionWrapper(
                       Teuchos::rcp(new ContactStructureGenAlpha(genalphaparams,actdis,solver,output))));
  }
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
  xparams->set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());

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
    INPAR::FSI::PartitionedCouplingMethod method
      = Teuchos::getIntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn,"PARTITIONED");
    xparams->set<bool>("structrobin",
                       ( (method==INPAR::FSI::DirichletRobin)
                         or (method==INPAR::FSI::RobinRobin) ));

    // THIS SHOULD GO, OR SHOULDN'T IT?
    xparams->set<double>("alpha s", fsidyn.get<double>("ALPHA_S"));

    // check if predictor fits to FSI algo
    int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
    if ( (coupling == fsi_iter_monolithicfluidsplit)
         or (coupling == fsi_iter_monolithiclagrange)
         or (coupling == fsi_iter_monolithicstructuresplit) )
    {
      if ((Teuchos::getIntegralValue<INPAR::STR::PredEnum>(*sdyn,"PREDICT")
          != INPAR::STR::pred_constdisvelacc) and
          (Teuchos::getIntegralValue<INPAR::STR::PredEnum>(*sdyn,"PREDICT")
          != INPAR::STR::pred_constdisvelaccpres))
      {
        dserror("only constant structure predictor with monolithic FSI possible");
      }
    }
  }

  // create a solver
  Teuchos::RCP<ParameterList> solveparams
    = Teuchos::rcp(new ParameterList());
  Teuchos::RCP<LINALG::Solver> solver
    = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->StructSolverParams(),
                                      actdis->Comm(),
                                      DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // create marching time integrator
  RCP<Structure> tmpstr;
  tmpstr = Teuchos::rcp(new StructureTimInt(ioflags, sdyn, xparams,
                                                actdis, solver, output));

  if (tmpstr->HaveConstraint())
    structure_ = rcp(new StructureNOXCorrectionWrapper(rcp(new StructureConstrained(tmpstr))));
  else
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    if (Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO")==fsi_iter_monolithicxfem)
    {
      // no NOX correction required
      structure_ = tmpstr;
    }
    else
      structure_ = rcp(new StructureNOXCorrectionWrapper(tmpstr));
  }

  // see you
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::Structure::Integrate()
{
  // a few parameters
  double time = GetTime();
  const double timeend = GetTimeEnd();
  const double timestepsize = GetTimeStepSize();
  int step = GetTimeStep();
  const int stepend = GetTimeNumStep();

  // loop ahead --- if timestepsize>0
  while ( (time < timeend) and (step < stepend) )
  {
    PrepareTimeStep();
    Solve();

    // update
    Update();
    time +=  timestepsize;
    step += 1;

    // talk to user
    fprintf(stdout,
            "Finalised: step %6d"
            " | nstep %6d"
            " | time %-14.8E"
            " | dt %-14.8E\n",
            step, stepend, time, timestepsize);
    // print a beautiful line made exactly of 80 dashes
    fprintf(stdout,
            "--------------------------------------------------------------"
            "------------------\n");
    // do it, print now!
    fflush(stdout);
    // talk to disk
    Output();
  }

  // Jump you f***ers
  return;
}

/*----------------------------------------------------------------------*/
#endif
