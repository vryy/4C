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
#include "adapter_structure_cmtstrugenalpha.H"
#include "adapter_structure_timint.H"
#include "adapter_structure_constr_merged.H"
#include "adapter_structure_wrapper.H"
#include "adapter_structure_lung.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_statmech.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "../drt_structure/strugenalpha.H"
#include "../drt_contact/contactstrugenalpha.H"
#include "../drt_contact/beam3contactstrugenalpha.H"
#include "../drt_contactnew/strugenalpha_cmt.H"
#include "../drt_statmech/statmech_time.H"

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
  const Teuchos::ParameterList& scontact = DRT::Problem::Instance()->MeshtyingAndContactParams();
  const Teuchos::ParameterList& statmech = DRT::Problem::Instance()->StatisticalMechanicsParams();

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

  INPAR::STR::ControlType controltype = Teuchos::getIntegralValue<INPAR::STR::ControlType>(sdyn,"CONTROLTYPE");
  genalphaparams->set<INPAR::STR::ControlType>("CONTROLTYPE",controltype);
  {
    vector<int> controlnode;
    std::istringstream contnode(Teuchos::getNumericStringParameter(sdyn,"CONTROLNODE"));
    std::string word;
    while (contnode >> word)
      controlnode.push_back(std::atoi(word.c_str()));
    if ((int)controlnode.size() != 3) dserror("Give proper values for CONTROLNODE in input file");
    genalphaparams->set("CONTROLNODE",controlnode[0]);
    genalphaparams->set("CONTROLDOF",controlnode[1]);
    genalphaparams->set("CONTROLCURVE",controlnode[2]);
  }

  {
    // use linearization of follower loads in Newton
    int loadlin = Teuchos::getIntegralValue<int>(sdyn,"LOADLIN");
    genalphaparams->set<bool>("LOADLIN",loadlin!=0);
  }

  genalphaparams->set<bool>  ("damping",(not (sdyn.get<std::string>("DAMPING") == "no"
                                              or sdyn.get<std::string>("DAMPING") == "No"
                                              or sdyn.get<std::string>("DAMPING") == "NO")));
  genalphaparams->set<double>("damping factor K",sdyn.get<double>("K_DAMP"));
  genalphaparams->set<double>("damping factor M",sdyn.get<double>("M_DAMP"));

  genalphaparams->set<double>("beta",sdyn.get<double>("BETA"));
#ifdef STRUGENALPHA_BE
  genalphaparams->set<double>("delta",sdyn.get<double>("DELTA"));
#endif
  genalphaparams->set<double>("gamma",sdyn.get<double>("GAMMA"));
  genalphaparams->set<double>("alpha m",sdyn.get<double>("ALPHA_M"));
  genalphaparams->set<double>("alpha f",sdyn.get<double>("ALPHA_F"));

  genalphaparams->set<double>("total time",0.0);
  genalphaparams->set<double>("delta time",prbdyn.get<double>("TIMESTEP"));
  genalphaparams->set<double>("max time",sdyn.get<double>("MAXTIME"));
  genalphaparams->set<int>   ("step",0);
  genalphaparams->set<int>   ("nstep",prbdyn.get<int>("NUMSTEP"));
  genalphaparams->set<int>   ("max iterations",sdyn.get<int>("MAXITER"));
  genalphaparams->set<int>   ("num iterations",-1);
  genalphaparams->set<string>("convcheck",sdyn.get<string>("CONV_CHECK"));
  genalphaparams->set<double>("tolerance displacements",sdyn.get<double>("TOLDISP"));
  genalphaparams->set<double>("tolerance residual",sdyn.get<double>("TOLRES"));
  genalphaparams->set<double>("tolerance constraint",sdyn.get<double>("TOLCONSTR"));

  genalphaparams->set<double>("UZAWAPARAM",sdyn.get<double>("UZAWAPARAM"));
  genalphaparams->set<double>("UZAWATOL",sdyn.get<double>("UZAWATOL"));
  genalphaparams->set<int>   ("UZAWAMAXITER",sdyn.get<int>("UZAWAMAXITER"));
  genalphaparams->set<INPAR::STR::ConSolveAlgo>("UZAWAALGO",getIntegralValue<INPAR::STR::ConSolveAlgo>(sdyn,"UZAWAALGO"));

  genalphaparams->set<bool>  ("io structural disp",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_DISP"));
  genalphaparams->set<int>   ("io disp every nstep",prbdyn.get<int>("UPRES"));

  genalphaparams->set<bool>  ("ADAPTCONV",getIntegralValue<int>(sdyn,"ADAPTCONV")==1);
  genalphaparams->set<double>("ADAPTCONV_BETTER",sdyn.get<double>("ADAPTCONV_BETTER"));

  INPAR::STR::StressType iostress = Teuchos::getIntegralValue<INPAR::STR::StressType>(ioflags,"STRUCT_STRESS");
  genalphaparams->set<INPAR::STR::StressType>("io structural stress", iostress);
  genalphaparams->set<int>   ("io stress every nstep",sdyn.get<int>("RESEVRYSTRS"));

  INPAR::STR::StrainType iostrain = Teuchos::getIntegralValue<INPAR::STR::StrainType>(ioflags,"STRUCT_STRAIN");
  genalphaparams->set<INPAR::STR::StrainType>("io structural strain", iostrain);

  genalphaparams->set<bool>  ("io surfactant",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_SURFACTANT"));

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


  // sanity checks and default flags
  if (genprob.probtyp == prb_fsi or genprob.probtyp == prb_fsi_lung)
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
        coupling == fsi_iter_monolithicstructuresplit or
        coupling == fsi_iter_lung_monolithicstructuresplit or
        coupling == fsi_iter_lung_monolithicfluidsplit or
        coupling == fsi_iter_constr_monolithicstructuresplit or
        coupling == fsi_iter_constr_monolithicfluidsplit)
    {
      if ((Teuchos::getIntegralValue<INPAR::STR::PredEnum>(sdyn,"PREDICT")!=INPAR::STR::pred_constdisvelacc)
          and (Teuchos::getIntegralValue<INPAR::STR::PredEnum>(sdyn,"PREDICT") != INPAR::STR::pred_constdisvelaccpres))
        dserror("only constant structure predictor with monolithic FSI possible");

#if 0
      // overwrite time integration flags
      genalphaparams->set<double>("gamma",fsidyn.get<double>("GAMMA"));
      genalphaparams->set<double>("alpha m",fsidyn.get<double>("ALPHA_M"));
      genalphaparams->set<double>("alpha f",fsidyn.get<double>("ALPHA_F"));
#endif
    }
  }

  // invoke contact or meshtying strugenalpha
  bool contact = false;
  bool mortarcontact = false;
  bool mortarmeshtying = false;
  bool beamcontact = false;
  bool oldcontact = false;
  INPAR::CONTACT::ApplicationType apptype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ApplicationType>(scontact,"APPLICATION");
  switch (apptype)
  {
    case INPAR::CONTACT::app_none:
      break;
    case INPAR::CONTACT::app_mortarcontact:
      mortarcontact = true;
      break;
    case INPAR::CONTACT::app_mortarmeshtying:
      mortarmeshtying = true;
      break;
    case INPAR::CONTACT::app_beamcontact:
      beamcontact = true;
      break;
    default:
      dserror("Cannot cope with choice of contact or meshtying type");
      break;
  }
  if (mortarcontact || mortarmeshtying || beamcontact) contact = true;

  // detect whether thermal bath is present
  bool thermalbath = false;
  switch (Teuchos::getIntegralValue<INPAR::STATMECH::ThermalBathType>(statmech,"THERMALBATH"))
  {
    case INPAR::STATMECH::thermalbath_none:
      thermalbath = false;
      break;
    case INPAR::STATMECH::thermalbath_uniform:
      thermalbath = true;
      break;
    case INPAR::STATMECH::thermalbath_shearflow:
      thermalbath = true;
      break;
    default:
      dserror("Cannot cope with choice of thermal bath");
      break;
  }

  Teuchos::RCP<StruGenAlpha> tintegrator = null;
  if (!mortarcontact && !mortarmeshtying && !beamcontact && 
      !oldcontact && !thermalbath)
    tintegrator = rcp(new StruGenAlpha(*genalphaparams,*actdis,*solver,*output));
  else
  {
    if (mortarcontact)
      tintegrator = rcp(new CONTACT::CmtStruGenAlpha(*genalphaparams,*actdis,*solver,*output,apptype));
    if (mortarmeshtying)
      tintegrator = rcp (new CONTACT::CmtStruGenAlpha(*genalphaparams,*actdis,*solver,*output,apptype));
    if (beamcontact)
      tintegrator = rcp(new CONTACT::Beam3ContactStruGenAlpha(*genalphaparams,*actdis,*solver,*output));
    if (thermalbath)
      tintegrator = rcp(new StatMechTime(*genalphaparams,*actdis,*solver,*output));
  }
  if (tintegrator == null) dserror("Failed to allocate strugenalpha derived time integrator");

  RCP<Structure> tmpstr;
  tmpstr = rcp(new StructureGenAlpha(genalphaparams,tintegrator,actdis,solver,output));

  if (tmpstr->HaveConstraint())
  {
    structure_ = rcp(new StructureNOXCorrectionWrapper(
                   rcp(new StructureConstrMerged(tmpstr))));
  }
  else
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling==fsi_iter_monolithicxfem) // no NOX correction required
      structure_ = tmpstr;
    else if (coupling == fsi_iter_lung_monolithicstructuresplit or
             coupling == fsi_iter_lung_monolithicfluidsplit)
      structure_ = rcp(new StructureLung(
                     rcp(new StructureNOXCorrectionWrapper(tmpstr))));
    else // everything else
      structure_ = rcp(new StructureNOXCorrectionWrapper(tmpstr));
  }

#if 0 // at some point, this is supposed to go away as it becomes part of 
      // ADAPTER::StructureGenAlpha
  if(contact)
  {
    tmpstr = null;
    structure_ = rcp(new StructureNOXCorrectionWrapper(
                   rcp(new CmtStructureGenAlpha(genalphaparams,actdis,solver,output,apptype))));
  }
#endif
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
  if (genprob.probtyp == prb_fsi or genprob.probtyp == prb_fsi_lung)
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
         or (coupling == fsi_iter_monolithicstructuresplit)
         or (coupling == fsi_iter_lung_monolithicstructuresplit)
         or (coupling == fsi_iter_lung_monolithicfluidsplit)
         or (coupling == fsi_iter_constr_monolithicfluidsplit)
         or (coupling == fsi_iter_constr_monolithicstructuresplit))
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
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");

  if (tmpstr->HaveConstraint())
    if (coupling == fsi_iter_constr_monolithicstructuresplit or 
        coupling == fsi_iter_constr_monolithicfluidsplit)
      structure_ = rcp(new StructureNOXCorrectionWrapper(tmpstr));
    else
      structure_ = rcp(new StructureNOXCorrectionWrapper(rcp(new StructureConstrMerged(tmpstr))));
  else
  {
    if (coupling==fsi_iter_monolithicxfem)
    {
      // no NOX correction required
      structure_ = tmpstr;
    }
    else if (coupling == fsi_iter_lung_monolithicstructuresplit or
             coupling == fsi_iter_lung_monolithicfluidsplit)
      structure_ = rcp(new StructureLung(rcp(new StructureNOXCorrectionWrapper(tmpstr))));
    else
      structure_ = rcp(new StructureNOXCorrectionWrapper(tmpstr));
  }

  // see you
  return;
}

/*----------------------------------------------------------------------*/
#endif
