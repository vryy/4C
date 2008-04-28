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
ADAPTER::StructureGenAlpha::StructureGenAlpha(Teuchos::RCP<Teuchos::ParameterList> params,
                                              Teuchos::RCP<DRT::Discretization> dis,
                                              Teuchos::RCP<LINALG::Solver> solver,
                                              Teuchos::RCP<IO::DiscretizationWriter> output)
  : structure_(*params, *dis, *solver, *output),
    dis_(dis),
    params_(params),
    solver_(solver),
    output_(output)
{
  UTILS::SetupNDimExtractor(*dis,"FSICoupling",interface_);
  structure_.SetFSISurface(&interface_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::InitialGuess()
{
  return Teuchos::rcp(&structure_.Getdu(),false);
  //return structure_.Dispm();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::RHS()
{
  return structure_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::Dispnp()
{
  double alphaf = structure_.AlphaF();
  Teuchos::RCP<Epetra_Vector> dispnp = Teuchos::rcp(new Epetra_Vector(*Dispn()));
  dispnp->Update(1./(1.-alphaf),*Dispnm(),-alphaf/(1.-alphaf));
  return dispnp;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::Dispn()
{
  return structure_.Disp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::Dispnm()
{
  return structure_.Dispm();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::StructureGenAlpha::DofRowMap()
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::StructureGenAlpha::SystemMatrix()
{
  return structure_.SystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::StructureGenAlpha::Discretization()
{
  return structure_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::StructureGenAlpha::DispIncrFactor()
{
  return structure_.DispIncrFactor();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::FRobin()
{
  return structure_.FRobin();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::FExtn()
{
  return structure_.FExtn();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::FluidCondRHS() const
// {
//   // structure part of the rhs to enforce
//   // u(n+1) dt = d(n+1) - d(n)

//   // extrapolate d(n+1) at the interface and substract d(n)

//   Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractCondVector(structure_.Dispm());
//   Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_.Disp ());

//   double alphaf = structure_.AlphaF();
//   idis->Update(1./(1.-alphaf), *idism, -alphaf/(1.-alphaf)-1.);
//   return idis;
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::MeshCondRHS() const
// {
//   // structure part of the rhs to enforce
//   // d(G,n+1) = d(n+1)

//   // extrapolate d(n+1) at the interface

//   Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractCondVector(structure_.Dispm());
//   Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_.Disp ());

//   double alphaf = structure_.AlphaF();
//   idis->Update(1./(1.-alphaf), *idism, -alphaf/(1.-alphaf));
//   return idis;
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::PrepareTimeStep()
{
  // Note: MFSI requires a constant predictor. Otherwise the fields will get
  // out of sync.

  std::string pred = params_->get<string>("predictor","consistent");
  if (pred=="constant")
  {
    structure_.ConstantPredictor();
  }
  else if (pred=="consistent")
  {
    structure_.ConsistentPredictor();
  }
  else
    dserror("predictor %s unknown", pred.c_str());

  if (sumdisi_!=Teuchos::null)
    sumdisi_->PutScalar(0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::Evaluate(Teuchos::RCP<const Epetra_Vector> disp)
{
  // Yes, this is complicated. But we have to be very careful
  // here. The field solver always expects an increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest
  // increment only.
  if (disp!=Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(new Epetra_Vector(*disp));
    if (sumdisi_!=Teuchos::null)
    {
      disi->Update(-1.0,*sumdisi_,1.0);
      sumdisi_->Update(1.0,*disp,0.0);
    }
    else
    {
      sumdisi_ = Teuchos::rcp(new Epetra_Vector(*disp));
    }
    structure_.Evaluate(disi);
  }
  else
  {
    structure_.Evaluate(Teuchos::null);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::Update()
{
  structure_.UpdateandOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::Output()
{
  // noop
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Epetra_Map& ADAPTER::StructureGenAlpha::DomainMap()
{
  return structure_.DomainMap();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::ReadRestart(int step)
{
  structure_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::Solve()
{
  std::string equil = params_->get<string>("equilibrium iteration","full newton");

  if (structure_.HaveConstraint())
  {
    structure_.FullNewtonLinearUzawa();
  }
  else if (equil=="full newton")
  {
    structure_.FullNewton();
  }
  else if (equil=="modified newton")
  {
    structure_.ModifiedNewton();
  }
  else if (equil=="nonlinear cg")
  {
    structure_.NonlinearCG();
  }
  else if (equil=="ptc")
  {
    structure_.PTC();
  }
  else
    dserror("Unknown type of equilibrium iteration '%s'", equil.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce)
{
  Teuchos::RCP<Epetra_Vector> relax = interface_.InsertCondVector(iforce);
  Teuchos::RCP<Epetra_Vector> idisi = structure_.LinearRelaxationSolve(relax);

  // we are just interested in the incremental interface displacements
  idisi = interface_.ExtractCondVector(idisi);
  return idisi;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::ExtractInterfaceDispn()
{
  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_.Disp());
  return idis;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::ExtractInterfaceDispnp()
{
  Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractCondVector(structure_.Dispm());
  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_.Disp());

  double alphaf = params_->get<double>("alpha f", 0.459);
  idis->Update(1./(1.-alphaf),*idism,-alphaf/(1.-alphaf));

  return idis;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::ExtractInterfaceForces()
{
  Teuchos::RCP<Epetra_Vector> iforce = interface_.ExtractCondVector(structure_.FExtn());

  return iforce;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::PredictInterfaceDispnp()
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_.Disp());

  switch (Teuchos::getIntegralValue<int>(fsidyn,"PREDICTOR"))
  {
  case 1:
    // d(n)
    // nothing to do
    break;
  case 2:
    // d(n)+dt*(1.5*v(n)-0.5*v(n-1))
    dserror("interface velocity v(n-1) not available");
    break;
  case 3:
  {
    // d(n)+dt*v(n)
    double dt            = params_->get<double>("delta time"             ,0.01);

    Teuchos::RCP<Epetra_Vector> ivel  = interface_.ExtractCondVector(structure_.Vel());

    idis->Update(dt,*ivel,1.0);
    break;
  }
  case 4:
  {
    // d(n)+dt*v(n)+0.5*dt^2*a(n)
    double dt            = params_->get<double>("delta time"             ,0.01);

    Teuchos::RCP<Epetra_Vector> ivel  = interface_.ExtractCondVector(structure_.Vel());
    Teuchos::RCP<Epetra_Vector> iacc  = interface_.ExtractCondVector(structure_.Acc());

    idis->Update(dt,*ivel,0.5*dt*dt,*iacc,1.0);
    break;
  }
  default:
    dserror("unknown interface displacement predictor '%s'",
            fsidyn.get<string>("PREDICTOR").c_str());
  }

  return idis;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce)
{
  // Play it save. In the first iteration everything is already set up
  // properly. However, all following iterations need to calculate the
  // stiffness matrix here. Furthermore we are bound to reset fextm_
  // before we add our special contribution.
  // So we calculate the stiffness anyway (and waste the available
  // stiffness in the first iteration).
  structure_.ApplyExternalForce(interface_,iforce);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> iforce,
                                                          Teuchos::RCP<Epetra_Vector> ifluidvel)
{
  // get robin parameter and timestep
  double alphas  = params_->get<double>("alpha s",-1.);
  double dt      = params_->get<double>("delta time",-1.);
  double alphaf  = params_->get<double>("alpha f", 0.459);

  if (alphas<0. or dt<0.)
    dserror("couldn't get robin parameter alpha_s or time step size");

  // the RobinRHS is going to be:
  //
  // RobinRHS =
  //     - (alpha_s/dt)*(dis(n))
  //     - alpha_s*(1-alpha_f)*(fluidvel(n+1))
  //     + (1-alpha_f)*(iforce(n+1))

  // Attention: We must not change iforce here, because we would
  // implicitely change fextn_, too. fextn_ is needed to set fext_
  // after successfully reaching timestep end.
  // This is why an additional robin force vector is needed.

  Teuchos::RCP<Epetra_Vector> idisn  = interface_.ExtractCondVector(structure_.Disp());
  Teuchos::RCP<Epetra_Vector> frobin = interface_.ExtractCondVector(structure_.FRobin());

  // save robin coupling values in frobin vector (except iforce which
  // is passed separately)
  frobin->Update(alphas/dt,*idisn,alphas*(1-alphaf),*ifluidvel,0.0);

  interface_.InsertCondVector(frobin,structure_.FRobin());
  structure_.ApplyExternalForce(interface_,iforce);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::StructureGenAlpha::CreateFieldTest()
{
  return Teuchos::rcp(new StruResultTest(structure_));
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

  const Teuchos::ParameterList& size     = DRT::Problem::Instance()->ProblemSizeParams();

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

    if (Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO") == fsi_iter_monolithic)
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

#endif
