/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_timint.cpp

\brief Structure field adapter

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_structure/strtimint_create.H"
#include "adapter_structure_timint.H"
#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*/
//! General problem data
//!
//! global variable GENPROB genprob is defined in global_control.c
//! \author m.gee \date 06/01
extern struct _GENPROB genprob;

/*----------------------------------------------------------------------*/
//! global variable *solv, vector of lenght numfld of structures SOLVAR
//! defined in solver_control.c
//!
//! \author m.gee \date 11/00
extern struct _SOLVAR *solv;

/*----------------------------------------------------------------------*/
//! File pointers
//!
//! This structure struct _FILES allfiles is defined in input_control_global.c
//! and the type is in standardtypes.h
//! It holds all file pointers and some variables needed for the FRSYSTEM
extern struct _FILES allfiles;


/*======================================================================*/
/* constructor */
ADAPTER::StructureTimInt::StructureTimInt(
  Teuchos::RCP<Teuchos::ParameterList> ioparams,
  Teuchos::RCP<Teuchos::ParameterList> sdynparams,
  Teuchos::RCP<Teuchos::ParameterList> xparams,
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: structure_(STR::TimIntImplCreate(*ioparams, *sdynparams, *xparams,
                                   discret, solver, output)),
  discret_(discret),
  ioparams_(ioparams),
  sdynparams_(sdynparams),
  xparams_(xparams),
  solver_(solver),
  output_(output)
{
  // make sure
  if (structure_ == Teuchos::null)
    dserror("Failed to create structural integrator");

  // set-up FSI interface
  DRT::UTILS::SetupNDimExtractor(*discret, "FSICoupling", interface_);
  structure_->SetSurfaceFSI(&interface_);
}


/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimInt::InitialGuess()
{
  return structure_->DisRes();
  //return structure_->Dispm();
}


/*----------------------------------------------------------------------*/
/* right-hand side alias the dynamic force residual */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimInt::RHS()
{
  return structure_->ForceRes();
}


/*----------------------------------------------------------------------*/
/* get current displacements D_{n+1} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimInt::Dispnp()
{
#ifdef INVERSEDESIGNCREATE
  dserror("check this");
  return Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
#else
  return structure_->DisNew();
#endif
}


/*----------------------------------------------------------------------*/
/* get last converged displacements D_{n} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimInt::Dispn()
{
#ifdef INVERSEDESIGNCREATE
  dserror("check this");
  return Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
#else
  return structure_->Dis();
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// UNWANTED
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimInt::Dispnm()
{
/*
#ifdef INVERSEDESIGNCREATE
  return Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
#else
  return structure_->Dispm();
#endif
*/
  return Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
}


/*----------------------------------------------------------------------*/
/* non-overlapping DOF map */
Teuchos::RCP<const Epetra_Map> ADAPTER::StructureTimInt::DofRowMap()
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}


/*----------------------------------------------------------------------*/
/* stiffness, i.e. force residual R_{n+1} differentiated
 * by displacements D_{n+1} */
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::StructureTimInt::SystemMatrix()
{
  return structure_->Stiff();
}


/*----------------------------------------------------------------------*/
/* get discretisation */
Teuchos::RCP<DRT::Discretization> ADAPTER::StructureTimInt::Discretization()
{
  return structure_->Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// UNWANTED
double ADAPTER::StructureTimInt::DispIncrFactor()
{
  return -1; //structure_->DispIncrFactor();
}

/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimInt::FRobin()
{
  //return structure_->GetForceRobinFSI();
  return LINALG::CreateVector(*discret_->DofRowMap(), true);
}

/*----------------------------------------------------------------------*/
/* External force F_{ext,n+1} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimInt::FExtn()
{
  return structure_->FextNew();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimInt::FluidCondRHS() const
// {
//   // structure part of the rhs to enforce
//   // u(n+1) dt = d(n+1) - d(n)

//   // extrapolate d(n+1) at the interface and substract d(n)

//   Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractCondVector(structure_->Dispm());
//   Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_->Disp ());

//   double alphaf = structure_->AlphaF();
//   idis->Update(1./(1.-alphaf), *idism, -alphaf/(1.-alphaf)-1.);
//   return idis;
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimInt::MeshCondRHS() const
// {
//   // structure part of the rhs to enforce
//   // d(G,n+1) = d(n+1)

//   // extrapolate d(n+1) at the interface

//   Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractCondVector(structure_->Dispm());
//   Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_->Disp ());

//   double alphaf = structure_->AlphaF();
//   idis->Update(1./(1.-alphaf), *idism, -alphaf/(1.-alphaf));
//   return idis;
// }


/*----------------------------------------------------------------------*/
/* prepare time step */
void ADAPTER::StructureTimInt::PrepareTimeStep()
{
  // Note: MFSI requires a constant predictor. Otherwise the fields will get
  // out of sync.

  // predict
  structure_->Predict();

  // initialise incremental displacements
  if (disinc_ != Teuchos::null)
    disinc_->PutScalar(0.0);

}


/*----------------------------------------------------------------------*/
/* build linear system stiffness matrix and rhs/force residual
 *
 * Monolithic FSI accesses the linearised structure problem. */
void ADAPTER::StructureTimInt::Evaluate(
  Teuchos::RCP<const Epetra_Vector> disp
)
{
  // Yes, this is complicated. But we have to be very careful
  // here. The field solver always expects an increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest
  // increment only.
  if (disp != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(new Epetra_Vector(*disp));
    if (disinc_ != Teuchos::null)
    {
      disi->Update(-1.0, *disinc_, 1.0);
      disinc_->Update(1.0, *disp, 0.0);
    }
    else
    {
      disinc_ = Teuchos::rcp(new Epetra_Vector(*disp));
    }
    structure_->UpdateIterIncrementally(disi);
    structure_->EvaluateForceStiffResidual();
  }
  else
  {
    // ORIGINAL ACTION BUILDS TANGENT AND APPLIES DBC ON SLE
    // SENSIBLE ???
    structure_->UpdateIterIncrementally(Teuchos::null);
    structure_->EvaluateForceStiffResidual();
  }
}


/*----------------------------------------------------------------------*/
/* update time step */
void ADAPTER::StructureTimInt::Update()
{
  structure_->UpdateStepAndTime();
}


/*----------------------------------------------------------------------*/
/* output */
void ADAPTER::StructureTimInt::Output()
{
  structure_->OutputStep();
}


/*----------------------------------------------------------------------*/
/* domain map */
const Epetra_Map& ADAPTER::StructureTimInt::DomainMap()
{
  return structure_->GetDomainMap();
}


/*----------------------------------------------------------------------*/
/* read restart */
void ADAPTER::StructureTimInt::ReadRestart(int step)
{
  structure_->ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/* find iteratively solution */
void ADAPTER::StructureTimInt::Solve()
{
  structure_->Solve();
}


/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimInt::RelaxationSolve(
  Teuchos::RCP<Epetra_Vector> iforce
)
{
  Teuchos::RCP<Epetra_Vector> relax = interface_.InsertCondVector(iforce);
  structure_->SetForceInterface(relax);
  Teuchos::RCP<Epetra_Vector> idisi = structure_->SolveRelaxationLinear();

  // we are just interested in the incremental interface displacements
  idisi = interface_.ExtractCondVector(idisi);
  return idisi;
}

/*----------------------------------------------------------------------*/
/* extract interface displacements D_{n} */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimInt::ExtractInterfaceDispn()
{
#ifdef INVERSEDESIGNCREATE
  dserror("Check this");
  return Teuchos::rcp(new Epetra_Vector(*interface_.CondMap()));
#else
  Teuchos::RCP<Epetra_Vector> idis
    = interface_.ExtractCondVector(structure_->Dis());
  return idis;
#endif
}

/*----------------------------------------------------------------------*/
/* extract interface displacements D_{n+1} */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimInt::ExtractInterfaceDispnp()
{
#ifdef INVERSEDESIGNCREATE
  dserror("Check this");
  return Teuchos::rcp(new Epetra_Vector(*interface_.CondMap()));
#else
  Teuchos::RCP<Epetra_Vector> idis
    = interface_.ExtractCondVector(structure_->DisNew());
  return idis;
#endif
}

/*----------------------------------------------------------------------*/
/* extract external forces at interface F_{ext,n+1} */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimInt::ExtractInterfaceForces()
{
  Teuchos::RCP<Epetra_Vector> iforce
    = interface_.ExtractCondVector(structure_->FextNew());
  return iforce;
}

/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimInt::PredictInterfaceDispnp()
{
  const Teuchos::ParameterList& fsidyn
    = DRT::Problem::Instance()->FSIDynamicParams();

  Teuchos::RCP<Epetra_Vector> idis;

  switch (Teuchos::getIntegralValue<int>(fsidyn,"PREDICTOR"))
  {
  case 1:
  {
    // d(n)
    // respect Dirichlet conditions at the interface (required for pseudo-rigid body)
    idis  = interface_.ExtractCondVector(structure_->DisNew());
    break;
  }
  case 2:
    // d(n)+dt*(1.5*v(n)-0.5*v(n-1))
    dserror("interface velocity v(n-1) not available");
    break;
  case 3:
  {
    // d(n)+dt*v(n)
    double dt = sdynparams_->get<double>("TIMESTEP");

    idis = interface_.ExtractCondVector(structure_->Dis());
    Teuchos::RCP<Epetra_Vector> ivel
      = interface_.ExtractCondVector(structure_->Vel());

    idis->Update(dt,* ivel, 1.0);
    break;
  }
  case 4:
  {
    // d(n)+dt*v(n)+0.5*dt^2*a(n)
    double dt = sdynparams_->get<double>("TIMESTEP");

    idis = interface_.ExtractCondVector(structure_->Dis());
    Teuchos::RCP<Epetra_Vector> ivel
      = interface_.ExtractCondVector(structure_->Vel());
    Teuchos::RCP<Epetra_Vector> iacc
      = interface_.ExtractCondVector(structure_->Acc());

    idis->Update(dt, *ivel, 0.5*dt*dt, *iacc, 1.0);
    break;
  }
  default:
    dserror("unknown interface displacement predictor '%s'",
            fsidyn.get<string>("PREDICTOR").c_str());
  }

  return idis;
}


/*----------------------------------------------------------------------*/
/* */
void ADAPTER::StructureTimInt::ApplyInterfaceForces(
  Teuchos::RCP<Epetra_Vector> iforce
)
{
/*
  // Play it save. In the first iteration everything is already set up
  // properly. However, all following iterations need to calculate the
  // stiffness matrix here. Furthermore we are bound to reset fextm_
  // before we add our special contribution.
  // So we calculate the stiffness anyway (and waste the available
  // stiffness in the first iteration).
  structure_->ApplyExternalForce(interface_,iforce);
*/
  // This will add the provided interface force onto the residual forces
  // The sign convention of the interface force is external-force-like.
  structure_->SetForceInterface(interface_, iforce);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimInt::ApplyInterfaceRobinValue(
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> ifluidvel
)
{
  dserror("Not impl.");
/*
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

  Teuchos::RCP<Epetra_Vector> idisn  = interface_.ExtractCondVector(structure_->Disp());
  Teuchos::RCP<Epetra_Vector> frobin = interface_.ExtractCondVector(structure_->FRobin());

  // save robin coupling values in frobin vector (except iforce which
  // is passed separately)
  frobin->Update(alphas/dt,*idisn,alphas*(1-alphaf),*ifluidvel,0.0);

  interface_.InsertCondVector(frobin,structure_->FRobin());
  structure_->ApplyExternalForce(interface_,iforce);
*/
}


/*----------------------------------------------------------------------*/
/* structural result test */
Teuchos::RCP<DRT::ResultTest> ADAPTER::StructureTimInt::CreateFieldTest()
{
  return Teuchos::rcp(new StruResultTest(*structure_));
}


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
