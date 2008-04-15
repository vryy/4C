#ifdef CCADISCRET

#include "adapter_fluid.H"

// further includes for FluidBaseAlgorithm:
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::Fluid::~Fluid()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidAdapter::FluidAdapter(Teuchos::RCP<DRT::Discretization> dis,
                                 Teuchos::RCP<LINALG::Solver> solver,
                                 Teuchos::RCP<ParameterList> params,
                                 Teuchos::RCP<IO::DiscretizationWriter> output,
                                 bool isale)
  : fluid_(dis, *solver, *params, *output, isale),
    dis_(dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  UTILS::SetupNDimExtractor(*dis,"FSICoupling",interface_);
  UTILS::SetupNDimExtractor(*dis,"FREESURFCoupling",freesurface_);

  fluid_.SetFreeSurface(&freesurface_);

  // build inner velocity map
  // dofs at the interface are excluded
  // we use only velocity dofs and only those without Dirichlet constraint

  Teuchos::RCP<const Epetra_Map> velmap = fluid_.VelocityRowMap();
  Teuchos::RCP<Epetra_Vector> dirichtoggle = fluid_.Dirichlet();
  Teuchos::RCP<const Epetra_Map> fullmap = DofRowMap();

  int numvelids = velmap->NumMyElements();
  std::vector<int> velids;
  velids.reserve(numvelids);
  for (int i=0; i<numvelids; ++i)
  {
    int gid = velmap->GID(i);
    if (not interface_.CondMap()->MyGID(gid) and (*dirichtoggle)[fullmap->LID(gid)]==0.)
    {
      velids.push_back(gid);
    }
  }

  innervelmap_ = Teuchos::rcp(new Epetra_Map(-1,velids.size(), &velids[0], 0, velmap->Comm()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidAdapter::InitialGuess() const
{
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidAdapter::RHS() const
{
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidAdapter::Velnp() const
{
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidAdapter::Veln() const
{
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidAdapter::Dispnp() const
{
  return fluid_.Dispnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidAdapter::DofRowMap() const
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::FluidAdapter::SystemMatrix() const
{
  return fluid_.SystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidAdapter::Discretization()
{
  return fluid_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAdapter::StructCondRHS() const
// {
//   return interface_.ExtractCondVector(Velnp());
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::PrepareTimeStep()
{
  fluid_.PrepareTimeStep();

  // we add the whole fluid mesh displacement later on?
  //fluid_.Dispnp()->PutScalar(0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::Evaluate(Teuchos::RCP<const Epetra_Vector> vel) const
{
  if (vel!=Teuchos::null)
  {
    fluid_.Evaluate(vel);
  }
  else
  {
    fluid_.Evaluate(Teuchos::null);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::Update()
{
  fluid_.TimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::Output()
{
  fluid_.Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::NonlinearSolve()
{
  fluid_.NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidAdapter::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidAdapter::PressureRowMap()
{
  return fluid_.PressureRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidAdapter::ResidualScaling() const
{
  return fluid_.ResidualScaling();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidAdapter::TimeScaling() const
{
  return 1./fluid_.Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::ReadRestart(int step)
{
  fluid_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidAdapter::Time()
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::FluidAdapter::Step()
{
  return fluid_.Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::LiftDrag()
{
  fluid_.LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAdapter::ExtractInterfaceForces()
{
  return interface_.ExtractCondVector(fluid_.TrueResidual());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  interface_.InsertCondVector(ivel,fluid_.Velnp());

  // mark all interface velocities as dirichlet values
  // this is very easy, but there are two dangers:
  // - We change ivel here. It must not be used afterwards.
  // - The algorithm must support the sudden change of dirichtoggle_
  ivel->PutScalar(1.0);
  interface_.InsertCondVector(ivel,fluid_.Dirichlet());

  //----------------------- compute an inverse of the dirichtoggle vector
  fluid_.InvDirichlet()->PutScalar(1.0);
  fluid_.InvDirichlet()->Update(-1.0,*fluid_.Dirichlet(),1.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp) const
{
  meshmap_.InsertCondVector(fluiddisp,fluid_.Dispnp());

  // new grid velocity
  fluid_.UpdateGridv();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel) const
{
  meshmap_.InsertCondVector(gridvel,fluid_.GridVel());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::ConvertInterfaceUnknown(Teuchos::RCP<Epetra_Vector> fcx) const
{
  // get interface velocity at t(n)
  Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractCondVector(Veln());

  // We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( Delta u(n+1,i+1) + u(n) ) * dt
  //
  fcx->Update(-1.,*veln,TimeScaling());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::FluidAdapter::Itemax() const
{
  return fluid_.Itemax();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdapter::SetItemax(int itemax)
{
  fluid_.SetItemax(itemax);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAdapter::IntegrateInterfaceShape()
{
  return interface_.ExtractCondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAdapter::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
  interface_.InsertCondVector(ivel,relax);
  fluid_.LinearRelaxationSolve(relax);
  return ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidAdapter::CreateFieldTest()
{
  return Teuchos::rcp(new FluidResultTest(fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidAdapter::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres) const
{
   return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres);
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::XFluidAdapter::XFluidAdapter(Teuchos::RCP<DRT::Discretization> dis,
                                 Teuchos::RCP<LINALG::Solver> solver,
                                 Teuchos::RCP<ParameterList> params,
                                 Teuchos::RCP<IO::DiscretizationWriter> output,
                                 bool isale)
  : fluid_(dis, *solver, *params, *output, isale),
    dis_(dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  UTILS::SetupNDimExtractor(*dis,"FSICoupling",interface_);
  UTILS::SetupNDimExtractor(*dis,"FREESURFCoupling",freesurface_);

  fluid_.SetFreeSurface(&freesurface_);

  // build inner velocity map
  // dofs at the interface are excluded
  // we use only velocity dofs and only those without Dirichlet constraint

  Teuchos::RCP<const Epetra_Map> velmap = fluid_.VelocityRowMap();
  Teuchos::RCP<Epetra_Vector> dirichtoggle = fluid_.Dirichlet();
  Teuchos::RCP<const Epetra_Map> fullmap = DofRowMap();

  int numvelids = velmap->NumMyElements();
  std::vector<int> velids;
  velids.reserve(numvelids);
  for (int i=0; i<numvelids; ++i)
  {
    int gid = velmap->GID(i);
    if (not interface_.CondMap()->MyGID(gid) and (*dirichtoggle)[fullmap->LID(gid)]==0.)
    {
      velids.push_back(gid);
    }
  }

  innervelmap_ = Teuchos::rcp(new Epetra_Map(-1,velids.size(), &velids[0], 0, velmap->Comm()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidAdapter::InitialGuess() const
{
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidAdapter::RHS() const
{
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidAdapter::Velnp() const
{
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidAdapter::Veln() const
{
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidAdapter::Dispnp() const
{
  return fluid_.Dispnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluidAdapter::DofRowMap() const
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::XFluidAdapter::SystemMatrix() const
{
  return fluid_.SystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::XFluidAdapter::Discretization()
{
  return fluid_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidAdapter::StructCondRHS() const
// {
//   return interface_.ExtractCondVector(Velnp());
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::PrepareTimeStep()
{
  fluid_.PrepareTimeStep();

  // we add the whole fluid mesh displacement later on?
  //fluid_.Dispnp()->PutScalar(0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::Evaluate(Teuchos::RCP<const Epetra_Vector> vel) const
{
  if (vel!=Teuchos::null)
  {
    fluid_.Evaluate(vel);
  }
  else
  {
    fluid_.Evaluate(Teuchos::null);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::Update()
{
  fluid_.TimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::Output()
{
  fluid_.Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::NonlinearSolve()
{
  fluid_.NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluidAdapter::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluidAdapter::PressureRowMap()
{
  return fluid_.PressureRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluidAdapter::ResidualScaling() const
{
  return fluid_.ResidualScaling();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluidAdapter::TimeScaling() const
{
  return 1./fluid_.Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::ReadRestart(int step)
{
  fluid_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluidAdapter::Time()
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::XFluidAdapter::Step()
{
  return fluid_.Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::LiftDrag()
{
  fluid_.LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidAdapter::ExtractInterfaceForces()
{
  return interface_.ExtractCondVector(fluid_.TrueResidual());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  interface_.InsertCondVector(ivel,fluid_.Velnp());

  // mark all interface velocities as dirichlet values
  // this is very easy, but there are two dangers:
  // - We change ivel here. It must not be used afterwards.
  // - The algorithm must support the sudden change of dirichtoggle_
  ivel->PutScalar(1.0);
  interface_.InsertCondVector(ivel,fluid_.Dirichlet());

  //----------------------- compute an inverse of the dirichtoggle vector
  fluid_.InvDirichlet()->PutScalar(1.0);
  fluid_.InvDirichlet()->Update(-1.0,*fluid_.Dirichlet(),1.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp) const
{
  meshmap_.InsertCondVector(fluiddisp,fluid_.Dispnp());

  // new grid velocity
  fluid_.UpdateGridv();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel) const
{
  meshmap_.InsertCondVector(gridvel,fluid_.GridVel());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::ConvertInterfaceUnknown(Teuchos::RCP<Epetra_Vector> fcx) const
{
  // get interface velocity at t(n)
  Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractCondVector(Veln());

  // We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( Delta u(n+1,i+1) + u(n) ) * dt
  //
  fcx->Update(-1.,*veln,TimeScaling());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::XFluidAdapter::Itemax() const
{
  return fluid_.Itemax();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidAdapter::SetItemax(int itemax)
{
  fluid_.SetItemax(itemax);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidAdapter::IntegrateInterfaceShape()
{
  return interface_.ExtractCondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidAdapter::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
  interface_.InsertCondVector(ivel,relax);
  fluid_.LinearRelaxationSolve(relax);
  return ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::XFluidAdapter::CreateFieldTest()
{
  return Teuchos::rcp(new FluidResultTest(fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidAdapter::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres) const
{
   return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres);
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

ADAPTER::FluidGenAlphaAdapter::FluidGenAlphaAdapter(
  Teuchos::RCP<DRT::Discretization>      dis,
  Teuchos::RCP<LINALG::Solver>           solver,
  Teuchos::RCP<ParameterList>            params,
  Teuchos::RCP<IO::DiscretizationWriter> output,
  bool                                   isale)
  : fluid_ (dis, *solver, *params, *output, isale),
    dis_   (dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  UTILS::SetupNDimExtractor(*dis,"FSICoupling",interface_);
  UTILS::SetupNDimExtractor(*dis,"FREESURFCoupling",freesurface_);

  // build inner velocity map
  // dofs at the interface are excluded
  // we use only velocity dofs and only those without Dirichlet constraint

  Teuchos::RCP<const Epetra_Map>    velmap       = fluid_.VelocityRowMap();
  Teuchos::RCP<const Epetra_Vector> dirichtoggle = fluid_.Dirichlet();
  Teuchos::RCP<const Epetra_Map>    fullmap      = DofRowMap();

  int numvelids = velmap->NumMyElements();
  std::vector<int> velids;
  velids.reserve(numvelids);
  for (int i=0; i<numvelids; ++i)
  {
    int gid = velmap->GID(i);
    if (not interface_.CondMap()->MyGID(gid) and (*dirichtoggle)[fullmap->LID(gid)]==0.)
    {
      velids.push_back(gid);
    }
  }

  innervelmap_ = Teuchos::rcp(new Epetra_Map(-1,velids.size(), &velids[0], 0, velmap->Comm()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidGenAlphaAdapter::InitialGuess() const
{
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidGenAlphaAdapter::RHS() const
{
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidGenAlphaAdapter::Velnp() const
{
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidGenAlphaAdapter::Veln() const
{
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidGenAlphaAdapter::Dispnp() const
{
  return fluid_.Dispnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidGenAlphaAdapter::DofRowMap() const
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::FluidGenAlphaAdapter::SystemMatrix() const
{
  return fluid_.SysMat();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidGenAlphaAdapter::Discretization()
{
  return fluid_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::PrepareTimeStep()
{
  fluid_.GenAlphaIncreaseTimeAndStep();

  fluid_.GenAlphaEchoToScreen("print time algorithm info");
  fluid_.GenAlphaPredictNewSolutionValues();
  fluid_.GenAlphaApplyDirichletAndNeumann();
  fluid_.GenAlphaCalcInitialAccelerations();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::Evaluate(Teuchos::RCP<const Epetra_Vector> dacc) const
{
  if (dacc!=Teuchos::null)
  {
    fluid_.ExternIncrementOfVelnp(dacc);
  }

  fluid_.GenAlphaComputeIntermediateSol();
  fluid_.GenAlphaAssembleResidualAndMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::Update()
{
  fluid_.GenAlphaTimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::Output()
{
  fluid_.GenAlphaOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::NonlinearSolve()
{
  fluid_.DoGenAlphaPredictorCorrectorIteration();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidGenAlphaAdapter::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidGenAlphaAdapter::PressureRowMap()
{
  return fluid_.PressureRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidGenAlphaAdapter::ResidualScaling() const
{
  return fluid_.ResidualScaling();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidGenAlphaAdapter::TimeScaling() const
{
  double dt = fluid_.Dt();
  double gamma = fluid_.Gamma();
  return 1./(dt*dt*gamma);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::ReadRestart(int step)
{
  fluid_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidGenAlphaAdapter::Time()
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::FluidGenAlphaAdapter::Step()
{
  return fluid_.Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::LiftDrag()
{
  fluid_.LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidGenAlphaAdapter::ExtractInterfaceForces()
{
  return interface_.ExtractCondVector(fluid_.TrueResidual());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  interface_.InsertCondVector(ivel,fluid_.Velnp());

  // mark all interface velocities as dirichlet values
  // this is very easy, but there are two dangers:
  // - We change ivel here. It must not be used afterwards.
  // - The algorithm must support the sudden change of dirichtoggle_
  ivel->PutScalar(1.0);
  interface_.InsertCondVector(ivel,fluid_.Dirichlet());

  //----------------------- compute an inverse of the dirichtoggle vector
  fluid_.InvDirichlet()->PutScalar(1.0);
  fluid_.InvDirichlet()->Update(-1.0,*fluid_.Dirichlet(),1.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp) const
{
  meshmap_.InsertCondVector(fluiddisp,fluid_.Dispnp());

  // new grid velocity
  fluid_.UpdateGridv();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel) const
{
  meshmap_.InsertCondVector(gridvel,fluid_.GridVel());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::ConvertInterfaceUnknown(Teuchos::RCP<Epetra_Vector> fcx) const
{
  // We convert Delta d(n+1,i+1) to Delta a(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) =   ( Delta u(n+1,i+1)              + u(n) ) * dt
  //
  //                  = ( ( Delta a(n+1,i+1) * gamma * dt + u(n) ) * dt
  //

  double dt = fluid_.Dt();
  double gamma = fluid_.Gamma();

  // get interface velocity at t(n)
  Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractCondVector(fluid_.Veln());

  // reduce to Delta u(n+1,i+1)
  fcx->Update(-1.,*veln,1./dt);

  // reduce to Delta a(n+1,i+1)
  fcx->Scale(1./dt/gamma);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::FluidGenAlphaAdapter::Itemax() const
{
  return fluid_.Itemax();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlphaAdapter::SetItemax(int itemax)
{
  fluid_.SetItemax(itemax);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidGenAlphaAdapter::IntegrateInterfaceShape()
{
  return interface_.ExtractCondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidGenAlphaAdapter::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
  interface_.InsertCondVector(ivel,relax);
  fluid_.LinearRelaxationSolve(relax);
  return ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidGenAlphaAdapter::CreateFieldTest()
{
  return Teuchos::rcp(new FluidResultTest(fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidGenAlphaAdapter::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres) const
{
   return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres);
}


#endif  // #ifdef CCADISCRET
