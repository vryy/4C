/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_genalpha.cpp

\brief Fluid field adapter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_fluid_genalpha.H"

ADAPTER::FluidGenAlpha::FluidGenAlpha(
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
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidGenAlpha::InitialGuess()
{
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidGenAlpha::RHS()
{
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidGenAlpha::Velnp()
{
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidGenAlpha::Veln()
{
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidGenAlpha::Dispnp()
{
  return fluid_.Dispnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidGenAlpha::DofRowMap()
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::FluidGenAlpha::SystemMatrix()
{
  return fluid_.SysMat();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidGenAlpha::Discretization()
{
  return fluid_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImpl::StructCondRHS() const
// {
//   return interface_.ExtractCondVector(Velnp());
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlpha::PrepareTimeStep()
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
void ADAPTER::FluidGenAlpha::Evaluate(Teuchos::RCP<const Epetra_Vector> dacc)
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
void ADAPTER::FluidGenAlpha::Update()
{
  fluid_.GenAlphaTimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlpha::Output()
{
  fluid_.GenAlphaOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlpha::NonlinearSolve()
{
  fluid_.DoGenAlphaPredictorCorrectorIteration();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidGenAlpha::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidGenAlpha::PressureRowMap()
{
  return fluid_.PressureRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlpha::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidGenAlpha::ResidualScaling() const
{
  return fluid_.ResidualScaling();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidGenAlpha::TimeScaling() const
{
  double dt = fluid_.Dt();
  double gamma = fluid_.Gamma();
  return 1./(dt*dt*gamma);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlpha::ReadRestart(int step)
{
  fluid_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidGenAlpha::Time()
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::FluidGenAlpha::Step()
{
  return fluid_.Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlpha::LiftDrag()
{
  fluid_.LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidGenAlpha::ExtractInterfaceForces()
{
  return interface_.ExtractCondVector(fluid_.TrueResidual());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidGenAlpha::ExtractInterfaceFluidVelocity()
{
  return interface_.ExtractCondVector(fluid_.Velnp());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlpha::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
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
void ADAPTER::FluidGenAlpha::ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<Epetra_Vector> iforce)
{
  dserror ("You must not use robin-BC with FluidGenAlphaIntegration! It has not been implementet yet!");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlpha::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp)
{
  meshmap_.InsertCondVector(fluiddisp,fluid_.Dispnp());

  // new grid velocity
  fluid_.UpdateGridv();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlpha::ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel)
{
  meshmap_.InsertCondVector(gridvel,fluid_.GridVel());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlpha::ConvertInterfaceUnknown(Teuchos::RCP<Epetra_Vector> fcx)
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
int ADAPTER::FluidGenAlpha::Itemax() const
{
  return fluid_.Itemax();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidGenAlpha::SetItemax(int itemax)
{
  fluid_.SetItemax(itemax);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidGenAlpha::IntegrateInterfaceShape()
{
  return interface_.ExtractCondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidGenAlpha::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
  interface_.InsertCondVector(ivel,relax);
  fluid_.LinearRelaxationSolve(relax);
  return ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidGenAlpha::CreateFieldTest()
{
  return Teuchos::rcp(new FluidResultTest(fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidGenAlpha::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres)
{
   return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres);
}


#endif  // #ifdef CCADISCRET
