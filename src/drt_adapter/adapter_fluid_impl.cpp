/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_impl.cpp

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

#include "adapter_fluid_impl.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidImpl::FluidImpl(
        Teuchos::RCP<DRT::Discretization> dis,
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

  const int numvelids = velmap->NumMyElements();
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
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidImpl::InitialGuess() const
{
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidImpl::RHS() const
{
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidImpl::Velnp() const
{
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidImpl::Veln() const
{
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidImpl::Dispnp() const
{
  return fluid_.Dispnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidImpl::DofRowMap() const
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::FluidImpl::SystemMatrix() const
{
  return fluid_.SystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidImpl::Discretization()
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
void ADAPTER::FluidImpl::PrepareTimeStep()
{
  fluid_.PrepareTimeStep();

  // we add the whole fluid mesh displacement later on?
  //fluid_.Dispnp()->PutScalar(0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImpl::Evaluate(Teuchos::RCP<const Epetra_Vector> vel) const
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
void ADAPTER::FluidImpl::Update()
{
  fluid_.TimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImpl::Output()
{
  fluid_.Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImpl::NonlinearSolve()
{
  fluid_.NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidImpl::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidImpl::PressureRowMap()
{
  return fluid_.PressureRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidImpl::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidImpl::ResidualScaling() const
{
  return fluid_.ResidualScaling();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidImpl::TimeScaling() const
{
  return 1./fluid_.Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImpl::ReadRestart(int step)
{
  fluid_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidImpl::Time()
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::FluidImpl::Step()
{
  return fluid_.Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImpl::LiftDrag()
{
  fluid_.LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImpl::ExtractInterfaceForces()
{
  return interface_.ExtractCondVector(fluid_.TrueResidual());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImpl::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
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
void ADAPTER::FluidImpl::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp) const
{
  meshmap_.InsertCondVector(fluiddisp,fluid_.Dispnp());

  // new grid velocity
  fluid_.UpdateGridv();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImpl::ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel) const
{
  meshmap_.InsertCondVector(gridvel,fluid_.GridVel());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImpl::ConvertInterfaceUnknown(Teuchos::RCP<Epetra_Vector> fcx) const
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
int ADAPTER::FluidImpl::Itemax() const
{
  return fluid_.Itemax();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImpl::SetItemax(int itemax)
{
  fluid_.SetItemax(itemax);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImpl::IntegrateInterfaceShape()
{
  return interface_.ExtractCondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImpl::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
  interface_.InsertCondVector(ivel,relax);
  fluid_.LinearRelaxationSolve(relax);
  return ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidImpl::CreateFieldTest()
{
  return Teuchos::rcp(new FluidResultTest(fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidImpl::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres) const
{
   return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres);
}


#endif  // #ifdef CCADISCRET
