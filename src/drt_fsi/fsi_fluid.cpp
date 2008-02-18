#ifdef CCADISCRET

#include "fsi_fluid.H"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Fluid::~Fluid()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidAdapter::FluidAdapter(Teuchos::RCP<DRT::Discretization> dis,
                                 Teuchos::RCP<LINALG::Solver> solver,
                                 Teuchos::RCP<ParameterList> params,
                                 Teuchos::RCP<IO::DiscretizationWriter> output)
  : fluid_(dis, *solver, *params, *output, true),
    dis_(dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  FSI::UTILS::SetupInterfaceExtractor(*dis,"FSICoupling",interface_);

  // build inner velocity map
  // dofs at the interface are excluded
  // we use only velocity dofs and only those without Dirichlet constraint

  Teuchos::RCP<Epetra_Map> velmap = fluid_.VelocityRowMap();
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
Teuchos::RCP<Epetra_Vector> FSI::FluidAdapter::InitialGuess()
{
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAdapter::RHS()
{
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAdapter::Velnp()
{
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAdapter::Veln()
{
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FSI::FluidAdapter::DofRowMap()
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> FSI::FluidAdapter::SysMat()
{
  return fluid_.SysMat();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> FSI::FluidAdapter::Discretization()
{
  return fluid_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAdapter::StructCondRHS()
{
  return interface_.ExtractCondVector(Velnp());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::PrepareTimeStep()
{
  fluid_.PrepareTimeStep();

  // we add the whole fluid mesh displacement later on?
  //fluid_.Dispnp()->PutScalar(0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
{
  // Yes, this is complicated. But we have to be very careful
  // here. The field solver always expects an increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest
  // increment only.
  if (vel!=Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> incvel = Teuchos::rcp(new Epetra_Vector(*vel));
    incvel->Update(-1.0,*fluid_.Velnp(),1.0);
    fluid_.Evaluate(incvel);
  }
  else
  {
    fluid_.Evaluate(Teuchos::null);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::Update()
{
  fluid_.TimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::Output()
{
  fluid_.Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::NonlinearSolve()
{
  fluid_.NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> FSI::FluidAdapter::InterfaceMap()
{
  return interface_.CondMap();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> FSI::FluidAdapter::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> FSI::FluidAdapter::PressureRowMap()
{
  return fluid_.PressureRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::FluidAdapter::SetMeshMap(Teuchos::RCP<Epetra_Map> mm)
{
  meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::FluidAdapter::ResidualScaling()
{
  return fluid_.ResidualScaling();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::ReadRestart(int step)
{
  fluid_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::FluidAdapter::Time()
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int FSI::FluidAdapter::Step()
{
  return fluid_.Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::LiftDrag()
{
  fluid_.LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAdapter::ExtractInterfaceForces()
{
  return interface_.ExtractCondVector(fluid_.TrueResidual());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
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
void FSI::FluidAdapter::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp)
{
  meshmap_.InsertCondVector(fluiddisp,fluid_.Dispnp());

  // new grid velocity
  fluid_.UpdateGridv();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel)
{
  meshmap_.InsertCondVector(gridvel,fluid_.GridVel());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int FSI::FluidAdapter::Itemax() const
{
  return fluid_.Itemax();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::SetItemax(int itemax)
{
  fluid_.SetItemax(itemax);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAdapter::IntegrateInterfaceShape()
{
  return interface_.ExtractCondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAdapter::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
  interface_.InsertCondVector(ivel,relax);
  fluid_.LinearRelaxationSolve(relax);
  return ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> FSI::FluidAdapter::CreateFieldTest()
{
  return Teuchos::rcp(new FluidResultTest(fluid_));
}


#endif
