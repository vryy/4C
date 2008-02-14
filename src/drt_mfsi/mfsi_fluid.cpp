#ifdef CCADISCRET

#include "mfsi_fluid.H"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MFSI::Fluid::~Fluid()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MFSI::FluidAdapter::FluidAdapter(Teuchos::RCP<DRT::Discretization> dis,
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::FluidAdapter::InitialGuess()
{
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::FluidAdapter::RHS()
{
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::FluidAdapter::Velnp()
{
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::FluidAdapter::Veln()
{
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> MFSI::FluidAdapter::DofRowMap()
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> MFSI::FluidAdapter::SysMat()
{
  return fluid_.SysMat();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> MFSI::FluidAdapter::Discretization()
{
  return fluid_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::FluidAdapter::StructCondRHS()
{
  return interface_.ExtractCondVector(Velnp());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::FluidAdapter::PrepareTimeStep()
{
  fluid_.PrepareTimeStep();

  // we add the whole fluid mesh displacement later on?
  //fluid_.Dispnp()->PutScalar(0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::FluidAdapter::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
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
void MFSI::FluidAdapter::Update()
{
  fluid_.TimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::FluidAdapter::Output()
{
  fluid_.Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::FluidAdapter::SetInterfaceMap(Teuchos::RefCountPtr<Epetra_Map> im)
{
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
Teuchos::RCP<Epetra_Map> MFSI::FluidAdapter::InterfaceMap()
{
  return interface_.CondMap();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> MFSI::FluidAdapter::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> MFSI::FluidAdapter::PressureRowMap()
{
  return fluid_.PressureRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MFSI::FluidAdapter::SetMeshMap(Teuchos::RCP<Epetra_Map> mm)
{
  meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::FluidAdapter::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp)
{
  Teuchos::RCP<Epetra_Vector> deltadispnp = meshmap_.InsertCondVector(fluiddisp);

  //fluid_.Dispnp()->Update(1.0, *deltadispnp, 1.0, *fluid_.Dispn(), 0.0);
  fluid_.Dispnp()->Update(1.0, *deltadispnp, 0.0);

  // new grid velocity
  fluid_.UpdateGridv();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MFSI::FluidAdapter::ResidualScaling()
{
  return fluid_.ResidualScaling();
}


#endif
