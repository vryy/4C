#ifdef CCADISCRET

#include "fsi_structure.H"
#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Structure::~Structure()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::StructureAdapter::StructureAdapter(Teuchos::RCP<Teuchos::ParameterList> params,
                                         Teuchos::RCP<DRT::Discretization> dis,
                                         Teuchos::RCP<LINALG::Solver> solver,
                                         Teuchos::RCP<IO::DiscretizationWriter> output)
  : structure_(*params, *dis, *solver, *output),
    dis_(dis),
    params_(params),
    solver_(solver),
    output_(output)
{
  FSI::UTILS::SetupInterfaceExtractor(*dis,"FSICoupling",interface_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::StructureAdapter::InitialGuess() const
{
  return Teuchos::rcp(&structure_.Getdu(),false);
  //return structure_.Dispm();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::StructureAdapter::RHS() const
{
  return structure_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::StructureAdapter::Dispnp() const
{
  double alphaf = structure_.AlphaF();
  Teuchos::RCP<Epetra_Vector> dispnp = Teuchos::rcp(new Epetra_Vector(*Dispn()));
  dispnp->Update(1./(1.-alphaf),*Dispnm(),-alphaf/(1.-alphaf));
  return dispnp;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::StructureAdapter::Dispn() const
{
  return structure_.Disp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::StructureAdapter::Dispnm() const
{
  return structure_.Dispm();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FSI::StructureAdapter::DofRowMap() const
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> FSI::StructureAdapter::SystemMatrix() const
{
  return structure_.SystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> FSI::StructureAdapter::Discretization()
{
  return structure_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::StructureAdapter::DispIncrFactor()
{
  return structure_.DispIncrFactor();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> FSI::StructureAdapter::FluidCondRHS() const
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
// Teuchos::RCP<Epetra_Vector> FSI::StructureAdapter::MeshCondRHS() const
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
void FSI::StructureAdapter::PrepareTimeStep()
{
  // Note: FSI requires a constant predictor. Otherwise the fields will get
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::StructureAdapter::Evaluate(Teuchos::RCP<const Epetra_Vector> disp) const
{
  // Yes, this is complicated. But we have to be very careful
  // here. The field solver always expects an increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest
  // increment only.
  if (disp!=Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(new Epetra_Vector(*disp));
    disi->Update(-1.0,*Dispnp(),1.0);
    structure_.Evaluate(disi);
  }
  else
  {
    structure_.Evaluate(Teuchos::null);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::StructureAdapter::Update()
{
  structure_.UpdateandOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::StructureAdapter::Output()
{
  // noop
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Epetra_Map& FSI::StructureAdapter::DomainMap()
{
  return structure_.DomainMap();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::StructureAdapter::ReadRestart(int step)
{
  structure_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::StructureAdapter::Solve()
{
  std::string equil = params_->get<string>("equilibrium iteration","full newton");
  if (equil=="full newton")
  {
    structure_.FullNewton();
  }
  else if (equil=="modified newton")
  {
    structure_.ModifiedNewton();
  }
  else if (equil=="matrixfree newton")
  {
    structure_.MatrixFreeNewton();
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
Teuchos::RCP<Epetra_Vector> FSI::StructureAdapter::RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce)
{
  Teuchos::RCP<Epetra_Vector> relax = interface_.InsertCondVector(iforce);
  Teuchos::RCP<Epetra_Vector> idisi = structure_.LinearRelaxationSolve(relax);

  // we are just interested in the incremental interface displacements
  idisi = interface_.ExtractCondVector(idisi);
  return idisi;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::StructureAdapter::ExtractInterfaceDispn()
{
  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_.Disp());
  return idis;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::StructureAdapter::ExtractInterfaceDispnp()
{
  Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractCondVector(structure_.Dispm());
  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_.Disp());

  double alphaf = params_->get<double>("alpha f", 0.459);
  idis->Update(1./(1.-alphaf),*idism,-alphaf/(1.-alphaf));

  return idis;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::StructureAdapter::PredictInterfaceDispnp()
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
void FSI::StructureAdapter::ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce)
{
  // Play it save. In the first iteration everything is already set up
  // properly. However, all following iterations need to calculate the
  // stiffness matrix here. Furthermore we are bound to reset fextm_
  // before we add our special contribution.
  // So we calculate the stiffness anyway (and waste the available
  // stiffness in the first iteration).

  structure_.ApplyExternalForce(interface_,iforce);
}


#endif
