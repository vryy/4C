/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_strugenalpha.cpp

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
#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>



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
  DRT::UTILS::SetupNDimExtractor(*dis,"FSICoupling",interface_);
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
#ifdef INVERSEDESIGNCREATE
  return Teuchos::rcp(new Epetra_Vector(*dis_->DofRowMap()));
#else
  double alphaf = structure_.AlphaF();
  Teuchos::RCP<Epetra_Vector> dispnp = Teuchos::rcp(new Epetra_Vector(*Dispn()));
  dispnp->Update(1./(1.-alphaf),*Dispnm(),-alphaf/(1.-alphaf));
  return dispnp;
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::Dispn()
{
#ifdef INVERSEDESIGNCREATE
  return Teuchos::rcp(new Epetra_Vector(*dis_->DofRowMap()));
#else
  return structure_.Disp();
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::Dispnm()
{
#ifdef INVERSEDESIGNCREATE
  return Teuchos::rcp(new Epetra_Vector(*dis_->DofRowMap()));
#else
  return structure_.Dispm();
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::StructureGenAlpha::DofRowMap()
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
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
  structure_.Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::Output()
{
  structure_.Output();
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
  std::string equil = params_->get<string>("equilibrium iteration","undefined solution algorithm");

  if (structure_.HaveConstraint())
  {
    structure_.FullNewtonLinearUzawa();
  }
  else if (equil=="full newton")
  {
    structure_.FullNewton();
  }
  else if (equil=="line search newton")
  {
    structure_.LineSearchNewton();
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
#ifdef INVERSEDESIGNCREATE
  return Teuchos::rcp(new Epetra_Vector(*interface_.CondMap()));
#else
  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_.Disp());
  return idis;
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::ExtractInterfaceDispnp()
{
#ifdef INVERSEDESIGNCREATE
  return Teuchos::rcp(new Epetra_Vector(*interface_.CondMap()));
#else
  Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractCondVector(structure_.Dispm());
  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_.Disp());

  double alphaf = params_->get<double>("alpha f", 0.459);
  idis->Update(1./(1.-alphaf),*idism,-alphaf/(1.-alphaf));

  return idis;
#endif
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

  Teuchos::RCP<Epetra_Vector> idis;

  switch (Teuchos::getIntegralValue<int>(fsidyn,"PREDICTOR"))
  {
  case 1:
  {
    // d(n)
    // respect Dirichlet conditions at the interface (required for pseudo-rigid body)
    idis  = interface_.ExtractCondVector(structure_.Dispn());
    break;
  }
  case 2:
    // d(n)+dt*(1.5*v(n)-0.5*v(n-1))
    dserror("interface velocity v(n-1) not available");
    break;
  case 3:
  {
    // d(n)+dt*v(n)
    double dt            = params_->get<double>("delta time"             ,0.01);

    idis  = interface_.ExtractCondVector(structure_.Disp());
    Teuchos::RCP<Epetra_Vector> ivel  = interface_.ExtractCondVector(structure_.Vel());

    idis->Update(dt,*ivel,1.0);
    break;
  }
  case 4:
  {
    // d(n)+dt*v(n)+0.5*dt^2*a(n)
    double dt            = params_->get<double>("delta time"             ,0.01);

    idis  = interface_.ExtractCondVector(structure_.Disp());
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


#endif
