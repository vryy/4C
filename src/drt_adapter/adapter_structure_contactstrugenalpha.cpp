/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_strugenalpha.cpp

\brief Structure field adapter

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_structure.H"
#include "adapter_structure_contactstrugenalpha.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ContactStructureGenAlpha::ContactStructureGenAlpha(Teuchos::RCP<Teuchos::ParameterList> params,
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
Teuchos::RCP<const Epetra_Vector> ADAPTER::ContactStructureGenAlpha::InitialGuess()
{
#if defined(INVERSEDESIGNCREATE) || defined(PRESTRESS)
  return Teuchos::rcp(new Epetra_Vector(structure_.Getdu().Map(),true));
#else
  return Teuchos::rcp(&structure_.Getdu(),false);
#endif
  //return structure_.Dispm();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ContactStructureGenAlpha::RHS()
{
  return structure_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ContactStructureGenAlpha::Dispnp()
{
#if defined(INVERSEDESIGNCREATE) || defined(PRESTRESS)
  return Teuchos::rcp(new Epetra_Vector(*dis_->DofRowMap(),true));
#else
  double alphaf = structure_.AlphaF();
  Teuchos::RCP<Epetra_Vector> dispnp = Teuchos::rcp(new Epetra_Vector(*Dispn()));
  dispnp->Update(1./(1.-alphaf),*Dispnm(),-alphaf/(1.-alphaf));
  return dispnp;
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ContactStructureGenAlpha::Dispn()
{
#if defined(INVERSEDESIGNCREATE) || defined(PRESTRESS)
  return Teuchos::rcp(new Epetra_Vector(*dis_->DofRowMap(),true));
#else
  return structure_.Disp();
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ContactStructureGenAlpha::Dispnm()
{
#if defined(INVERSEDESIGNCREATE) || defined(PRESTRESS)
  return Teuchos::rcp(new Epetra_Vector(*dis_->DofRowMap(),true));
#else
  return structure_.Dispm();
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::ContactStructureGenAlpha::DofRowMap()
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::ContactStructureGenAlpha::SystemMatrix()
{
  return structure_.SystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::ContactStructureGenAlpha::Discretization()
{
  return structure_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ContactStructureGenAlpha::FRobin()
{
  return structure_.FRobin();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ContactStructureGenAlpha::FExtn()
{
  return structure_.FExtn();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::ContactStructureGenAlpha::FluidCondRHS() const
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
// Teuchos::RCP<Epetra_Vector> ADAPTER::ContactStructureGenAlpha::MeshCondRHS() const
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
void ADAPTER::ContactStructureGenAlpha::PrepareTimeStep()
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ContactStructureGenAlpha::Evaluate(Teuchos::RCP<const Epetra_Vector> disp)
{
  structure_.Evaluate(disp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ContactStructureGenAlpha::Update()
{
  structure_.Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ContactStructureGenAlpha::Output()
{
  structure_.Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Epetra_Map& ADAPTER::ContactStructureGenAlpha::DomainMap()
{
  return structure_.DomainMap();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ContactStructureGenAlpha::ReadRestart(int step)
{
  structure_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ContactStructureGenAlpha::Solve()
{
  std::string equil = params_->get<string>("equilibrium iteration","undefined solution algorithm");

  if (equil=="full newton")
  {
    RCP<CONTACT::Manager> contactmanager = structure_.getContactmanager();
    bool semismooth = (contactmanager->GetStrategy().Params()).get<bool>("semismooth newton",false);

    if(semismooth)
    {
      // initialize convergence status and step no. for active set
      contactmanager->GetStrategy().ActiveSetConverged() = false;
      contactmanager->GetStrategy().ActiveSetSteps() = 1;
      cout << "SEMISMOOTH" << endl;
      structure_.SemiSmoothNewton();
    }
    else
      dserror("only semismooth newton implemented");
  }
  else
    dserror("Unknown type of equilibrium iteration '%s'", equil.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ContactStructureGenAlpha::RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce)
{
  Teuchos::RCP<Epetra_Vector> relax = interface_.InsertCondVector(iforce);
  Teuchos::RCP<Epetra_Vector> idisi = structure_.LinearRelaxationSolve(relax);

  // we are just interested in the incremental interface displacements
  idisi = interface_.ExtractCondVector(idisi);
  return idisi;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ContactStructureGenAlpha::ExtractInterfaceDispn()
{
#if defined(INVERSEDESIGNCREATE) || defined(PRESTRESS)
  return Teuchos::rcp(new Epetra_Vector(*interface_.CondMap(),true));
#else
  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(structure_.Disp());
  return idis;
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ContactStructureGenAlpha::ExtractInterfaceDispnp()
{
#if defined(INVERSEDESIGNCREATE) || defined(PRESTRESS)
  return Teuchos::rcp(new Epetra_Vector(*interface_.CondMap(),true));
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
Teuchos::RCP<Epetra_Vector> ADAPTER::ContactStructureGenAlpha::ExtractInterfaceForces()
{
  Teuchos::RCP<Epetra_Vector> iforce = interface_.ExtractCondVector(structure_.FExtn());

  return iforce;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ContactStructureGenAlpha::PredictInterfaceDispnp()
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  Teuchos::RCP<Epetra_Vector> idis;

  switch (Teuchos::getIntegralValue<int>(fsidyn,"PREDICTOR"))
  {
  case 1:
  {
    // d(n)
    // respect Dirichlet conditions at the interface (required for pseudo-rigid body)
#if defined(INVERSEDESIGNCREATE) || defined(PRESTRESS)
    idis = Teuchos::rcp(new Epetra_Vector(*interface_.CondMap(),true));
#else
    idis  = interface_.ExtractCondVector(structure_.Dispn());
#endif
    break;
  }
  case 2:
    // d(n)+dt*(1.5*v(n)-0.5*v(n-1))
    dserror("interface velocity v(n-1) not available");
    break;
  case 3:
  {
#if defined(INVERSEDESIGNCREATE) || defined(PRESTRESS)
    idis = Teuchos::rcp(new Epetra_Vector(*interface_.CondMap(),true));
#else
    // d(n)+dt*v(n)
    double dt            = params_->get<double>("delta time"             ,0.01);

    idis  = interface_.ExtractCondVector(structure_.Disp());
    Teuchos::RCP<Epetra_Vector> ivel  = interface_.ExtractCondVector(structure_.Vel());

    idis->Update(dt,*ivel,1.0);
#endif
    break;
  }
  case 4:
  {
#if defined(INVERSEDESIGNCREATE) || defined(PRESTRESS)
    idis = Teuchos::rcp(new Epetra_Vector(*interface_.CondMap(),true));
#else
    // d(n)+dt*v(n)+0.5*dt^2*a(n)
    double dt            = params_->get<double>("delta time"             ,0.01);

    idis  = interface_.ExtractCondVector(structure_.Disp());
    Teuchos::RCP<Epetra_Vector> ivel  = interface_.ExtractCondVector(structure_.Vel());
    Teuchos::RCP<Epetra_Vector> iacc  = interface_.ExtractCondVector(structure_.Acc());

    idis->Update(dt,*ivel,0.5*dt*dt,*iacc,1.0);
#endif
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
void ADAPTER::ContactStructureGenAlpha::ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce)
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
void ADAPTER::ContactStructureGenAlpha::ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> iforce,
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
Teuchos::RCP<DRT::ResultTest> ADAPTER::ContactStructureGenAlpha::CreateFieldTest()
{
  return Teuchos::rcp(new StruResultTest(structure_));
}

#endif

