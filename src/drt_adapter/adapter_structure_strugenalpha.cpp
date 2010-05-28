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
#include "../drt_lib/drt_condition_utils.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureGenAlpha::StructureGenAlpha(Teuchos::RCP<Teuchos::ParameterList> params,
                                              Teuchos::RCP<StruGenAlpha> tintegrator,
                                              Teuchos::RCP<DRT::Discretization> dis,
                                              Teuchos::RCP<LINALG::Solver> solver,
                                              Teuchos::RCP<IO::DiscretizationWriter> output)
  : structure_(tintegrator),
    dis_(dis),
    params_(params),
    solver_(solver),
    output_(output)
{
  //setup fsi-Interface
  interface_.Setup(*dis_, *dis_->DofRowMap());
  structure_->SetFSISurface(&interface_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::InitialGuess()
{
  // prestressing business
  double time = 0.0;
  double pstime = -1.0;
  const ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  INPAR::STR::PreStress pstype = Teuchos::getIntegralValue<INPAR::STR::PreStress>(pslist,"PRESTRESS");
  if (pstype != INPAR::STR::prestress_none)
  {
    time   = structure_->GetTime();
    pstime = pslist.get<double>("PRESTRESSTIME");
  }

  if (pstype != INPAR::STR::prestress_none && time <= pstime)
    return Teuchos::rcp(new Epetra_Vector(structure_->Getdu().Map(),true));
  else
    return Teuchos::rcp(&structure_->Getdu(),false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::RHS()
{
  return structure_->Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::Dispnp()
{
  // prestressing business
  double time = 0.0;
  double pstime = -1.0;
  const ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  INPAR::STR::PreStress pstype = Teuchos::getIntegralValue<INPAR::STR::PreStress>(pslist,"PRESTRESS");
  if (pstype != INPAR::STR::prestress_none)
  {
    time = structure_->GetTime();
    pstime = pslist.get<double>("PRESTRESSTIME");
  }

  if (pstype != INPAR::STR::prestress_none && time <= pstime)
    return Teuchos::rcp(new Epetra_Vector(*dis_->DofRowMap(),true));
  else
  {
    double alphaf = structure_->AlphaF();
    Teuchos::RCP<Epetra_Vector> dispnp = Teuchos::rcp(new Epetra_Vector(*Dispn()));
    dispnp->Update(1./(1.-alphaf),*Dispnm(),-alphaf/(1.-alphaf));
    return dispnp;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::Dispn()
{
  // prestressing business
  double time = 0.0;
  double pstime = -1.0;
  const ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  INPAR::STR::PreStress pstype = Teuchos::getIntegralValue<INPAR::STR::PreStress>(pslist,"PRESTRESS");
  if (pstype != INPAR::STR::prestress_none)
  {
    time = structure_->GetTime();
    pstime = pslist.get<double>("PRESTRESSTIME");
  }

  if (pstype != INPAR::STR::prestress_none && time <= pstime)
    return Teuchos::rcp(new Epetra_Vector(*dis_->DofRowMap(),true));
  else
    return structure_->Disp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::Dispnm()
{
  // prestressing business
  double time = 0.0;
  double pstime = -1.0;
  const ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  INPAR::STR::PreStress pstype = Teuchos::getIntegralValue<INPAR::STR::PreStress>(pslist,"PRESTRESS");
  if (pstype != INPAR::STR::prestress_none)
  {
    time = structure_->GetTime();
    pstime = pslist.get<double>("PRESTRESSTIME");
  }

  if (pstype != INPAR::STR::prestress_none && time <= pstime)
    return Teuchos::rcp(new Epetra_Vector(*dis_->DofRowMap(),true));
  else
    return structure_->Dispm();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::StructureGenAlpha::DofRowMap()
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}


/*----------------------------------------------------------------------*/
/* non-overlapping DOF map */
Teuchos::RCP<const Epetra_Map> ADAPTER::StructureGenAlpha::DofRowMap(unsigned nds)
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap(nds);
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::StructureGenAlpha::SystemMatrix()
{
  return structure_->SystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::StructureGenAlpha::BlockSystemMatrix()
{
  return structure_->BlockSystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::UseBlockMatrix()
{
  structure_->UseBlockMatrix(Interface(),Interface());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::TSIMatrix()
{
 // structure_->TSIMatrix();
  dserror("no application here");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MORTAR::ManagerBase> ADAPTER::StructureGenAlpha::ContactManager()
{
  // no contact with tsi in old time integration
  if (structure_->HaveContactMeshtying())
    dserror("TSI and contact only in new time integration");

  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::StructureGenAlpha::Discretization()
{
  return structure_->Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::FRobin()
{
  return structure_->FRobin();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureGenAlpha::FExtn()
{
  return structure_->FExtn();
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
    structure_->ConstantPredictor();
  }
  else if (pred=="consistent")
  {
    structure_->ConsistentPredictor();
  }
  else
    dserror("predictor %s unknown", pred.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  structure_->Evaluate(disiterinc);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::Update()
{
  structure_->Update();
  structure_->UpdateElement();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::Output()
{
  structure_->Output();

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Epetra_Map& ADAPTER::StructureGenAlpha::DomainMap()
{
  return structure_->DomainMap();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::ReadRestart(int step)
{
  structure_->ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::Solve()
{
  std::string equil = params_->get<string>("equilibrium iteration","undefined solution algorithm");

  if (structure_->HaveContactMeshtying())
  {
    structure_->CmtNonlinearSolve();
  }
  else if (structure_->HaveConstraint())
  {
    structure_->FullNewtonLinearUzawa();
  }
  else if (equil=="full newton")
  {
    structure_->FullNewton();
  }
  else if (equil=="line search newton")
  {
    structure_->LineSearchNewton();
  }
  else if (equil=="modified newton")
  {
    structure_->ModifiedNewton();
  }
  else if (equil=="nonlinear cg")
  {
    structure_->NonlinearCG();
  }
  else if (equil=="ptc")
  {
    structure_->PTC();
  }
  else
    dserror("Unknown type of equilibrium iteration '%s'", equil.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce)
{
  Teuchos::RCP<Epetra_Vector> relax = interface_.InsertFSICondVector(iforce);
  Teuchos::RCP<Epetra_Vector> idisi = structure_->LinearRelaxationSolve(relax);

  // we are just interested in the incremental interface displacements
  idisi = interface_.ExtractFSICondVector(idisi);
  return idisi;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::ExtractInterfaceDispn()
{
  // prestressing business
  double time = 0.0;
  double pstime = -1.0;
  const ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  INPAR::STR::PreStress pstype = Teuchos::getIntegralValue<INPAR::STR::PreStress>(pslist,"PRESTRESS");
  if (pstype != INPAR::STR::prestress_none)
  {
    time = structure_->GetTime();
    pstime = pslist.get<double>("PRESTRESSTIME");
  }

  if (pstype != INPAR::STR::prestress_none && time <= pstime)
    return Teuchos::rcp(new Epetra_Vector(*interface_.FSICondMap(),true));
  else
  {
    Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractFSICondVector(structure_->Disp());
    return idis;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::ExtractInterfaceDispnp()
{
  // prestressing business
  double time = 0.0;
  double pstime = -1.0;
  const ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  INPAR::STR::PreStress pstype = Teuchos::getIntegralValue<INPAR::STR::PreStress>(pslist,"PRESTRESS");
  if (pstype != INPAR::STR::prestress_none)
  {
    time = structure_->GetTime();
    pstime = pslist.get<double>("PRESTRESSTIME");
  }

  if (pstype != INPAR::STR::prestress_none && time <= pstime)
    return Teuchos::rcp(new Epetra_Vector(*interface_.FSICondMap(),true));
  else
  {
    Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractFSICondVector(structure_->Dispm());
    Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractFSICondVector(structure_->Disp());
    double alphaf = params_->get<double>("alpha f", 0.459);
    idis->Update(1./(1.-alphaf),*idism,-alphaf/(1.-alphaf));
    return idis;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::ExtractInterfaceForces()
{
  Teuchos::RCP<Epetra_Vector> iforce = interface_.ExtractFSICondVector(structure_->FExtn());

  return iforce;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::PredictInterfaceDispnp()
{
  // prestressing business
  double time = 0.0;
  double pstime = -1.0;
  const ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  INPAR::STR::PreStress pstype = Teuchos::getIntegralValue<INPAR::STR::PreStress>(pslist,"PRESTRESS");
  if (pstype != INPAR::STR::prestress_none)
  {
    time = structure_->GetTime();
    pstime = pslist.get<double>("PRESTRESSTIME");
  }

  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  Teuchos::RCP<Epetra_Vector> idis;

  switch (Teuchos::getIntegralValue<int>(fsidyn,"PREDICTOR"))
  {
  case 1:
  {
    // d(n)
    // respect Dirichlet conditions at the interface (required for pseudo-rigid body)
    if (pstype != INPAR::STR::prestress_none && time <= pstime)
      idis = Teuchos::rcp(new Epetra_Vector(*interface_.FSICondMap(),true));
    else
      idis  = interface_.ExtractFSICondVector(structure_->Dispn());
    break;
  }
  case 2:
    // d(n)+dt*(1.5*v(n)-0.5*v(n-1))
    dserror("interface velocity v(n-1) not available");
    break;
  case 3:
  {
    if (pstype != INPAR::STR::prestress_none && time <= pstime)
      idis = Teuchos::rcp(new Epetra_Vector(*interface_.FSICondMap(),true));
    else
    {
      // d(n)+dt*v(n)
      double dt            = params_->get<double>("delta time"             ,0.01);
      idis  = interface_.ExtractFSICondVector(structure_->Disp());
      Teuchos::RCP<Epetra_Vector> ivel  = interface_.ExtractFSICondVector(structure_->Vel());
      idis->Update(dt,*ivel,1.0);
    }
    break;
  }
  case 4:
  {
    if (pstype != INPAR::STR::prestress_none && time <= pstime)
      idis = Teuchos::rcp(new Epetra_Vector(*interface_.FSICondMap(),true));
    else
    {
      // d(n)+dt*v(n)+0.5*dt^2*a(n)
      double dt            = params_->get<double>("delta time"             ,0.01);
      idis  = interface_.ExtractFSICondVector(structure_->Disp());
      Teuchos::RCP<Epetra_Vector> ivel  = interface_.ExtractFSICondVector(structure_->Vel());
      Teuchos::RCP<Epetra_Vector> iacc  = interface_.ExtractFSICondVector(structure_->Acc());
      idis->Update(dt,*ivel,0.5*dt*dt,*iacc,1.0);
    }
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
  structure_->ApplyExternalForce(interface_,iforce);
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

  Teuchos::RCP<Epetra_Vector> idisn  = interface_.ExtractFSICondVector(structure_->Disp());
  Teuchos::RCP<Epetra_Vector> frobin = interface_.ExtractFSICondVector(structure_->FRobin());

  // save robin coupling values in frobin vector (except iforce which
  // is passed separately)
  frobin->Update(alphas/dt,*idisn,alphas*(1-alphaf),*ifluidvel,0.0);

  interface_.InsertFSICondVector(frobin,structure_->FRobin());
  structure_->ApplyExternalForce(interface_,iforce);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::StructureGenAlpha::CreateFieldTest()
{
  return Teuchos::rcp(new StruResultTest(*structure_));
}


/*----------------------------------------------------------------------*
 | Apply current temperature  (for TSI)                      dano 03/10 |
 *----------------------------------------------------------------------*/
void ADAPTER::StructureGenAlpha::ApplyTemperatures(
  Teuchos::RCP<Epetra_Vector> itemp
  )
{
  // Play it save. In the first iteration everything is already set up
  // properly. However, all following iterations need to calculate the
  // stiffness matrix here. Furthermore we are bound to reset fextm_
  // before we add our special contribution.
  // So we calculate the stiffness anyway (and waste the available
  // stiffness in the first iteration).
  dserror("no application here");
}


/*----------------------------------------------------------------------*
 | Extract displacements needed for TSI                      dano 05/10 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::ExtractDispn()
{
  dserror("no application here");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 | Extract current displacements needed for TSI              dano 05/10 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureGenAlpha::ExtractDispnp()
{
  dserror("no application here");
  return Teuchos::null;
}


#endif
