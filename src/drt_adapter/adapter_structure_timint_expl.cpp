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
#include "adapter_structure_timint_expl.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*======================================================================*/
/* constructor */
ADAPTER::StructureTimIntExpl::StructureTimIntExpl(
  Teuchos::RCP<STR::TimIntExpl> stie,
  Teuchos::RCP<Teuchos::ParameterList> ioparams,
  Teuchos::RCP<Teuchos::ParameterList> sdynparams,
  Teuchos::RCP<Teuchos::ParameterList> xparams,
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<LINALG::Solver> contactsolver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: structure_(stie),
  discret_(discret),
  ioparams_(ioparams),
  sdynparams_(sdynparams),
  xparams_(xparams),
  solver_(solver),
  contactsolver_(contactsolver),
  output_(output)
{
  // make sure
  if (structure_ == Teuchos::null)
    dserror("Failed to create structural integrator");

  // set-up FSI interface
  //interface_.Setup(*discret_, *discret_->DofRowMap());
  //structure_->SetSurfaceFSI(&interface_);
}


/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimIntExpl::InitialGuess()
{
  dserror("not implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* right-hand side alias the dynamic force residual */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimIntExpl::RHS()
{
  dserror("not implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* get current displacements D_{n+1} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimIntExpl::Dispnp()
{
  return structure_->DisNew();
}


/*----------------------------------------------------------------------*/
/* get last converged displacements D_{n} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimIntExpl::Dispn()
{
  return structure_->Dis();
}


/*----------------------------------------------------------------------*/
/* non-overlapping DOF map */
Teuchos::RCP<const Epetra_Map> ADAPTER::StructureTimIntExpl::DofRowMap()
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}

/*----------------------------------------------------------------------*/
/* non-overlapping DOF map */
Teuchos::RCP<const Epetra_Map> ADAPTER::StructureTimIntExpl::DofRowMap(unsigned nds)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap(nds);
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}

/*----------------------------------------------------------------------*/
/* stiffness, i.e. force residual R_{n+1} differentiated
 * by displacements D_{n+1} */
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::StructureTimIntExpl::SystemMatrix()
{
  dserror("not implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* stiffness, i.e. force residual R_{n+1} differentiated
 * by displacements D_{n+1} */
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::StructureTimIntExpl::BlockSystemMatrix()
{
  dserror("not implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntExpl::UseBlockMatrix()
{
  dserror("not implemented");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntExpl::TSIMatrix()
{
  dserror("not implemented");
}

/*----------------------------------------------------------------------*/
/* get contact manager */
Teuchos::RCP<MORTAR::ManagerBase> ADAPTER::StructureTimIntExpl::ContactManager()
{
  return structure_->ContactManager();
}

/*----------------------------------------------------------------------*/
/* get locsys manager */
Teuchos::RCP<DRT::UTILS::LocsysManager> ADAPTER::StructureTimIntExpl::LocsysManager()
{
  return structure_->LocsysManager();
}

/*----------------------------------------------------------------------*/
/* get discretisation */
Teuchos::RCP<DRT::Discretization> ADAPTER::StructureTimIntExpl::Discretization()
{
  return structure_->Discretization();
}


/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimIntExpl::FRobin()
{
  dserror("not implemented");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/* External force F_{ext,n+1} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimIntExpl::FExtn()
{
  dserror("not implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::FluidCondRHS() const
// {
//   // structure part of the rhs to enforce
//   // u(n+1) dt = d(n+1) - d(n)

//   // extrapolate d(n+1) at the interface and substract d(n)

//   Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractFSICondVector(structure_->Dispm());
//   Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractFSICondVector(structure_->Disp ());

//   double alphaf = structure_->AlphaF();
//   idis->Update(1./(1.-alphaf), *idism, -alphaf/(1.-alphaf)-1.);
//   return idis;
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::MeshCondRHS() const
// {
//   // structure part of the rhs to enforce
//   // d(G,n+1) = d(n+1)

//   // extrapolate d(n+1) at the interface

//   Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractFSICondVector(structure_->Dispm());
//   Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractFSICondVector(structure_->Disp ());

//   double alphaf = structure_->AlphaF();
//   idis->Update(1./(1.-alphaf), *idism, -alphaf/(1.-alphaf));
//   return idis;
// }


/*----------------------------------------------------------------------*/
/* prepare time step */
void ADAPTER::StructureTimIntExpl::PrepareTimeStep()
{
  dserror("not implemented");
}


/*----------------------------------------------------------------------*/
/* build linear system stiffness matrix and rhs/force residual
 *
 * Monolithic FSI accesses the linearised structure problem. */
void ADAPTER::StructureTimIntExpl::Evaluate(
  Teuchos::RCP<const Epetra_Vector>
)
{
  dserror("not implemented");
}


/*----------------------------------------------------------------------*/
/* update time step */
void ADAPTER::StructureTimIntExpl::Update()
{
  structure_->UpdateStepState();
  structure_->UpdateStepTime();
  structure_->UpdateStepElement();
  return;
}


/*----------------------------------------------------------------------*/
/* output */
void ADAPTER::StructureTimIntExpl::Output()
{
  structure_->OutputStep();
}


/*----------------------------------------------------------------------*/
/* domain map */
const Epetra_Map& ADAPTER::StructureTimIntExpl::DomainMap()
{
  return structure_->GetDomainMap();
}


/*----------------------------------------------------------------------*/
/* read restart */
void ADAPTER::StructureTimIntExpl::ReadRestart(int step)
{
  structure_->ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/* find iteratively solution */
void ADAPTER::StructureTimIntExpl::Solve()
{
  dserror("not implemented");
}


/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::RelaxationSolve(
  Teuchos::RCP<Epetra_Vector> iforce
)
{
  dserror("not implemented");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/* extract interface displacements D_{n} */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::ExtractInterfaceDispn()
{
  dserror("not implemented");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/* extract interface displacements D_{n+1} */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::ExtractInterfaceDispnp()
{
  dserror("not implemented");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/* extract external forces at interface F_{ext,n+1} */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::ExtractInterfaceForces()
{
  dserror("not implemented");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::PredictInterfaceDispnp()
{
  dserror("not implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* */
void ADAPTER::StructureTimIntExpl::ApplyInterfaceForces(
  Teuchos::RCP<Epetra_Vector> iforce
)
{
  dserror("not implemented");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntExpl::ApplyInterfaceRobinValue(
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> ifluidvel
)
{
  dserror("Not impl.");
}


/*----------------------------------------------------------------------*/
/* structural result test */
Teuchos::RCP<DRT::ResultTest> ADAPTER::StructureTimIntExpl::CreateFieldTest()
{
  return Teuchos::rcp(new StruResultTest(*structure_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntExpl::Integrate()
{
  structure_->Integrate();
}


/*----------------------------------------------------------------------*/
/* extract displacements D_{n} needed for TSI                dano 05/10 */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::ExtractDispnp()
{
  return structure_->DisNew();
}


/*----------------------------------------------------------------------*/
/* extract displacements D_{n+1} needed for TSI               dano 05/10 */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::ExtractDispn()
{
  return structure_->Dis();
}


/*----------------------------------------------------------------------*/
/* extract velocities V_{n} needed for TSI                   dano 08/10 */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::ExtractVeln()
{
  return structure_->Vel();
}


/*----------------------------------------------------------------------*/
/* extract velocities V_{n+1} needed for TSI                 dano 06/10 */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::ExtractVelnp()
{
  return structure_->VelNew();
}


/*----------------------------------------------------------------------*
 | Extract midpoint velocities                                          |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::ExtractVelaf()
{
  return structure_->Velaf();
}


/*----------------------------------------------------------------------*/
/* apply current temperatures (FSI like)                     dano 03/10 */
void ADAPTER::StructureTimIntExpl::ApplyTemperatures(
  Teuchos::RCP<const Epetra_Vector> temp
  )
{
  dserror("not implemented");
}

/*----------------------------------------------------------------------*
 | Current material displacements (structure with ale)       mgit 05/11 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntExpl::DispMat()
{
  dserror("not implemented");
  return null;
}

/*----------------------------------------------------------------------*
 | Apply material displacements to                           mgit 05/11 |
 | structure field (structure with ale)                                 |
 *----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntExpl::ApplyDisMat(
  Teuchos::RCP<Epetra_Vector> dismat
  )
{
  dserror("not implemented");
}

/*----------------------------------------------------------------------*/
/* prepare partition step                                    dano 12/10 */
/* (iterative staggered partitioned schemes)                            */
void ADAPTER::StructureTimIntExpl::PreparePartitionStep()
{
  dserror("not implemented");
}


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
