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
#include "adapter_structure_timint_adaptive.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*======================================================================*/
/* constructor */
ADAPTER::StructureTimIntAda::StructureTimIntAda(
  Teuchos::RCP<STR::TimAda> sta,
  Teuchos::RCP<STR::TimInt> sti,
  Teuchos::RCP<Teuchos::ParameterList> ioparams,
  Teuchos::RCP<Teuchos::ParameterList> sdynparams,
  Teuchos::RCP<Teuchos::ParameterList> xparams,
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: structure_(sta),
  sti_(sti),
  discret_(discret),
  ioparams_(ioparams),
  sdynparams_(sdynparams),
  xparams_(xparams),
  solver_(solver),
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
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimIntAda::InitialGuess()
{
  dserror("not implemented"); return Teuchos::null;
  //return structure_->Dispm();
}


/*----------------------------------------------------------------------*/
/* right-hand side alias the dynamic force residual */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimIntAda::RHS()
{
  // this expects a _negative_ (Newton-ready) residual with blanked
  // Dirichlet DOFs. We did it in #Evaluate.
  dserror("not implemented"); return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* get current displacements D_{n+1} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimIntAda::Dispnp()
{
  dserror("not implemented"); return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* get last converged displacements D_{n} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimIntAda::Dispn()
{
  dserror("not implemented"); return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* non-overlapping DOF map */
Teuchos::RCP<const Epetra_Map> ADAPTER::StructureTimIntAda::DofRowMap()
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}

/*----------------------------------------------------------------------*/
/* non-overlapping DOF map */
Teuchos::RCP<const Epetra_Map> ADAPTER::StructureTimIntAda::DofRowMap(unsigned nds)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap(nds);
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}

/*----------------------------------------------------------------------*/
/* stiffness, i.e. force residual R_{n+1} differentiated
 * by displacements D_{n+1} */
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::StructureTimIntAda::SystemMatrix()
{
//  cout<<*structure_->SystemMatrix()<<endl;
  dserror("not implemented"); return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* stiffness, i.e. force residual R_{n+1} differentiated
 * by displacements D_{n+1} */
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::StructureTimIntAda::BlockSystemMatrix()
{
  dserror("not implemented"); return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntAda::UseBlockMatrix()
{
  dserror("not implemented");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntAda::TSIMatrix()
{
  dserror("not implemented");
}

/*----------------------------------------------------------------------*/
/* get contact manager */
Teuchos::RCP<MORTAR::ManagerBase> ADAPTER::StructureTimIntAda::ContactManager()
{
  dserror("not implemented"); return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/* get discretisation */
Teuchos::RCP<DRT::Discretization> ADAPTER::StructureTimIntAda::Discretization()
{
  dserror("not implemented"); return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimIntAda::FRobin()
{
  //return structure_->GetForceRobinFSI();
  return LINALG::CreateVector(*discret_->DofRowMap(), true);
}

/*----------------------------------------------------------------------*/
/* External force F_{ext,n+1} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureTimIntAda::FExtn()
{
  dserror("not implemented"); return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntAda::FluidCondRHS() const
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
// Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntAda::MeshCondRHS() const
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
void ADAPTER::StructureTimIntAda::PrepareTimeStep()
{
  // Note: MFSI requires a constant predictor. Otherwise the fields will get
  // out of sync.

  // predict
  dserror("not implemented");
}


/*----------------------------------------------------------------------*/
/* build linear system stiffness matrix and rhs/force residual
 *
 * Monolithic FSI accesses the linearised structure problem. */
void ADAPTER::StructureTimIntAda::Evaluate(
  Teuchos::RCP<const Epetra_Vector> disiterinc
)
{
  dserror("not implemented");
}


/*----------------------------------------------------------------------*/
/* update time step */
void ADAPTER::StructureTimIntAda::Update()
{
  dserror("not implemented");
  return;
}


/*----------------------------------------------------------------------*/
/* output */
void ADAPTER::StructureTimIntAda::Output()
{
  structure_->OutputPeriod();
}


/*----------------------------------------------------------------------*/
/* domain map */
const Epetra_Map& ADAPTER::StructureTimIntAda::DomainMap()
{
  dserror("not implemented");
  return *discret_->DofRowMap(); // make it compile!
}


/*----------------------------------------------------------------------*/
/* read restart */
void ADAPTER::StructureTimIntAda::ReadRestart(int step)
{
  sti_->ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/* find iteratively solution */
void ADAPTER::StructureTimIntAda::Solve()
{
  dserror("not implemented");
}


/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntAda::RelaxationSolve(
  Teuchos::RCP<Epetra_Vector> iforce
)
{
  dserror("not implemented");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/* extract interface displacements D_{n} */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntAda::ExtractInterfaceDispn()
{
  dserror("not implemented");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/* extract interface displacements D_{n+1} */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntAda::ExtractInterfaceDispnp()
{
  dserror("not implemented");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/* extract external forces at interface F_{ext,n+1} */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntAda::ExtractInterfaceForces()
{
  dserror("not implemented");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntAda::PredictInterfaceDispnp()
{
  dserror("not implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* */
void ADAPTER::StructureTimIntAda::ApplyInterfaceForces(
  Teuchos::RCP<Epetra_Vector> iforce
)
{
  dserror("not implemented");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntAda::ApplyInterfaceRobinValue(
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> ifluidvel
)
{
  dserror("Not impl.");
}


/*----------------------------------------------------------------------*/
/* structural result test */
Teuchos::RCP<DRT::ResultTest> ADAPTER::StructureTimIntAda::CreateFieldTest()
{
  return Teuchos::rcp(new StruResultTest(*sti_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimIntAda::Integrate()
{
  structure_->Integrate();
}


/*----------------------------------------------------------------------*/
/* apply the current temperatures (FSI like)                 dano 03/10 */
void ADAPTER::StructureTimIntAda::ApplyTemperatures(
  Teuchos::RCP<Epetra_Vector> itemp
)
{
  dserror("not implemented");
}


/*----------------------------------------------------------------------*/
/* extract displacements needed for coupling in TSI */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntAda::ExtractDispn()
{
  dserror("not implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* extract displacements D_{n+1} needed for coupling in TSI*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureTimIntAda::ExtractDispnp()
{
  dserror("not implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
