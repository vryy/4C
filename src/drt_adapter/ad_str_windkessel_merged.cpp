/*----------------------------------------------------------------------*/
/*!
\file ad_str_windkessel_merged.cpp

\brief Adapter Layer for Structures with 0D Windkessel coupling

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

***************************************************************************************************************
NOT SURE WHETHER WE WILL NEED THIS ADAPTER IF USING THE WK-STRUCTURE COUPLING WITHIN ANOTHER COUPLED PROBLEM!!!
IT EMERGED FROM THE CONSTRAINT FRAMEWORK AND WAS CONSTRUCTED ANALOGOUSLY!!!
***************************************************************************************************************
*/
/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_structure/strtimint_create.H"
#include "ad_str_windkessel_merged.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_constraint/windkessel_manager.H"
#include "../drt_structure/stru_aux.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*======================================================================*/
/* constructor */
ADAPTER::StructureWindkesselMerged::StructureWindkesselMerged
(
  Teuchos::RCP<Structure> stru
)
: StructureWrapper(stru)
{
  // make sure
  if (structure_ == Teuchos::null)
    dserror("Failed to create the underlying structural adapter");

  // build merged dof row map
  dofrowmap_ = LINALG::MergeMap(*(structure_->DofRowMap()),
                                *(structure_->GetWindkesselManager()->GetWindkesselMap()),
                                false);

  // set up interface between merged and single maps
  windkmerger_ = Teuchos::rcp(new LINALG::MapExtractor);
  windkmerger_->Setup(*dofrowmap_,
                    structure_->DofRowMap(),
                    structure_->GetWindkesselManager()->GetWindkesselMap());

  //setup fsi-Interface
  interface_ = Teuchos::rcp(new STR::AUX::MapExtractor);
  interface_->Setup(*Discretization(), *dofrowmap_);
}



/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureWindkesselMerged::InitialGuess()
{
  //get initial guesses from structure and windkesselmanager
  Teuchos::RCP<const Epetra_Vector> strucGuess = structure_->InitialGuess();
  Teuchos::RCP<const Epetra_Vector> wkGuess = Teuchos::rcp(new Epetra_Vector(*(structure_->GetWindkesselManager()->GetWindkesselMap()),true));

  //merge stuff together
  Teuchos::RCP<Epetra_Vector> mergedGuess = Teuchos::rcp(new Epetra_Vector(*dofrowmap_,true));
  windkmerger_->AddCondVector(strucGuess,mergedGuess);
  windkmerger_->AddOtherVector(wkGuess,mergedGuess);

  return mergedGuess;
}

/*----------------------------------------------------------------------*/
/* right-hand side alias the dynamic force residual */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureWindkesselMerged::RHS()
{
  //get rhs-vectors from structure and windkesselmanager
  Teuchos::RCP<const Epetra_Vector> struRHS = structure_->RHS();
  Teuchos::RCP<const Epetra_Vector> wkRHS = structure_->GetWindkesselManager()->GetWindkesselRHS();

  //merge stuff together
  Teuchos::RCP<Epetra_Vector> mergedRHS = Teuchos::rcp(new Epetra_Vector(*dofrowmap_,true));
  windkmerger_->AddCondVector(struRHS,mergedRHS);
  windkmerger_->AddOtherVector(-1.0,wkRHS,mergedRHS);

  return mergedRHS;
}


/*----------------------------------------------------------------------*/
/* get current displacements D_{n+1} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureWindkesselMerged::Dispnp()
{
  //get current state from structure and windkesselmanager
  Teuchos::RCP<const Epetra_Vector> strudis = structure_->Dispnp();
  Teuchos::RCP<const Epetra_Vector> wkdof = structure_->GetWindkesselManager()->GetDofVector();

  //merge stuff together
  Teuchos::RCP<Epetra_Vector> mergedstat = Teuchos::rcp(new Epetra_Vector(*dofrowmap_,true));
  windkmerger_->AddCondVector(strudis,mergedstat);
  windkmerger_->AddOtherVector(wkdof,mergedstat);

  return mergedstat;
}


/*----------------------------------------------------------------------*/
/* get last converged displacements D_{n} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureWindkesselMerged::Dispn()
{
  //get last converged state from structure and windkesselmanager
  Teuchos::RCP<const Epetra_Vector> strudis = structure_->Dispn();
  Teuchos::RCP<const Epetra_Vector> pres = structure_->GetWindkesselManager()->GetDofVectorOld();

  //merge stuff together
  Teuchos::RCP<Epetra_Vector> mergedstat = Teuchos::rcp(new Epetra_Vector(*dofrowmap_,true));
  windkmerger_->AddCondVector(strudis,mergedstat);
  windkmerger_->AddOtherVector(pres,mergedstat);

  return mergedstat;
}

/*----------------------------------------------------------------------*/
/* non-overlapping DOF map */
Teuchos::RCP<const Epetra_Map> ADAPTER::StructureWindkesselMerged::DofRowMap()
{
  return dofrowmap_;
}


/*----------------------------------------------------------------------*/
/* stiffness, i.e. force residual R_{n+1} differentiated
 * by displacements D_{n+1} */
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::StructureWindkesselMerged::SystemMatrix()
{
  //create empty large matrix and get small ones from structure and Windkessel functions
  Teuchos::RCP<LINALG::SparseMatrix> mergedmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_, 81));
  Teuchos::RCP<LINALG::SparseMatrix> strustiff = structure_->SystemMatrix();
  strustiff->Complete();
  
  Teuchos::RCP<LINALG::SparseOperator> mat_dwindk_dd = structure_->GetWindkesselManager()->GetMatDwindkDd();
  mat_dwindk_dd->Complete();

  Teuchos::RCP<LINALG::SparseOperator> mat_dstruct_dwkdof = structure_->GetWindkesselManager()->GetMatDstructDp();
  mat_dstruct_dwkdof->Complete();

  Teuchos::RCP<LINALG::SparseMatrix> windkstiff = structure_->GetWindkesselManager()->GetWindkesselStiffness();
  windkstiff->Complete();

  // Add matrices together
  mergedmatrix -> Add(*strustiff,false,1.0,0.0);
  mergedmatrix -> Add(*mat_dwindk_dd,false,1.0,1.0);
  mergedmatrix -> Add(*mat_dstruct_dwkdof,false,1.0,1.0);
  mergedmatrix -> Add(*windkstiff,false,1.0,1.0);
  mergedmatrix -> Complete(*dofrowmap_,*dofrowmap_);

  mergedmatrix -> ApplyDirichlet( *(structure_->GetDBCMapExtractor()->CondMap()));

  return mergedmatrix;
}


/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::StructureWindkesselMerged::BlockSystemMatrix()
{
  dserror("constrained BlockSparseMatrix never to be implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* build linear system stiffness matrix and rhs/force residual
 *
 * Monolithic FSI accesses the linearised structure problem. */
void ADAPTER::StructureWindkesselMerged::Evaluate(
  Teuchos::RCP<const Epetra_Vector> dispstepinc
)
{
  // 'initialize' structural displacement as null-pointer
  Teuchos::RCP<Epetra_Vector> dispstructstepinc = Teuchos::null;

  // Compute residual increments, update total increments and update pressure increments
  if (dispstepinc != Teuchos::null)
  {
    // Extract increments for pressures and do update
    Teuchos::RCP<Epetra_Vector> presincr = windkmerger_->ExtractOtherVector(dispstepinc);
    structure_->UpdateIterIncrWindkessel(presincr);
    dispstructstepinc = windkmerger_->ExtractCondVector(dispstepinc);
  }
  // Hand down incremental displacements,
  // structure_ will compute the residual increments on its own
  structure_->Evaluate(dispstructstepinc);

}


/*----------------------------------------------------------------------*/
/* domain map */
const Epetra_Map& ADAPTER::StructureWindkesselMerged::DomainMap()
{
  return *(LINALG::MergeMap(structure_->DomainMap(),
                            *(structure_->GetWindkesselManager()->GetWindkesselMap()),
                            false));
}

/*======================================================================*/
void ADAPTER::StructureWindkesselMerged::SetPressure(Teuchos::RCP<Epetra_Vector> couppres)
{
  const Epetra_BlockMap& condmap = couppres->Map();

  for (int i=0; i<condmap.NumMyElements(); ++i)
  {
    int condID = condmap.GID(i);
    DRT::Condition* cond = coupcond_[condID];
    std::vector<double> newval(6,0.0);
    newval[0] = (*couppres)[i];
    cond->Add("val",newval);
  }
}

