/*----------------------------------------------------------------------*/
/*!
\file ad_scatra_wrapper_cellmigration.cpp
\brief Cell Migration specific wrapper for the scatra time integrator.
\level 1
\maintainer Andreas Rauch
*/
/*----------------------------------------------------------------------*/


#include "ad_scatra_wrapper_cellmigration.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../drt_mat/biochemo_mechano_cell_activefiber.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
ADAPTER::AdapterScatraWrapperCellMigration::AdapterScatraWrapperCellMigration(Teuchos::RCP<ScatraInterface> scatra)
:   AdapterScatraWrapper(scatra)
{
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  rates_ = LINALG::CreateVector(*dofrowmap,true);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::AdapterScatraWrapperCellMigration::EvaluateAdditionalSolutionDependingModels(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,      //!< system matrix
    Teuchos::RCP<Epetra_Vector>          rhs
)
{
  // get rcp to cell discretization
  Teuchos::RCP<DRT::Discretization> celldis = DRT::Problem::Instance()->GetDis("cell");

  // Evaluate biochemical transport and reaction at focal adhesion sites
  // only if BioChemoMechanoCellActiveFiber is used in cell.
  DRT::Condition* FocalAdhesionCondition = Discretization()->GetCondition("CellFocalAdhesion");
  if(FocalAdhesionCondition != NULL and
     Teuchos::rcp_dynamic_cast<MAT::BioChemoMechanoCellActiveFiber>(
         celldis->gElement(celldis->ElementRowMap()->GID(0))->Material() ) != Teuchos::null
     )
  {
    // declare ParameterList
    Teuchos::ParameterList params;

    // zero rates vector
    rates_->Scale(0.0);

    params.set<int>("action",SCATRA::calc_cell_mechanotransduction);
    params.set<int>("ndsdisp",NdsDisp());

    // add element parameters and set state vectors according to time-integration scheme
    AddTimeIntegrationSpecificVectors();
    Discretization()->Evaluate(params,Teuchos::null,Teuchos::null,rates_,Teuchos::null,Teuchos::null);

    /* -------------------------------------------------------------------------------------------*/

    params.set<int>("action",SCATRA::bd_calc_mechanotransduction);

    // add element parameters and set state vectors according to time-integration scheme
    Discretization()->EvaluateCondition(params,Teuchos::null,Teuchos::null,rhs,Teuchos::null,Teuchos::null,"CellFocalAdhesion",-1);

    // safety
    Discretization()->ClearState();

    // write result to right-hand-side
    rhs->Update(1.0,*rates_,1.0);

  } // at focal adhesions

  return;
}
