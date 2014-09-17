/*----------------------------------------------------------------------*/
/*!
 \file poro_monolithicsplit_nopenetration.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"

#include "../drt_structure/stru_aux.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../linalg/linalg_blocksparsematrix.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_assemblestrategy.H"

#include "poro_monolithicsplit_nopenetration.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROELAST::MonolithicSplitNoPenetration::MonolithicSplitNoPenetration(const Epetra_Comm& comm,
                                                              const Teuchos::ParameterList& timeparams)
  : MonolithicSplit(comm, timeparams)
{

  // Recovering of Lagrange multiplier happens on fluid field
  lambda_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap()));

  return;
}

/*----------------------------------------------------------------------*
 | setup system (called in porolast.cpp)                                 |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::SetupSystem()
{
  //use full maps of both fields. Only Lagrange multipliers are condensed
  {
    // create combined map
    std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;

    vecSpaces.push_back(StructureField()->DofRowMap());
    vecSpaces.push_back(FluidField()->DofRowMap());

    if (vecSpaces[0]->NumGlobalElements() == 0)
      dserror("No structure equation. Panic.");
    if (vecSpaces[1]->NumGlobalElements()==0)
      dserror("No fluid equation. Panic.");

    // full Poroelasticity-map
    fullmap_ = LINALG::MultiMapExtractor::MergeMaps(vecSpaces);
    // full Poroelasticity-blockmap
    blockrowdofmap_->Setup(*fullmap_, vecSpaces);
  }

  // Switch fluid to interface split block matrix
  FluidField()->UseBlockMatrix(true);

  //setup couling objects, system and coupling matrixes
  SetupCouplingAndMatrixes();

  // initialize no penetration coupling matrixes
  k_struct_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(StructureField()->Interface()->FSICondMap()), 81, true, true));

  k_fluid_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(StructureField()->Interface()->FSICondMap()), 81, true, true));

  k_lambda_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(StructureField()->Interface()->FSICondMap()), 81, true, true));

  k_D_ = Teuchos::rcp(new LINALG::SparseMatrix(
          *(StructureField()->Interface()->FSICondMap()), 81, true, true));

  nopenetration_rhs_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap(), true));

  // build map of dofs subjected to a DBC of whole problem
  BuildCombinedDBCMap();
} // SetupSystem()

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::SetupRHS( bool firstcall)
{
  // only Lagrange multipliers are condensed -> use unchanged maps from single fields
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicSplitNoPenetration::SetupRHS");

  // create full monolithic rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  SetupVector(*rhs_,
              StructureField()->RHS(),
              FluidField()->RHS());
}

/*----------------------------------------------------------------------*
 | setup vector of the structure and fluid field            vuong 01/12|
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::SetupVector(Epetra_Vector &f,
                                        Teuchos::RCP<const Epetra_Vector> sv,
                                        Teuchos::RCP<const Epetra_Vector> fv)
{
  // extract dofs of the two fields
  // and put the structural/fluid field vector into the global vector f
  // noticing the block number

  Extractor()->InsertVector(*sv, 0, f);
  Extractor()->InsertVector(*fv, 1, f);
}

/*----------------------------------------------------------------------*
 | setup vector of the structure and fluid field            vuong 01/12|
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::RecoverLagrangeMultiplier()
{
  //Todo: do something!
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicSplitNoPenetration::SetupSystemMatrix");

  //Teuchos::RCP<LINALG::BlockSparseMatrixBase> s = StructureField()->BlockSystemMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField()->SystemMatrix();
  if (s==Teuchos::null)
    dserror("expect structure matrix");
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> f = FluidField()->BlockSystemMatrix();
  if (f==Teuchos::null)
    dserror("expect fluid block matrix");

  //just to play it save ...
  mat.Zero();

  /*----------------------------------------------------------------------*/

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  /*----------------------------------------------------------------------*/
  // structural part k_sf (3nxn)
  // build mechanical-fluid block

  // create empty matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> k_sf = StructFluidCouplingBlockMatrix();

  // call the element and calculate the matrix block
  ApplyStrCouplMatrix(k_sf);

  /*----------------------------------------------------------------------*/
  // fluid part k_fs ( (3n+1)x3n )
  // build fluid-mechanical block

  // create empty matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> k_fs = FluidStructCouplingBlockMatrix();

  // call the element and calculate the matrix block
  ApplyFluidCouplMatrix(k_fs);

  /*----------------------------------------------------------------------*/

  k_fs->Complete();
  k_sf->Complete();

  /*----------------------------------------------------------------------*/
  // pure structural part
  mat.Assign(0,0,View,*s);

  // structure coupling part
  mat.Matrix(0,1).Add(k_sf->Matrix(0,0),false,1.0,0.0);
  mat.Matrix(0,1).Add(k_sf->Matrix(1,0),false,1.0,1.0);

  /*----------------------------------------------------------------------*/
  // pure fluid part
  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  //f->UnComplete();

  //Todo: fill the right blocks ...
  mat.Matrix(1,1).Add(f->Matrix(0,0),false,1.,0.0);
  Teuchos::RCP<LINALG::SparseMatrix> eye = LINALG::Eye(*FluidField()->Interface()->FSICondMap());
  mat.Matrix(1,1).Add(*eye,false,1.,1.0);


  // fluid coupling part
  mat.Matrix(1,0).Add(k_fs->Matrix(0,0),false,1.0,0.0);
  mat.Matrix(1,0).Add(k_fs->Matrix(0,1),false,1.0,1.0);
  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();

  return;
}

/*----------------------------------------------------------------------*
 |    evaluate fluid-structural system matrix at state                |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::ApplyFluidCouplMatrix(
    Teuchos::RCP< LINALG::SparseOperator> k_fs //!< off-diagonal tangent matrix term
  )
{
  //reset
  k_fluid_->Zero();
  k_D_->Zero();
  k_struct_->Zero();
  k_lambda_->Zero();

  //call base class
  Monolithic::ApplyFluidCouplMatrix(k_fs);

  //fill diagonal blocks
  {
    // create the parameters for the discretization
    Teuchos::ParameterList params;
    // action for elements
    params.set<int>("action", FLD::poro_splitnopenetration);
    params.set("total time", Time());
    params.set("delta time", Dt());
    params.set("timescale",FluidField()->ResidualScaling());
    params.set<int>("physical type", FluidField()->PhysicalType());

    FluidField()->Discretization()->ClearState();
    FluidField()->Discretization()->SetState(0,"dispnp",FluidField()->Dispnp());
    FluidField()->Discretization()->SetState(0,"gridv",FluidField()->GridVel());
    FluidField()->Discretization()->SetState(0,"velnp",FluidField()->Velnp());
    FluidField()->Discretization()->SetState(0,"scaaf",FluidField()->Scaaf());

//    FluidField()->Discretization()->SetState(0,"lambda",
//        FluidField()->Interface()->InsertFSICondVector(StructureToFluidAtInterface(lambda_)));

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of FluidField:
    // fluiddofset = 0, structdofset = 1
    DRT::AssembleStrategy fluidstrategy(
        1,              // structdofset for row
        0,              // fluiddofset for column
        k_fluid_,
        k_D_,
        nopenetration_rhs_,
        Teuchos::null,
        Teuchos::null
    );
    FluidField()->Discretization()->EvaluateCondition( params, fluidstrategy,"FSICoupling" );

    FluidField()->Discretization()->ClearState();
  }

  //fill off diagonal blocks
  {
    // create the parameters for the discretization
    Teuchos::ParameterList params;
    // action for elements
    params.set<int>("action", FLD::poro_splitnopenetration_OD);
    params.set("total time", Time());
    params.set("delta time", Dt());
    params.set("timescale",FluidField()->ResidualScaling());
    params.set<int>("physical type", FluidField()->PhysicalType());

    FluidField()->Discretization()->ClearState();
    FluidField()->Discretization()->SetState(0,"dispnp",FluidField()->Dispnp());
    FluidField()->Discretization()->SetState(0,"gridv",FluidField()->GridVel());
    FluidField()->Discretization()->SetState(0,"velnp",FluidField()->Velnp());
    FluidField()->Discretization()->SetState(0,"scaaf",FluidField()->Scaaf());

    FluidField()->Discretization()->SetState(0,"lambda",
        FluidField()->Interface()->InsertFSICondVector(StructureToFluidAtInterface(lambda_)));

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of FluidField:
    // fluiddofset = 0, structdofset = 1
    DRT::AssembleStrategy fluidstrategy(
        1,              // structdofset for row
        1,              // structdofset for column
        k_struct_,           // fluid-mechanical matrix
        k_lambda_,
        Teuchos::null ,
        Teuchos::null,
        Teuchos::null
    );
    FluidField()->Discretization()->EvaluateCondition( params, fluidstrategy,"FSICoupling" );

    FluidField()->Discretization()->ClearState();
  }

  //Todo: transform rows/cols  and put in system matrix
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate mechanical-fluid system matrix at state                    |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplitNoPenetration::ApplyStrCouplMatrix(
    Teuchos::RCP<LINALG::SparseOperator> k_sf //!< off-diagonal tangent matrix term
)
{
  //call base class
  Monolithic::ApplyStrCouplMatrix(k_sf);

  //done.
  return;
}
