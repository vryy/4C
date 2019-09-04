/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problem with constraints

\level 2

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "fsi_constrmonolithic_structuresplit.H"
#include "fsi_matrixtransform.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_io/io_control.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/ad_fld_fluid.H"

#include "../linalg/linalg_utils.H"

#include "../drt_constraint/constraint_manager.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"
#include "../drt_adapter/ad_ale_fsi.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::ConstrMonolithicStructureSplit::ConstrMonolithicStructureSplit(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : ConstrMonolithic(comm, timeparams)
{
  // ---------------------------------------------------------------------------
  // FSI specific check of Dirichlet boundary conditions
  // ---------------------------------------------------------------------------
  // Create intersection of slave DOFs that hold a Dirichlet boundary condition
  // and are located at the FSI interface
  std::vector<Teuchos::RCP<const Epetra_Map>> intersectionmaps;
  intersectionmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(StructureField()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  // Check whether the intersection is empty
  if (intersectionmap->NumGlobalElements() != 0)
  {
    // It is not allowed, that slave DOFs at the interface hold a Dirichlet
    // boundary condition. Thus --> ToDo: Error message

    // We do not have to care whether ALE interface DOFs carry DBCs in the
    // input file since they do not occur in the monolithic system and, hence,
    // do not cause a conflict.

    std::stringstream errormsg;
    errormsg << "  "
                "+---------------------------------------------------------------------------------"
                "------------+"
             << std::endl
             << "  |                DIRICHLET BOUNDARY CONDITIONS ON SLAVE SIDE OF FSI INTERFACE   "
                "              |"
             << std::endl
             << "  "
                "+---------------------------------------------------------------------------------"
                "------------+"
             << std::endl
             << "  | NOTE: The slave side of the interface is not allowed to carry Dirichlet "
                "boundary conditions.|"
             << std::endl
             << "  |                                                                               "
                "              |"
             << std::endl
             << "  | This is a structure split scheme. Hence, master and slave field are chosen as "
                "follows:      |"
             << std::endl
             << "  |     MASTER  = FLUID                                                           "
                "              |"
             << std::endl
             << "  |     SLAVE   = STRUCTURE                                                       "
                "              |"
             << std::endl
             << "  |                                                                               "
                "              |"
             << std::endl
             << "  | Dirichlet boundary conditions were detected on slave interface degrees of "
                "freedom. Please   |"
             << std::endl
             << "  | remove Dirichlet boundary conditions from the slave side of the FSI "
                "interface.              |"
             << std::endl
             << "  | Only the master side of the FSI interface is allowed to carry Dirichlet "
                "boundary conditions.|"
             << std::endl
             << "  "
                "+---------------------------------------------------------------------------------"
                "------------+"
             << std::endl;

    std::cout << errormsg.str();
  }
  // ---------------------------------------------------------------------------

  sggtransform_ = Teuchos::rcp(new UTILS::MatrixRowColTransform);
  sgitransform_ = Teuchos::rcp(new UTILS::MatrixRowTransform);
  sigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  aigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  scgitransform_ = Teuchos::rcp(new UTILS::MatrixRowTransform);
  csigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::SetupSystem()
{
  GeneralSetup();

  // create combined map
  CreateCombinedDofRowMap();

  FluidField()->UseBlockMatrix(false);

  // Use splitted structure matrix
  StructureField()->UseBlockMatrix();

  Teuchos::RCP<Epetra_Map> emptymap =
      Teuchos::rcp(new Epetra_Map(-1, 0, NULL, 0, StructureField()->Discretization()->Comm()));
  Teuchos::RCP<LINALG::MapExtractor> extractor;
  extractor->Setup(*conman_->GetConstraintMap(), emptymap, conman_->GetConstraintMap());
  conman_->UseBlockMatrix(extractor, StructureField()->Interface());
  sconT_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *StructureField()->Interface(), *extractor, 81, false, true));

  // build ale system matrix in splitted system
  AleField()->CreateSystemMatrix(AleField()->Interface());

  CreateSystemMatrix(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::CreateCombinedDofRowMap()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(StructureField()->Interface()->OtherMap());
  vecSpaces.push_back(FluidField()->DofRowMap());
  vecSpaces.push_back(AleField()->Interface()->OtherMap());

  vecSpaces.push_back(conman_->GetConstraintMap());

  if (vecSpaces[0]->NumGlobalElements() == 0)
    dserror("No inner structural equations. Splitting not possible. Panic.");

  SetDofRowMaps(vecSpaces);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::SetupRHSResidual(Epetra_Vector& f)
{
  const double scale = FluidField()->ResidualScaling();

  SetupVector(f, StructureField()->RHS(), FluidField()->RHS(), AleField()->RHS(),
      conman_->GetError(), scale);

  // add additional ale residual
  Extractor().AddVector(*aleresidual_, 2, f);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::SetupRHSLambda(Epetra_Vector& f)
{
  // ToDo: We still need to implement this.

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::SetupRHSFirstiter(Epetra_Vector& f)
{
  // additional rhs term for ALE equations
  // -dt Aig u(n)
  //
  //    1/dt Delta d(n+1) = theta Delta u(n+1) + u(n)
  //
  // And we are concerned with the u(n) part here.

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField()->BlockSystemMatrix();
  if (a == Teuchos::null) dserror("expect ale block matrix");

  const LINALG::SparseMatrix& aig = a->Matrix(0, 1);

  Teuchos::RCP<const Epetra_Vector> fveln = FluidField()->ExtractInterfaceVeln();
  Teuchos::RCP<const Epetra_Vector> sveln = FluidToStruct(fveln);
  Teuchos::RCP<const Epetra_Vector> aveln = StructToAle(sveln);
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap()));
  aig.Apply(*aveln, *rhs);

  rhs->Scale(-1. * Dt());

  Extractor().AddVector(*rhs, 2, f);

  // structure
  Teuchos::RCP<Epetra_Vector> veln = StructureField()->Interface()->InsertFSICondVector(sveln);
  rhs = Teuchos::rcp(new Epetra_Vector(veln->Map()));

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> s = StructureField()->BlockSystemMatrix();
  s->Apply(*veln, *rhs);

  rhs->Scale(-1. * Dt());

  veln = StructureField()->Interface()->ExtractOtherVector(rhs);
  Extractor().AddVector(*veln, 0, f);

  veln = StructureField()->Interface()->ExtractFSICondVector(rhs);
  veln = FluidField()->Interface()->InsertFSICondVector(StructToFluid(veln));

  const double scale = FluidField()->ResidualScaling();

  veln->Scale(1. / scale);

  Extractor().AddVector(*veln, 1, f);

  // shape derivatives
  Teuchos::RCP<const LINALG::BlockSparseMatrixBase> mmm = FluidField()->ShapeDerivatives();
  if (mmm != Teuchos::null)
  {
    const LINALG::SparseMatrix& fmig = mmm->Matrix(0, 1);
    const LINALG::SparseMatrix& fmgg = mmm->Matrix(1, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap()));
    fmig.Apply(*fveln, *rhs);
    veln = FluidField()->Interface()->InsertOtherVector(rhs);

    rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap()));
    fmgg.Apply(*fveln, *rhs);
    FluidField()->Interface()->InsertFSICondVector(rhs, veln);

    veln->Scale(-1. * Dt());

    Extractor().AddVector(*veln, 1, f);
  }

  //--------------------------------------------------------------------------------
  // constraint
  //--------------------------------------------------------------------------------
  LINALG::SparseOperator& tmp = *conman_->GetConstrMatrix();
  LINALG::BlockSparseMatrixBase& scon = dynamic_cast<LINALG::BlockSparseMatrixBase&>(tmp);
  for (int rowblock = 0; rowblock < scon.Rows(); ++rowblock)
  {
    for (int colblock = 0; colblock < scon.Cols(); ++colblock)
    {
      sconT_->Matrix(colblock, rowblock).Add(scon.Matrix(rowblock, colblock), true, 1.0, 0.0);
    }
  }
  sconT_->Complete();

  LINALG::SparseMatrix& csig = sconT_->Matrix(0, 1);

  rhs = Teuchos::rcp(new Epetra_Vector(csig.RowMap()));
  csig.Apply(*sveln, *rhs);
  rhs->Scale(-1. * Dt());
  Extractor().AddVector(*rhs, 3, f);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::SetupSystemMatrix");

  // extract Jacobian matrices and put them into composite system
  // matrix W

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  // const ADAPTER::Coupling& coupsa = StructureAleCoupling();

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> s = StructureField()->BlockSystemMatrix();
  if (s == Teuchos::null) dserror("expect structure block matrix");
  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField()->SystemMatrix();
  if (f == Teuchos::null) dserror("expect fluid matrix");
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField()->BlockSystemMatrix();
  if (a == Teuchos::null) dserror("expect ale block matrix");

  LINALG::SparseMatrix& aii = a->Matrix(0, 0);
  LINALG::SparseMatrix& aig = a->Matrix(0, 1);

  /*----------------------------------------------------------------------*/

  double scale = FluidField()->ResidualScaling();
  double timescale = FluidField()->TimeScaling();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  f->UnComplete();

  mat.Assign(0, 0, LINALG::View, s->Matrix(0, 0));

  (*sigtransform_)(s->FullRowMap(), s->FullColMap(), s->Matrix(0, 1), 1. / timescale,
      ADAPTER::CouplingMasterConverter(coupsf), mat.Matrix(0, 1));
  (*sggtransform_)(s->Matrix(1, 1), 1. / (scale * timescale),
      ADAPTER::CouplingMasterConverter(coupsf), ADAPTER::CouplingMasterConverter(coupsf), *f, true,
      true);
  (*sgitransform_)(
      s->Matrix(1, 0), 1. / scale, ADAPTER::CouplingMasterConverter(coupsf), mat.Matrix(1, 0));

  mat.Assign(1, 1, LINALG::View, *f);

  (*aigtransform_)(a->FullRowMap(), a->FullColMap(), aig, 1. / timescale,
      ADAPTER::CouplingSlaveConverter(*icoupfa_), mat.Matrix(2, 1));
  mat.Assign(2, 2, LINALG::View, aii);

  /*----------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField()->ShapeDerivatives();
  if (mmm != Teuchos::null)
  {
    LINALG::SparseMatrix& fmii = mmm->Matrix(0, 0);
    LINALG::SparseMatrix& fmig = mmm->Matrix(0, 1);
    LINALG::SparseMatrix& fmgi = mmm->Matrix(1, 0);
    LINALG::SparseMatrix& fmgg = mmm->Matrix(1, 1);

    mat.Matrix(1, 1).Add(fmgg, false, 1. / timescale, 1.0);
    mat.Matrix(1, 1).Add(fmig, false, 1. / timescale, 1.0);

    const ADAPTER::Coupling& coupfa = FluidAleCoupling();

    (*fmgitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmgi, 1.,
        ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, false);

    (*fmiitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmii, 1.,
        ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, true);
  }


  /*----------------------------------------------------------------------*/
  // structure constraint part

  LINALG::SparseOperator& tmp = *conman_->GetConstrMatrix();
  LINALG::BlockSparseMatrixBase& scon = dynamic_cast<LINALG::BlockSparseMatrixBase&>(tmp);
  for (int rowblock = 0; rowblock < scon.Rows(); ++rowblock)
  {
    for (int colblock = 0; colblock < scon.Cols(); ++colblock)
    {
      sconT_->Matrix(colblock, rowblock).Add(scon.Matrix(rowblock, colblock), true, 1.0, 0.0);
    }
  }
  sconT_->Complete();


  scon.Complete();

  mat.Assign(0, 3, LINALG::View, scon.Matrix(0, 0));

  (*scgitransform_)(
      scon.Matrix(1, 0), 1. / scale, ADAPTER::CouplingMasterConverter(coupsf), mat.Matrix(1, 3));

  mat.Assign(3, 0, LINALG::View, sconT_->Matrix(0, 0));

  (*csigtransform_)(*coupsf.MasterDofMap(), sconT_->Matrix(0, 1).ColMap(), sconT_->Matrix(0, 1),
      1. / timescale, ADAPTER::CouplingMasterConverter(coupsf), mat.Matrix(3, 1), true);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();

  // Finally, take care of Dirichlet boundary conditions
  mat.ApplyDirichlet(*(dbcmaps_->CondMap()), true);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::InitialGuess");

  Teuchos::RCP<Epetra_Vector> ConstraintInitialGuess =
      Teuchos::rcp(new Epetra_Vector(*(conman_->GetConstraintMap()), true));

  SetupVector(*ig, StructureField()->InitialGuess(), FluidField()->InitialGuess(),
      AleField()->InitialGuess(), ConstraintInitialGuess, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::SetupVector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av, Teuchos::RCP<const Epetra_Vector> cv, double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<Epetra_Vector> sov = StructureField()->Interface()->ExtractOtherVector(sv);
  Teuchos::RCP<Epetra_Vector> aov = AleField()->Interface()->ExtractOtherVector(av);

  if (fluidscale != 0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> scv = StructureField()->Interface()->ExtractFSICondVector(sv);
    Teuchos::RCP<Epetra_Vector> modfv =
        FluidField()->Interface()->InsertFSICondVector(StructToFluid(scv));
    modfv->Update(1.0, *fv, 1. / fluidscale);

    Extractor().InsertVector(*modfv, 1, f);
  }
  else
  {
    Extractor().InsertVector(*fv, 1, f);
  }

  Extractor().InsertVector(*sov, 0, f);
  Extractor().InsertVector(*aov, 2, f);

  Epetra_Vector modcv = *cv;
  modcv.Scale(-1.0);
  Extractor().InsertVector(modcv, 3, f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::ConstrMonolithicStructureSplit::ExtractFieldVectors");

  fx = Extractor().ExtractVector(x, 1);

  // process structure unknowns

  Teuchos::RCP<Epetra_Vector> fcx = FluidField()->Interface()->ExtractFSICondVector(fx);
  FluidField()->VelocityToDisplacement(fcx);
  Teuchos::RCP<const Epetra_Vector> sox = Extractor().ExtractVector(x, 0);
  Teuchos::RCP<Epetra_Vector> scx = FluidToStruct(fcx);

  Teuchos::RCP<Epetra_Vector> s = StructureField()->Interface()->InsertOtherVector(sox);
  StructureField()->Interface()->InsertFSICondVector(scx, s);
  sx = s;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x, 2);
  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);

  Teuchos::RCP<Epetra_Vector> a = AleField()->Interface()->InsertOtherVector(aox);
  AleField()->Interface()->InsertFSICondVector(acx, a);

  ax = a;
}
