/*----------------------------------------------------------------------*/
/*!
\file fsi_lungmonolithic_fluidsplit.cpp
\brief Volume-coupled FSI (fluid-split)

\maintainer Matthias Mayr

\level 3
*/
/*----------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "fsi_lungmonolithic_fluidsplit.H"
#include "fsi_matrixtransform.H"
#include "fsi_lung_overlapprec.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_adapter/ad_str_lung.H"
#include "../drt_adapter/ad_fld_lung.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_io/io_control.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_adapter/ad_ale_fsi.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::LungMonolithicFluidSplit::LungMonolithicFluidSplit(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : LungMonolithic(comm, timeparams)
{
  // ---------------------------------------------------------------------------
  // FSI specific check of Dirichlet boundary conditions
  // ---------------------------------------------------------------------------
  // Create intersection of slave DOFs that hold a Dirichlet boundary condition
  // and are located at the FSI interface
  std::vector<Teuchos::RCP<const Epetra_Map>> intersectionmaps;
  intersectionmaps.push_back(FluidField()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(FluidField()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  // Check whether the intersection is empty
  if (intersectionmap->NumGlobalElements() != 0)
  {
    // It is not allowed, that slave DOFs at the interface hold a Dirichlet
    // boundary condition. Thus --> ToDO: Error message

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
             << "  | This is a fluid split scheme. Hence, master and slave field are chosen as "
                "follows:          |"
             << std::endl
             << "  |     MASTER  = STRUCTURE                                                       "
                "              |"
             << std::endl
             << "  |     SLAVE   = FLUID                                                           "
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

  fggtransform_ = Teuchos::rcp(new UTILS::MatrixRowColTransform);
  fgitransform_ = Teuchos::rcp(new UTILS::MatrixRowTransform);
  fgGtransform_ = Teuchos::rcp(new UTILS::MatrixRowTransform);
  figtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fGgtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  fmiitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmGitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new UTILS::MatrixRowColTransform);
  fmigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmGgtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmggtransform_ = Teuchos::rcp(new UTILS::MatrixRowColTransform);
  fmiGtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmGGtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmgGtransform_ = Teuchos::rcp(new UTILS::MatrixRowColTransform);

  addfmGGtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  addfmGgtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  fcgitransform_ = Teuchos::rcp(new UTILS::MatrixRowTransform);

  aigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  aiGtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  caiGtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::SetupSystem()
{
  GeneralSetup();

  // create combined map
  CreateCombinedDofRowMap();

  FluidField()->UseBlockMatrix(true);

  // build ale system matrix in splitted system
  AleField()->CreateSystemMatrix(AleField()->Interface());


  //-----------------------------------------------------------------------------
  // create block system matrix
  //-----------------------------------------------------------------------------
  CreateSystemMatrix(false);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::CreateCombinedDofRowMap()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(StructureField()->DofRowMap());
  vecSpaces.push_back(FluidField()->DofRowMap());
  vecSpaces.push_back(AleField()->Interface()->OtherMap());
  vecSpaces.push_back(ConstrMap_);

  if (vecSpaces[0]->NumGlobalElements() == 0)
    dserror("No inner structural equations. Splitting not possible. Panic.");

  SetDofRowMaps(vecSpaces);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::SetupRHSResidual(Epetra_Vector& f)
{
  Teuchos::RCP<Epetra_Vector> structureRHS =
      Teuchos::rcp(new Epetra_Vector(*StructureField()->Discretization()->DofRowMap()));
  structureRHS->Update(1.0, *StructureField()->RHS(), 1.0, *AddStructRHS_, 0.0);

  Teuchos::RCP<Epetra_Vector> fluidRHS =
      Teuchos::rcp(new Epetra_Vector(*FluidField()->Discretization()->DofRowMap()));
  fluidRHS->Update(1.0, *FluidField()->RHS(), 1.0, *AddFluidRHS_, 0.0);

  const double scale = FluidField()->ResidualScaling();

  SetupVector(f, structureRHS, fluidRHS, AleField()->RHS(), ConstrRHS_, scale);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::SetupRHSLambda(Epetra_Vector& f)
{
  // ToDo: We still need to implement this.

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::SetupRHSFirstiter(Epetra_Vector& f)
{
  Teuchos::RCP<Epetra_Vector> structureRHS =
      Teuchos::rcp(new Epetra_Vector(*StructureField()->Discretization()->DofRowMap()));
  structureRHS->Update(1.0, *StructureField()->RHS(), 1.0, *AddStructRHS_, 0.0);

  const double scale = FluidField()->ResidualScaling();

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf = FluidField()->BlockSystemMatrix();

  LINALG::SparseMatrix& fig = blockf->Matrix(0, 1);
  LINALG::SparseMatrix& fgg = blockf->Matrix(1, 1);
  LINALG::SparseMatrix& fGg = blockf->Matrix(3, 1);

  Teuchos::RCP<Epetra_Vector> fveln = FluidField()->ExtractInterfaceVeln();
  double timescale = FluidField()->TimeScaling();

  Teuchos::RCP<ADAPTER::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidLung>(FluidField());
  Teuchos::RCP<Epetra_Vector> rhs =
      Teuchos::rcp(new Epetra_Vector(*fluidfield->InnerSplit()->FullMap()));
  Teuchos::RCP<Epetra_Vector> tmp = Teuchos::rcp(new Epetra_Vector(fig.RowMap()));
  fig.Apply(*fveln, *tmp);
  fluidfield->InnerSplit()->InsertOtherVector(tmp, rhs);
  tmp = Teuchos::rcp(new Epetra_Vector(fGg.RowMap()));
  fGg.Apply(*fveln, *tmp);
  fluidfield->InnerSplit()->InsertCondVector(tmp, rhs);
  rhs->Scale(timescale * Dt());
  rhs = fluidfield->FSIInterface()->InsertOtherVector(rhs);
  Extractor().AddVector(*rhs, 1, f);

  rhs = Teuchos::rcp(new Epetra_Vector(fgg.RowMap()));
  fgg.Apply(*fveln, *rhs);
  rhs->Scale(scale * timescale * Dt());
  rhs = FluidToStruct(rhs);
  rhs = StructureField()->Interface()->InsertFSICondVector(rhs);
  Extractor().AddVector(*rhs, 0, f);

  //--------------------------------------------------------------------------------
  // constraint fluid
  //--------------------------------------------------------------------------------
  // split in two blocks according to inner and fsi structure dofs
  Teuchos::RCP<Epetra_Map> emptymap =
      Teuchos::rcp(new Epetra_Map(-1, 0, NULL, 0, StructureField()->Discretization()->Comm()));
  LINALG::MapExtractor extractor;
  extractor.Setup(*ConstrMap_, emptymap, ConstrMap_);

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> constrfluidblocks =
      ConstrFluidMatrix_->Split<LINALG::DefaultBlockMatrixStrategy>(
          *fluidfield->FSIInterface(), extractor);
  constrfluidblocks->Complete();

  LINALG::SparseMatrix& cfig = constrfluidblocks->Matrix(0, 1);
  rhs = Teuchos::rcp(new Epetra_Vector(cfig.RowMap()));
  cfig.Apply(*fveln, *rhs);
  rhs->Scale(timescale * Dt());
  Extractor().AddVector(*rhs, 3, f);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::SetupSystemMatrix");

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  const ADAPTER::Coupling& coupsa = StructureAleCoupling();
  const ADAPTER::Coupling& coupfa = FluidAleCoupling();

  /*----------------------------------------------------------------------*/

  double scale = FluidField()->ResidualScaling();
  double timescale = FluidField()->TimeScaling();

  /*----------------------------------------------------------------------*/
  // structure part

  LINALG::SparseMatrix s = *StructureField()->SystemMatrix();
  s.UnComplete();
  s.Add(AddStructConstrMatrix_->Matrix(0, 0), false, 1.0, 1.0);

  mat.Assign(0, 0, LINALG::View, s);

  /*----------------------------------------------------------------------*/
  // structure constraint part

  LINALG::SparseMatrix scon = AddStructConstrMatrix_->Matrix(0, 1);
  mat.Assign(0, 3, LINALG::View, scon);

  /*----------------------------------------------------------------------*/
  // fluid part

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf = FluidField()->BlockSystemMatrix();

  LINALG::SparseMatrix& fii = blockf->Matrix(0, 0);
  LINALG::SparseMatrix& fig = blockf->Matrix(0, 1);
  LINALG::SparseMatrix& fiG = blockf->Matrix(0, 3);
  LINALG::SparseMatrix& fgi = blockf->Matrix(1, 0);
  LINALG::SparseMatrix& fgg = blockf->Matrix(1, 1);
  LINALG::SparseMatrix& fgG = blockf->Matrix(1, 3);
  LINALG::SparseMatrix& fGi = blockf->Matrix(3, 0);
  LINALG::SparseMatrix& fGg = blockf->Matrix(3, 1);
  LINALG::SparseMatrix& fGG = blockf->Matrix(3, 3);

  // mat.Matrix(0,0).UnComplete();
  (*fggtransform_)(fgg, scale * timescale, ADAPTER::CouplingSlaveConverter(coupsf),
      ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(0, 0), true, true);
  (*fgitransform_)(fgi, scale, ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(0, 1));
  (*fgGtransform_)(fgG, scale, ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(0, 1), true);

  (*figtransform_)(blockf->FullRowMap(), blockf->FullColMap(), fig, timescale,
      ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(1, 0));
  (*fGgtransform_)(blockf->FullRowMap(), blockf->FullColMap(), fGg, timescale,
      ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(1, 0), true, true);

  mat.Matrix(1, 1).Add(fii, false, 1.0, 0.0);
  mat.Matrix(1, 1).Add(fiG, false, 1.0, 1.0);
  mat.Matrix(1, 1).Add(fGi, false, 1.0, 1.0);
  mat.Matrix(1, 1).Add(fGG, false, 1.0, 1.0);

  Teuchos::RCP<LINALG::SparseMatrix> eye = LINALG::Eye(*FluidField()->Interface()->FSICondMap());
  mat.Matrix(1, 1).Add(*eye, false, 1.0, 1.0);

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField()->ShapeDerivatives();

  if (mmm != Teuchos::null)
  {
    LINALG::SparseMatrix& fmii = mmm->Matrix(0, 0);
    LINALG::SparseMatrix& fmig = mmm->Matrix(0, 1);
    LINALG::SparseMatrix& fmiG = mmm->Matrix(0, 3);
    LINALG::SparseMatrix& fmgi = mmm->Matrix(1, 0);
    LINALG::SparseMatrix& fmgg = mmm->Matrix(1, 1);
    LINALG::SparseMatrix& fmgG = mmm->Matrix(1, 3);
    LINALG::SparseMatrix& fmGi = mmm->Matrix(3, 0);
    LINALG::SparseMatrix& fmGg = mmm->Matrix(3, 1);
    LINALG::SparseMatrix& fmGG = mmm->Matrix(3, 3);

    // We cannot copy the pressure value. It is not used anyway. So no exact
    // match here.
    (*fmiitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmii, 1.,
        ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, false);
    (*fmGitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmGi, 1.,
        ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, true);
    (*fmgitransform_)(fmgi, scale, ADAPTER::CouplingSlaveConverter(coupsf),
        ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(0, 2), false, false);

    (*fmigtransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmig, 1.,
        ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(1, 0), true, true);
    (*fmGgtransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmGg, 1.,
        ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(1, 0), true, true);

    (*fmggtransform_)(fmgg, scale, ADAPTER::CouplingSlaveConverter(coupsf),
        ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(0, 0), true, true);

    (*fmiGtransform_)(*coupfsout_->MasterDofMap(), fmiG.ColMap(), fmiG, 1.,
        ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(1, 0), true, true);
    (*fmGGtransform_)(*coupfsout_->MasterDofMap(), fmGG.ColMap(), fmGG, 1.,
        ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(1, 0), true, true);
    (*fmgGtransform_)(fmgG, scale, ADAPTER::CouplingSlaveConverter(coupsf),
        ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(0, 0), true, true);

    LINALG::SparseMatrix& addfmGg = AddFluidShapeDerivMatrix_->Matrix(3, 1);
    LINALG::SparseMatrix& addfmGG = AddFluidShapeDerivMatrix_->Matrix(3, 3);

    (*addfmGGtransform_)(AddFluidShapeDerivMatrix_->FullRowMap(),
        AddFluidShapeDerivMatrix_->FullColMap(), addfmGG, 1.,
        ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(1, 0), true, true);

    (*addfmGgtransform_)(AddFluidShapeDerivMatrix_->FullRowMap(),
        AddFluidShapeDerivMatrix_->FullColMap(), addfmGg, 1.,
        ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(1, 0), true, true);
  }

  /*----------------------------------------------------------------------*/
  // fluid constraint part

  // split in two blocks according to inner and fsi dofs

  Teuchos::RCP<ADAPTER::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidLung>(FluidField());

  Teuchos::RCP<Epetra_Map> emptymap =
      Teuchos::rcp(new Epetra_Map(-1, 0, NULL, 0, FluidField()->Discretization()->Comm()));
  LINALG::MapExtractor extractor;
  extractor.Setup(*ConstrMap_, emptymap, ConstrMap_);

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> fluidconstrblocks =
      FluidConstrMatrix_->Split<LINALG::DefaultBlockMatrixStrategy>(
          extractor, *fluidfield->FSIInterface());

  fluidconstrblocks->Complete();

  LINALG::SparseMatrix& fcii = fluidconstrblocks->Matrix(0, 0);
  LINALG::SparseMatrix& fcgi = fluidconstrblocks->Matrix(1, 0);

  // fcii cannot be simply assigned here (in case of fsi amg which is default)
  // due to non-matching maps
  mat.Matrix(1, 3).Zero();
  mat.Matrix(1, 3).Add(fcii, false, 1.0, 0.0);

  // add interface part to structure block

  mat.Matrix(0, 3).UnComplete();
  (*fcgitransform_)(fcgi, scale, ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(0, 3), true);

  /*----------------------------------------------------------------------*/
  // ale part

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField()->BlockSystemMatrix();

  if (a == Teuchos::null) dserror("expect ale block matrix");

  a->Complete();

  LINALG::SparseMatrix& aii = a->Matrix(0, 0);
  LINALG::SparseMatrix& aig = a->Matrix(0, 1);
  LINALG::SparseMatrix& aiG = a->Matrix(0, 3);

  mat.Assign(2, 2, LINALG::View, aii);

  (*aigtransform_)(a->FullRowMap(), a->FullColMap(), aig, 1.,
      ADAPTER::CouplingSlaveConverter(coupsa), mat.Matrix(2, 0));

  (*aiGtransform_)(a->FullRowMap(), a->FullColMap(), aiG, 1.,
      ADAPTER::CouplingSlaveConverter(*coupsaout_), mat.Matrix(2, 0), true, true);

  /*----------------------------------------------------------------------*/
  // constraint part -> structure

  mat.Assign(3, 0, LINALG::View, AddStructConstrMatrix_->Matrix(1, 0));

  /*----------------------------------------------------------------------*/
  // constraint part -> fluid

  // split in two blocks according to inner and fsi structure dofs

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> constrfluidblocks =
      ConstrFluidMatrix_->Split<LINALG::DefaultBlockMatrixStrategy>(
          *fluidfield->FSIInterface(), extractor);
  constrfluidblocks->Complete();

  LINALG::SparseMatrix& cfii = constrfluidblocks->Matrix(0, 0);

  // cfii cannot be simply assigned here (in case of fsi amg which is default)
  // due to non-matching maps
  mat.Matrix(3, 1).Zero();
  mat.Matrix(3, 1).Add(cfii, false, 1.0, 0.0);

  /*----------------------------------------------------------------------*/
  // constraint part -> "ale"

  LINALG::SparseMatrix& caiG = ConstrAleMatrix_->Matrix(0, 3);
  (*caiGtransform_)(*coupfsout_->MasterDofMap(), caiG.ColMap(), caiG, 1.0,
      ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(3, 0), true, true);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::InitialGuess");

  Teuchos::RCP<Epetra_Vector> ConstraintInitialGuess =
      Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));

  SetupVector(*ig, StructureField()->InitialGuess(), FluidField()->InitialGuess(),
      AleField()->InitialGuess(), ConstraintInitialGuess, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::SetupVector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av, Teuchos::RCP<const Epetra_Vector> cv, double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<ADAPTER::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidLung>(FluidField());
  Teuchos::RCP<Epetra_Vector> fov = fluidfield->FSIInterface()->ExtractOtherVector(fv);
  fov = fluidfield->FSIInterface()->InsertOtherVector(fov);
  Teuchos::RCP<Epetra_Vector> aov = AleField()->Interface()->ExtractOtherVector(av);

  if (fluidscale != 0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> fcv = FluidField()->Interface()->ExtractFSICondVector(fv);
    Teuchos::RCP<Epetra_Vector> modsv =
        StructureField()->Interface()->InsertFSICondVector(FluidToStruct(fcv));
    modsv->Update(1.0, *sv, fluidscale);

    Extractor().InsertVector(*modsv, 0, f);
  }
  else
  {
    Extractor().InsertVector(*sv, 0, f);
  }

  Extractor().InsertVector(*fov, 1, f);
  Extractor().InsertVector(*aov, 2, f);
  Extractor().InsertVector(*cv, 3, f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::LungMonolithicFluidSplit::ExtractFieldVectors");

  Teuchos::RCP<ADAPTER::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidLung>(FluidField());

  // We have overlap at the interface. Thus we need the interface part of the
  // structure vector and append it to the fluid and ale vector. (With the
  // right translation.)

  sx = Extractor().ExtractVector(x, 0);
  Teuchos::RCP<const Epetra_Vector> scx = StructureField()->Interface()->ExtractFSICondVector(sx);

  // process fluid unknowns

  Teuchos::RCP<const Epetra_Vector> fox = Extractor().ExtractVector(x, 1);
  fox = fluidfield->FSIInterface()->ExtractOtherVector(fox);
  Teuchos::RCP<Epetra_Vector> fcx = StructToFluid(scx);

  FluidField()->DisplacementToVelocity(fcx);

  Teuchos::RCP<Epetra_Vector> f = fluidfield->FSIInterface()->InsertOtherVector(fox);
  FluidField()->Interface()->InsertFSICondVector(fcx, f);
  fx = f;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x, 2);
  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);
  Teuchos::RCP<Epetra_Vector> a = AleField()->Interface()->InsertOtherVector(aox);
  AleField()->Interface()->InsertVector(acx, 1, a);

  Teuchos::RCP<Epetra_Vector> scox = StructureField()->Interface()->ExtractLungASICondVector(sx);
  Teuchos::RCP<Epetra_Vector> acox = StructToAleOutflow(scox);
  AleField()->Interface()->InsertVector(acox, 3, a);

  ax = a;
}
