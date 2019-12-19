/*----------------------------------------------------------------------*/
/*! \file
\brief Volume-coupled FSI (structure-split)

\maintainer Matthias Mayr

\level 3
*/
/*----------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "fsi_lungmonolithic_structuresplit.H"
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
FSI::LungMonolithicStructureSplit::LungMonolithicStructureSplit(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : LungMonolithic(comm, timeparams)
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

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::SetupSystem()
{
  GeneralSetup();
  // create combined map
  CreateCombinedDofRowMap();

  FluidField()->UseBlockMatrix(false);

  // build ale system matrix in splitted system
  AleField()->CreateSystemMatrix(AleField()->Interface());

  //-----------------------------------------------------------------------------
  // create block system matrix
  //-----------------------------------------------------------------------------
  CreateSystemMatrix(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::CreateCombinedDofRowMap()
{
  const Teuchos::RCP<ADAPTER::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());

  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(structfield->FSIInterface()->OtherMap());
  vecSpaces.push_back(FluidField()->DofRowMap());
  // remaining (not coupled) dofs of ale field
  vecSpaces.push_back(AleField()->Interface()->Map(0));
  // additional volume constraints
  vecSpaces.push_back(ConstrMap_);

  if (vecSpaces[0]->NumGlobalElements() == 0)
    dserror("No inner structural equations. Splitting not possible. Panic.");

  SetDofRowMaps(vecSpaces);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::SetupRHSResidual(Epetra_Vector& f)
{
  Teuchos::RCP<Epetra_Vector> structureRHS =
      Teuchos::rcp(new Epetra_Vector(*StructureField()->Discretization()->DofRowMap()));
  structureRHS->Update(1.0, *StructureField()->RHS(), 1.0, *AddStructRHS_, 0.0);

  const double fluidscale = FluidField()->ResidualScaling();

  Teuchos::RCP<Epetra_Vector> fluidRHS =
      Teuchos::rcp(new Epetra_Vector(*FluidField()->Discretization()->DofRowMap()));
  fluidRHS->Update(1.0, *FluidField()->RHS(), 1.0, *AddFluidRHS_, 0.0);

  SetupVector(f, structureRHS, fluidRHS, AleField()->RHS(), ConstrRHS_, fluidscale);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::SetupRHSLambda(Epetra_Vector& f)
{
  // ToDo: We still need to implement this.

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::SetupRHSFirstiter(Epetra_Vector& f)
{
  const Teuchos::RCP<ADAPTER::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());
  const Teuchos::RCP<ADAPTER::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidLung>(FluidField());

  const double fluidscale = FluidField()->ResidualScaling();

  // additional rhs term for ALE equations
  // -dt Aig u(n)
  //
  //    1/dt Delta d(n+1) = theta Delta u(n+1) + u(n)
  //
  // And we are concerned with the u(n) part here.

  //--------------------------------------------------------------------------------
  // ale
  //--------------------------------------------------------------------------------
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField()->BlockSystemMatrix();
  if (a == Teuchos::null) dserror("expect ale block matrix");

  LINALG::SparseMatrix& aig = a->Matrix(0, 1);

  Teuchos::RCP<Epetra_Vector> fveln = FluidField()->ExtractInterfaceVeln();
  Teuchos::RCP<Epetra_Vector> sveln = FluidToStruct(fveln);
  Teuchos::RCP<Epetra_Vector> aveln = StructToAle(sveln);
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap()));
  aig.Apply(*aveln, *rhs);

  rhs->Scale(-1. * Dt());

  Extractor().AddVector(*rhs, 2, f);

  //--------------------------------------------------------------------------------
  // structure
  //--------------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> veln = StructureField()->Interface()->InsertFSICondVector(sveln);
  rhs = Teuchos::rcp(new Epetra_Vector(veln->Map()));

  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField()->SystemMatrix();
  s->Apply(*veln, *rhs);

  Teuchos::RCP<Epetra_Vector> addrhs = Teuchos::rcp(new Epetra_Vector(veln->Map()));
  AddStructConstrMatrix_->Matrix(0, 0).Apply(*veln, *addrhs);

  rhs->Update(1.0, *addrhs, 1.0);
  rhs->Scale(-1. * Dt());

  veln = structfield->FSIInterface()->ExtractOtherVector(rhs);
  Extractor().AddVector(*veln, 0, f);

  veln = StructureField()->Interface()->ExtractFSICondVector(rhs);
  veln = FluidField()->Interface()->InsertFSICondVector(StructToFluid(veln));

  veln->Scale(1. / fluidscale);
  Extractor().AddVector(*veln, 1, f);

  //--------------------------------------------------------------------------------
  // constraint structure
  //--------------------------------------------------------------------------------
  // split in two blocks according to inner and fsi structure dofs
  Teuchos::RCP<Epetra_Map> emptymap =
      Teuchos::rcp(new Epetra_Map(-1, 0, NULL, 0, StructureField()->Discretization()->Comm()));
  LINALG::MapExtractor extractor;
  extractor.Setup(*ConstrMap_, emptymap, ConstrMap_);

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> constrstructblocks =
      AddStructConstrMatrix_->Matrix(1, 0).Split<LINALG::DefaultBlockMatrixStrategy>(
          *structfield->FSIInterface(), extractor);
  constrstructblocks->Complete();

  LINALG::SparseMatrix& csig = constrstructblocks->Matrix(0, 1);

  rhs = Teuchos::rcp(new Epetra_Vector(csig.RowMap()));
  csig.Apply(*sveln, *rhs);
  rhs->Scale(-1. * Dt());
  Extractor().AddVector(*rhs, 3, f);

  //--------------------------------------------------------------------------------
  // constraint ale
  //--------------------------------------------------------------------------------
  LINALG::SparseMatrix& caig = ConstrAleMatrix_->Matrix(0, 1);
  caig.Apply(*fveln, *rhs);
  rhs->Scale(-1. * Dt());
  Extractor().AddVector(*rhs, 3, f);

  //--------------------------------------------------------------------------------
  // shape derivatives
  //--------------------------------------------------------------------------------
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = fluidfield->ShapeDerivatives();

  if (mmm != Teuchos::null)
  {
    LINALG::SparseMatrix& fmig = mmm->Matrix(0, 1);
    LINALG::SparseMatrix& fmgg = mmm->Matrix(1, 1);
    LINALG::SparseMatrix& fmGg = mmm->Matrix(3, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap()));
    fmig.Apply(*fveln, *rhs);
    veln = FluidField()->Interface()->InsertVector(rhs, 0);

    rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap()));
    fmgg.Apply(*fveln, *rhs);
    FluidField()->Interface()->InsertVector(rhs, 1, veln);

    rhs = Teuchos::rcp(new Epetra_Vector(fmGg.RowMap()));
    fmGg.Apply(*fveln, *rhs);
    FluidField()->Interface()->InsertVector(rhs, 3, veln);

    veln->Scale(-1. * Dt());
    Extractor().AddVector(*veln, 1, f);

    LINALG::SparseMatrix& afmgg = AddFluidShapeDerivMatrix_->Matrix(1, 1);
    LINALG::SparseMatrix& afmGg = AddFluidShapeDerivMatrix_->Matrix(3, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(afmgg.RowMap()));
    afmgg.Apply(*fveln, *rhs);
    veln = FluidField()->Interface()->InsertVector(rhs, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(afmGg.RowMap()));
    afmGg.Apply(*fveln, *rhs);
    FluidField()->Interface()->InsertVector(rhs, 3, veln);

    veln->Scale(-1. * Dt());
    Extractor().AddVector(*veln, 1, f);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::SetupSystemMatrix");

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();

  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField()->SystemMatrix();

  double scale = FluidField()->ResidualScaling();
  double timescale = FluidField()->TimeScaling();

  /*----------------------------------------------------------------------*/
  // fluid part
  mat.Assign(1, 1, LINALG::View, *f);

  // fluid linearization with respect to mesh motion block
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField()->ShapeDerivatives();
  const ADAPTER::Coupling& coupfa = FluidAleCoupling();

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

    LINALG::SparseMatrix& addfmGg = AddFluidShapeDerivMatrix_->Matrix(3, 1);
    LINALG::SparseMatrix& addfmGG = AddFluidShapeDerivMatrix_->Matrix(3, 3);

    mat.Matrix(1, 1).Add(fmig, false, 1. / timescale, 1.0);
    mat.Matrix(1, 1).Add(fmgg, false, 1. / timescale, 1.0);
    mat.Matrix(1, 1).Add(fmGg, false, 1. / timescale, 1.0);

    fmiitransform_(mmm->FullRowMap(), mmm->FullColMap(), fmii, 1.,
        ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, false);
    fmgitransform_(mmm->FullRowMap(), mmm->FullColMap(), fmgi, 1.,
        ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, true);
    fmGitransform_(mmm->FullRowMap(), mmm->FullColMap(), fmGi, 1.,
        ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, true);

    fmgGtransform_(mmm->FullRowMap(), mmm->FullColMap(), fmgG, 1.,
        ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(1, 0), true, false);
    fmiGtransform_(mmm->FullRowMap(), mmm->FullColMap(), fmiG, 1.,
        ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(1, 0), true, true);
    fmGGtransform_(mmm->FullRowMap(), mmm->FullColMap(), fmGG, 1.,
        ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(1, 0), true, true);

    mat.Matrix(1, 1).Add(addfmGg, false, 1. / timescale, 1.0);

    addfmGGtransform_(AddFluidShapeDerivMatrix_->FullRowMap(),
        AddFluidShapeDerivMatrix_->FullColMap(), addfmGG, 1.,
        ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(1, 0), true, true);
  }

  /*----------------------------------------------------------------------*/
  // structure and additional structure part

  const Teuchos::RCP<ADAPTER::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());

  LINALG::SparseMatrix& s = *StructureField()->SystemMatrix();
  s.Add(AddStructConstrMatrix_->Matrix(0, 0), false, 1.0, 1.0);

  siitransform_(s, *structfield->FSIInterface()->Map(0), *structfield->FSIInterface()->Map(0), 1.0,
      NULL, NULL, mat.Matrix(0, 0), true, false);

  ADAPTER::CouplingMasterConverter converter(coupsf);
  sigtransform_(s, *structfield->FSIInterface()->Map(0), *structfield->FSIInterface()->Map(1),
      1. / timescale, NULL, &converter, mat.Matrix(0, 1), true, false);

  sgitransform_(s, *structfield->FSIInterface()->Map(1), *structfield->FSIInterface()->Map(0),
      1. / scale, &converter, NULL, mat.Matrix(1, 0), true, true);

  sggtransform_(s, *structfield->FSIInterface()->Map(1), *structfield->FSIInterface()->Map(1),
      1. / (scale * timescale), &converter, &converter, mat.Matrix(1, 1), true, true);

  /*----------------------------------------------------------------------*/
  // fluid constraint part
  // add into new matrix instead of view -> avoids matrix UnComplete() operations in the
  // second and subsequent time steps
  mat.Matrix(1, 3).Add(*FluidConstrMatrix_, false, 1., 0.);

  /*----------------------------------------------------------------------*/
  // structure constraint part

  // split in two blocks according to inner and fsi structure dofs
  LINALG::MapExtractor extractor;

  sciitransform_(AddStructConstrMatrix_->Matrix(0, 1), *structfield->FSIInterface()->Map(0),
      *ConstrMap_, 1., NULL, NULL, mat.Matrix(0, 3), true, false);

  // add interface part to fluid block
  scgitransform_(AddStructConstrMatrix_->Matrix(0, 1), *structfield->FSIInterface()->Map(1),
      *ConstrMap_, 1. / scale, &converter, NULL, mat.Matrix(1, 3), true, true);

  /*----------------------------------------------------------------------*/
  // ale part
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField()->BlockSystemMatrix();

  if (a == Teuchos::null) dserror("expect ale block matrix");

  a->Complete();

  LINALG::SparseMatrix& aii = a->Matrix(0, 0);
  LINALG::SparseMatrix& aig = a->Matrix(0, 1);
  LINALG::SparseMatrix& aiG = a->Matrix(0, 3);

  aiGtransform_(a->FullRowMap(), a->FullColMap(), aiG, 1.,
      ADAPTER::CouplingSlaveConverter(*coupsaout_), mat.Matrix(2, 0));
  aigtransform_(a->FullRowMap(), a->FullColMap(), aig, 1. / timescale,
      ADAPTER::CouplingSlaveConverter(*icoupfa_), mat.Matrix(2, 1));
  mat.Assign(2, 2, LINALG::View, aii);

  /*----------------------------------------------------------------------*/
  // constraint part -> fluid
  // add into new matrix instead of view -> avoids matrix UnComplete() operations in the
  // second and subsequent time steps

  mat.Matrix(3, 1).Add(*ConstrFluidMatrix_, false, 1., 0.);

  /*----------------------------------------------------------------------*/
  // constraint part -> structure
  // split in two blocks according to inner and fsi structure dofs

  csiitransform_(AddStructConstrMatrix_->Matrix(1, 0), *ConstrMap_,
      *structfield->FSIInterface()->Map(0), 1., NULL, NULL, mat.Matrix(3, 0), true, false);
  csigtransform_(AddStructConstrMatrix_->Matrix(1, 0), *ConstrMap_,
      *structfield->FSIInterface()->Map(1), 1. / timescale, NULL, &converter, mat.Matrix(3, 1),
      true, true);

  /*----------------------------------------------------------------------*/
  // constraint part -> ale

  LINALG::SparseMatrix& caiG = ConstrAleMatrix_->Matrix(0, 3);
  caiGtransform_(*coupfsout_->MasterDofMap(), caiG.ColMap(), caiG, 1.0,
      ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(3, 0), true, true);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::InitialGuess");

  Teuchos::RCP<Epetra_Vector> ConstraintInitialGuess =
      Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));

  SetupVector(*ig, StructureField()->InitialGuess(), FluidField()->InitialGuess(),
      AleField()->InitialGuess(), ConstraintInitialGuess, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::SetupVector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av, Teuchos::RCP<const Epetra_Vector> cv, double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  const Teuchos::RCP<ADAPTER::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());

  Teuchos::RCP<Epetra_Vector> sov = structfield->FSIInterface()->ExtractOtherVector(sv);

  Teuchos::RCP<Epetra_Vector> aov = AleField()->Interface()->ExtractVector(av, 0);

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
  Extractor().InsertVector(*cv, 3, f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  const Teuchos::RCP<ADAPTER::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());

  fx = Extractor().ExtractVector(x, 1);

  // process structure unknowns

  Teuchos::RCP<Epetra_Vector> fcx = FluidField()->Interface()->ExtractFSICondVector(fx);
  FluidField()->VelocityToDisplacement(fcx);
  Teuchos::RCP<const Epetra_Vector> sox = Extractor().ExtractVector(x, 0);
  Teuchos::RCP<Epetra_Vector> scx = FluidToStruct(fcx);

  Teuchos::RCP<Epetra_Vector> s = structfield->FSIInterface()->InsertOtherVector(sox);
  StructureField()->Interface()->InsertFSICondVector(scx, s);
  sx = s;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x, 2);
  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);
  Teuchos::RCP<Epetra_Vector> a = AleField()->Interface()->InsertVector(aox, 0);
  AleField()->Interface()->InsertVector(acx, 1, a);

  Teuchos::RCP<Epetra_Vector> scox = StructureField()->Interface()->ExtractLungASICondVector(sx);
  Teuchos::RCP<Epetra_Vector> acox = StructToAleOutflow(scox);
  AleField()->Interface()->InsertVector(acox, 3, a);

  ax = a;
}
