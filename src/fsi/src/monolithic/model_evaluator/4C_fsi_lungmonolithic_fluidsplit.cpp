/*----------------------------------------------------------------------*/
/*! \file
\brief Volume-coupled FSI (fluid-split)


\level 3
*/
/*----------------------------------------------------------------------*/
#include "4C_fsi_lungmonolithic_fluidsplit.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_lung.hpp"
#include "4C_adapter_str_lung.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_lung_overlapprec.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_structure_aux.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

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
  intersectionmaps.push_back(fluid_field()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(fluid_field()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      CORE::LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

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

  fggtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowColTransform);
  fgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  fg_gtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  figtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  f_ggtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);

  fmiitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fm_gitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowColTransform);
  fmigtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fm_ggtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fmggtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowColTransform);
  fmi_gtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fm_g_gtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fmg_gtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowColTransform);

  addfm_g_gtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  addfm_ggtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);

  fcgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);

  aigtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  ai_gtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);

  cai_gtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::SetupSystem()
{
  GeneralSetup();

  // create combined map
  create_combined_dof_row_map();

  fluid_field()->use_block_matrix(true);

  // build ale system matrix in splitted system
  ale_field()->create_system_matrix(ale_field()->Interface());


  //-----------------------------------------------------------------------------
  // create block system matrix
  //-----------------------------------------------------------------------------
  create_system_matrix(false);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::create_combined_dof_row_map()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(structure_field()->dof_row_map());
  vecSpaces.push_back(fluid_field()->dof_row_map());
  vecSpaces.push_back(ale_field()->Interface()->OtherMap());
  vecSpaces.push_back(ConstrMap_);

  if (vecSpaces[0]->NumGlobalElements() == 0)
    FOUR_C_THROW("No inner structural equations. Splitting not possible. Panic.");

  set_dof_row_maps(vecSpaces);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::setup_rhs_residual(Epetra_Vector& f)
{
  Teuchos::RCP<Epetra_Vector> structureRHS =
      Teuchos::rcp(new Epetra_Vector(*structure_field()->discretization()->dof_row_map()));
  structureRHS->Update(1.0, *structure_field()->RHS(), 1.0, *AddStructRHS_, 0.0);

  Teuchos::RCP<Epetra_Vector> fluidRHS =
      Teuchos::rcp(new Epetra_Vector(*fluid_field()->discretization()->dof_row_map()));
  fluidRHS->Update(1.0, *fluid_field()->RHS(), 1.0, *AddFluidRHS_, 0.0);

  const double scale = fluid_field()->residual_scaling();

  setup_vector(f, structureRHS, fluidRHS, ale_field()->RHS(), ConstrRHS_, scale);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::setup_rhs_lambda(Epetra_Vector& f)
{
  // ToDo: We still need to implement this.

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::setup_rhs_firstiter(Epetra_Vector& f)
{
  Teuchos::RCP<Epetra_Vector> structureRHS =
      Teuchos::rcp(new Epetra_Vector(*structure_field()->discretization()->dof_row_map()));
  structureRHS->Update(1.0, *structure_field()->RHS(), 1.0, *AddStructRHS_, 0.0);

  const double scale = fluid_field()->residual_scaling();

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> blockf = fluid_field()->BlockSystemMatrix();

  CORE::LINALG::SparseMatrix& fig = blockf->Matrix(0, 1);
  CORE::LINALG::SparseMatrix& fgg = blockf->Matrix(1, 1);
  CORE::LINALG::SparseMatrix& fGg = blockf->Matrix(3, 1);

  Teuchos::RCP<Epetra_Vector> fveln = fluid_field()->extract_interface_veln();
  double timescale = fluid_field()->TimeScaling();

  Teuchos::RCP<ADAPTER::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidLung>(fluid_field());
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
  extractor().AddVector(*rhs, 1, f);

  rhs = Teuchos::rcp(new Epetra_Vector(fgg.RowMap()));
  fgg.Apply(*fveln, *rhs);
  rhs->Scale(scale * timescale * Dt());
  rhs = fluid_to_struct(rhs);
  rhs = structure_field()->Interface()->InsertFSICondVector(rhs);
  extractor().AddVector(*rhs, 0, f);

  //--------------------------------------------------------------------------------
  // constraint fluid
  //--------------------------------------------------------------------------------
  // split in two blocks according to inner and fsi structure dofs
  Teuchos::RCP<Epetra_Map> emptymap =
      Teuchos::rcp(new Epetra_Map(-1, 0, nullptr, 0, structure_field()->discretization()->Comm()));
  CORE::LINALG::MapExtractor extractor_temp;
  extractor_temp.Setup(*ConstrMap_, emptymap, ConstrMap_);

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> constrfluidblocks =
      ConstrFluidMatrix_->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *fluidfield->FSIInterface(), extractor_temp);
  constrfluidblocks->Complete();

  CORE::LINALG::SparseMatrix& cfig = constrfluidblocks->Matrix(0, 1);
  rhs = Teuchos::rcp(new Epetra_Vector(cfig.RowMap()));
  cfig.Apply(*fveln, *rhs);
  rhs->Scale(timescale * Dt());
  extractor().AddVector(*rhs, 3, f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::setup_system_matrix(CORE::LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::setup_system_matrix");

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W

  const CORE::ADAPTER::Coupling& coupsf = structure_fluid_coupling();
  const CORE::ADAPTER::Coupling& coupsa = structure_ale_coupling();
  const CORE::ADAPTER::Coupling& coupfa = fluid_ale_coupling();

  /*----------------------------------------------------------------------*/

  double scale = fluid_field()->residual_scaling();
  double timescale = fluid_field()->TimeScaling();

  /*----------------------------------------------------------------------*/
  // structure part

  CORE::LINALG::SparseMatrix s = *structure_field()->SystemMatrix();
  s.UnComplete();
  s.Add(AddStructConstrMatrix_->Matrix(0, 0), false, 1.0, 1.0);

  mat.Assign(0, 0, CORE::LINALG::View, s);

  /*----------------------------------------------------------------------*/
  // structure constraint part

  CORE::LINALG::SparseMatrix scon = AddStructConstrMatrix_->Matrix(0, 1);
  mat.Assign(0, 3, CORE::LINALG::View, scon);

  /*----------------------------------------------------------------------*/
  // fluid part

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> blockf = fluid_field()->BlockSystemMatrix();

  CORE::LINALG::SparseMatrix& fii = blockf->Matrix(0, 0);
  CORE::LINALG::SparseMatrix& fig = blockf->Matrix(0, 1);
  CORE::LINALG::SparseMatrix& fiG = blockf->Matrix(0, 3);
  CORE::LINALG::SparseMatrix& fgi = blockf->Matrix(1, 0);
  CORE::LINALG::SparseMatrix& fgg = blockf->Matrix(1, 1);
  CORE::LINALG::SparseMatrix& fgG = blockf->Matrix(1, 3);
  CORE::LINALG::SparseMatrix& fGi = blockf->Matrix(3, 0);
  CORE::LINALG::SparseMatrix& fGg = blockf->Matrix(3, 1);
  CORE::LINALG::SparseMatrix& fGG = blockf->Matrix(3, 3);

  // mat.Matrix(0,0).UnComplete();
  (*fggtransform_)(fgg, scale * timescale, CORE::ADAPTER::CouplingSlaveConverter(coupsf),
      CORE::ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(0, 0), true, true);
  (*fgitransform_)(fgi, scale, CORE::ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(0, 1));
  (*fg_gtransform_)(
      fgG, scale, CORE::ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(0, 1), true);

  (*figtransform_)(blockf->FullRowMap(), blockf->FullColMap(), fig, timescale,
      CORE::ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(1, 0));
  (*f_ggtransform_)(blockf->FullRowMap(), blockf->FullColMap(), fGg, timescale,
      CORE::ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(1, 0), true, true);

  mat.Matrix(1, 1).Add(fii, false, 1.0, 0.0);
  mat.Matrix(1, 1).Add(fiG, false, 1.0, 1.0);
  mat.Matrix(1, 1).Add(fGi, false, 1.0, 1.0);
  mat.Matrix(1, 1).Add(fGG, false, 1.0, 1.0);

  Teuchos::RCP<CORE::LINALG::SparseMatrix> eye =
      CORE::LINALG::Eye(*fluid_field()->Interface()->FSICondMap());
  mat.Matrix(1, 1).Add(*eye, false, 1.0, 1.0);

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> mmm = fluid_field()->ShapeDerivatives();

  if (mmm != Teuchos::null)
  {
    CORE::LINALG::SparseMatrix& fmii = mmm->Matrix(0, 0);
    CORE::LINALG::SparseMatrix& fmig = mmm->Matrix(0, 1);
    CORE::LINALG::SparseMatrix& fmiG = mmm->Matrix(0, 3);
    CORE::LINALG::SparseMatrix& fmgi = mmm->Matrix(1, 0);
    CORE::LINALG::SparseMatrix& fmgg = mmm->Matrix(1, 1);
    CORE::LINALG::SparseMatrix& fmgG = mmm->Matrix(1, 3);
    CORE::LINALG::SparseMatrix& fmGi = mmm->Matrix(3, 0);
    CORE::LINALG::SparseMatrix& fmGg = mmm->Matrix(3, 1);
    CORE::LINALG::SparseMatrix& fmGG = mmm->Matrix(3, 3);

    // We cannot copy the pressure value. It is not used anyway. So no exact
    // match here.
    (*fmiitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmii, 1.,
        CORE::ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, false);
    (*fm_gitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmGi, 1.,
        CORE::ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, true);
    (*fmgitransform_)(fmgi, scale, CORE::ADAPTER::CouplingSlaveConverter(coupsf),
        CORE::ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(0, 2), false, false);

    (*fmigtransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmig, 1.,
        CORE::ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(1, 0), true, true);
    (*fm_ggtransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmGg, 1.,
        CORE::ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(1, 0), true, true);

    (*fmggtransform_)(fmgg, scale, CORE::ADAPTER::CouplingSlaveConverter(coupsf),
        CORE::ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(0, 0), true, true);

    (*fmi_gtransform_)(*coupfsout_->MasterDofMap(), fmiG.ColMap(), fmiG, 1.,
        CORE::ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(1, 0), true, true);
    (*fm_g_gtransform_)(*coupfsout_->MasterDofMap(), fmGG.ColMap(), fmGG, 1.,
        CORE::ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(1, 0), true, true);
    (*fmg_gtransform_)(fmgG, scale, CORE::ADAPTER::CouplingSlaveConverter(coupsf),
        CORE::ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(0, 0), true, true);

    CORE::LINALG::SparseMatrix& addfmGg = AddFluidShapeDerivMatrix_->Matrix(3, 1);
    CORE::LINALG::SparseMatrix& addfmGG = AddFluidShapeDerivMatrix_->Matrix(3, 3);

    (*addfm_g_gtransform_)(AddFluidShapeDerivMatrix_->FullRowMap(),
        AddFluidShapeDerivMatrix_->FullColMap(), addfmGG, 1.,
        CORE::ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(1, 0), true, true);

    (*addfm_ggtransform_)(AddFluidShapeDerivMatrix_->FullRowMap(),
        AddFluidShapeDerivMatrix_->FullColMap(), addfmGg, 1.,
        CORE::ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(1, 0), true, true);
  }

  /*----------------------------------------------------------------------*/
  // fluid constraint part

  // split in two blocks according to inner and fsi dofs

  Teuchos::RCP<ADAPTER::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidLung>(fluid_field());

  Teuchos::RCP<Epetra_Map> emptymap =
      Teuchos::rcp(new Epetra_Map(-1, 0, nullptr, 0, fluid_field()->discretization()->Comm()));
  CORE::LINALG::MapExtractor extractor;
  extractor.Setup(*ConstrMap_, emptymap, ConstrMap_);

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> fluidconstrblocks =
      FluidConstrMatrix_->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
          extractor, *fluidfield->FSIInterface());

  fluidconstrblocks->Complete();

  CORE::LINALG::SparseMatrix& fcii = fluidconstrblocks->Matrix(0, 0);
  CORE::LINALG::SparseMatrix& fcgi = fluidconstrblocks->Matrix(1, 0);

  // fcii cannot be simply assigned here (in case of fsi amg which is default)
  // due to non-matching maps
  mat.Matrix(1, 3).Zero();
  mat.Matrix(1, 3).Add(fcii, false, 1.0, 0.0);

  // add interface part to structure block

  mat.Matrix(0, 3).UnComplete();
  (*fcgitransform_)(
      fcgi, scale, CORE::ADAPTER::CouplingSlaveConverter(coupsf), mat.Matrix(0, 3), true);

  /*----------------------------------------------------------------------*/
  // ale part

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> a = ale_field()->BlockSystemMatrix();

  if (a == Teuchos::null) FOUR_C_THROW("expect ale block matrix");

  a->Complete();

  CORE::LINALG::SparseMatrix& aii = a->Matrix(0, 0);
  CORE::LINALG::SparseMatrix& aig = a->Matrix(0, 1);
  CORE::LINALG::SparseMatrix& aiG = a->Matrix(0, 3);

  mat.Assign(2, 2, CORE::LINALG::View, aii);

  (*aigtransform_)(a->FullRowMap(), a->FullColMap(), aig, 1.,
      CORE::ADAPTER::CouplingSlaveConverter(coupsa), mat.Matrix(2, 0));

  (*ai_gtransform_)(a->FullRowMap(), a->FullColMap(), aiG, 1.,
      CORE::ADAPTER::CouplingSlaveConverter(*coupsaout_), mat.Matrix(2, 0), true, true);

  /*----------------------------------------------------------------------*/
  // constraint part -> structure

  mat.Assign(3, 0, CORE::LINALG::View, AddStructConstrMatrix_->Matrix(1, 0));

  /*----------------------------------------------------------------------*/
  // constraint part -> fluid

  // split in two blocks according to inner and fsi structure dofs

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> constrfluidblocks =
      ConstrFluidMatrix_->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *fluidfield->FSIInterface(), extractor);
  constrfluidblocks->Complete();

  CORE::LINALG::SparseMatrix& cfii = constrfluidblocks->Matrix(0, 0);

  // cfii cannot be simply assigned here (in case of fsi amg which is default)
  // due to non-matching maps
  mat.Matrix(3, 1).Zero();
  mat.Matrix(3, 1).Add(cfii, false, 1.0, 0.0);

  /*----------------------------------------------------------------------*/
  // constraint part -> "ale"

  CORE::LINALG::SparseMatrix& caiG = ConstrAleMatrix_->Matrix(0, 3);
  (*cai_gtransform_)(*coupfsout_->MasterDofMap(), caiG.ColMap(), caiG, 1.0,
      CORE::ADAPTER::CouplingMasterConverter(*coupfsout_), mat.Matrix(3, 0), true, true);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::initial_guess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::initial_guess");

  Teuchos::RCP<Epetra_Vector> ConstraintInitialGuess =
      Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));

  setup_vector(*ig, structure_field()->initial_guess(), fluid_field()->initial_guess(),
      ale_field()->initial_guess(), ConstraintInitialGuess, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::setup_vector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av, Teuchos::RCP<const Epetra_Vector> cv, double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<ADAPTER::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidLung>(fluid_field());
  Teuchos::RCP<Epetra_Vector> fov = fluidfield->FSIInterface()->ExtractOtherVector(fv);
  fov = fluidfield->FSIInterface()->InsertOtherVector(fov);
  Teuchos::RCP<Epetra_Vector> aov = ale_field()->Interface()->ExtractOtherVector(av);

  if (fluidscale != 0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> fcv = fluid_field()->Interface()->ExtractFSICondVector(fv);
    Teuchos::RCP<Epetra_Vector> modsv =
        structure_field()->Interface()->InsertFSICondVector(fluid_to_struct(fcv));
    modsv->Update(1.0, *sv, fluidscale);

    extractor().InsertVector(*modsv, 0, f);
  }
  else
  {
    extractor().InsertVector(*sv, 0, f);
  }

  extractor().InsertVector(*fov, 1, f);
  extractor().InsertVector(*aov, 2, f);
  extractor().InsertVector(*cv, 3, f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::LungMonolithicFluidSplit::extract_field_vectors");

  Teuchos::RCP<ADAPTER::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidLung>(fluid_field());

  // We have overlap at the interface. Thus we need the interface part of the
  // structure vector and append it to the fluid and ale vector. (With the
  // right translation.)

  sx = extractor().ExtractVector(x, 0);
  Teuchos::RCP<const Epetra_Vector> scx = structure_field()->Interface()->ExtractFSICondVector(sx);

  // process fluid unknowns

  Teuchos::RCP<const Epetra_Vector> fox = extractor().ExtractVector(x, 1);
  fox = fluidfield->FSIInterface()->ExtractOtherVector(fox);
  Teuchos::RCP<Epetra_Vector> fcx = struct_to_fluid(scx);

  fluid_field()->displacement_to_velocity(fcx);

  Teuchos::RCP<Epetra_Vector> f = fluidfield->FSIInterface()->InsertOtherVector(fox);
  fluid_field()->Interface()->InsertFSICondVector(fcx, f);
  fx = f;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = extractor().ExtractVector(x, 2);
  Teuchos::RCP<Epetra_Vector> acx = struct_to_ale(scx);
  Teuchos::RCP<Epetra_Vector> a = ale_field()->Interface()->InsertOtherVector(aox);
  ale_field()->Interface()->InsertVector(acx, 1, a);

  Teuchos::RCP<Epetra_Vector> scox = structure_field()->Interface()->ExtractLungASICondVector(sx);
  Teuchos::RCP<Epetra_Vector> acox = struct_to_ale_outflow(scox);
  ale_field()->Interface()->InsertVector(acox, 3, a);

  ax = a;
}

FOUR_C_NAMESPACE_CLOSE
