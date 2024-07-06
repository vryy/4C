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
#include "4C_fem_discretization.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_lung_overlapprec.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
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
  intersectionmaps.push_back(fluid_field()->get_dbc_map_extractor()->cond_map());
  intersectionmaps.push_back(fluid_field()->interface()->fsi_cond_map());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      Core::LinAlg::MultiMapExtractor::intersect_maps(intersectionmaps);

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

  fggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
  fgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  fg_gtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  figtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  f_ggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);

  fmiitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fm_gitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
  fmigtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fm_ggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fmggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
  fmi_gtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fm_g_gtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fmg_gtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);

  addfm_g_gtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  addfm_ggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);

  fcgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);

  aigtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  ai_gtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);

  cai_gtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::setup_system()
{
  general_setup();

  // create combined map
  create_combined_dof_row_map();

  fluid_field()->use_block_matrix(true);

  // build ale system matrix in splitted system
  ale_field()->create_system_matrix(ale_field()->interface());


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
  vecSpaces.push_back(ale_field()->interface()->other_map());
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
  structureRHS->Update(1.0, *structure_field()->rhs(), 1.0, *AddStructRHS_, 0.0);

  Teuchos::RCP<Epetra_Vector> fluidRHS =
      Teuchos::rcp(new Epetra_Vector(*fluid_field()->discretization()->dof_row_map()));
  fluidRHS->Update(1.0, *fluid_field()->rhs(), 1.0, *AddFluidRHS_, 0.0);

  const double scale = fluid_field()->residual_scaling();

  setup_vector(f, structureRHS, fluidRHS, ale_field()->rhs(), ConstrRHS_, scale);

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
  structureRHS->Update(1.0, *structure_field()->rhs(), 1.0, *AddStructRHS_, 0.0);

  const double scale = fluid_field()->residual_scaling();

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockf = fluid_field()->block_system_matrix();

  Core::LinAlg::SparseMatrix& fig = blockf->matrix(0, 1);
  Core::LinAlg::SparseMatrix& fgg = blockf->matrix(1, 1);
  Core::LinAlg::SparseMatrix& fGg = blockf->matrix(3, 1);

  Teuchos::RCP<Epetra_Vector> fveln = fluid_field()->extract_interface_veln();
  double timescale = fluid_field()->time_scaling();

  Teuchos::RCP<Adapter::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<Adapter::FluidLung>(fluid_field());
  Teuchos::RCP<Epetra_Vector> rhs =
      Teuchos::rcp(new Epetra_Vector(*fluidfield->inner_split()->full_map()));
  Teuchos::RCP<Epetra_Vector> tmp = Teuchos::rcp(new Epetra_Vector(fig.row_map()));
  fig.Apply(*fveln, *tmp);
  fluidfield->inner_split()->insert_other_vector(tmp, rhs);
  tmp = Teuchos::rcp(new Epetra_Vector(fGg.row_map()));
  fGg.Apply(*fveln, *tmp);
  fluidfield->inner_split()->insert_cond_vector(tmp, rhs);
  rhs->Scale(timescale * dt());
  rhs = fluidfield->fsi_interface()->insert_other_vector(rhs);
  extractor().add_vector(*rhs, 1, f);

  rhs = Teuchos::rcp(new Epetra_Vector(fgg.row_map()));
  fgg.Apply(*fveln, *rhs);
  rhs->Scale(scale * timescale * dt());
  rhs = fluid_to_struct(rhs);
  rhs = structure_field()->interface()->insert_fsi_cond_vector(rhs);
  extractor().add_vector(*rhs, 0, f);

  //--------------------------------------------------------------------------------
  // constraint fluid
  //--------------------------------------------------------------------------------
  // split in two blocks according to inner and fsi structure dofs
  Teuchos::RCP<Epetra_Map> emptymap = Teuchos::rcp(
      new Epetra_Map(-1, 0, nullptr, 0, structure_field()->discretization()->get_comm()));
  Core::LinAlg::MapExtractor extractor_temp;
  extractor_temp.setup(*ConstrMap_, emptymap, ConstrMap_);

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> constrfluidblocks =
      ConstrFluidMatrix_->split<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *fluidfield->fsi_interface(), extractor_temp);
  constrfluidblocks->complete();

  Core::LinAlg::SparseMatrix& cfig = constrfluidblocks->matrix(0, 1);
  rhs = Teuchos::rcp(new Epetra_Vector(cfig.row_map()));
  cfig.Apply(*fveln, *rhs);
  rhs->Scale(timescale * dt());
  extractor().add_vector(*rhs, 3, f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::setup_system_matrix");

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W

  const Core::Adapter::Coupling& coupsf = structure_fluid_coupling();
  const Core::Adapter::Coupling& coupsa = structure_ale_coupling();
  const Core::Adapter::Coupling& coupfa = fluid_ale_coupling();

  /*----------------------------------------------------------------------*/

  double scale = fluid_field()->residual_scaling();
  double timescale = fluid_field()->time_scaling();

  /*----------------------------------------------------------------------*/
  // structure part

  Core::LinAlg::SparseMatrix s = *structure_field()->system_matrix();
  s.un_complete();
  s.add(AddStructConstrMatrix_->matrix(0, 0), false, 1.0, 1.0);

  mat.assign(0, 0, Core::LinAlg::View, s);

  /*----------------------------------------------------------------------*/
  // structure constraint part

  Core::LinAlg::SparseMatrix scon = AddStructConstrMatrix_->matrix(0, 1);
  mat.assign(0, 3, Core::LinAlg::View, scon);

  /*----------------------------------------------------------------------*/
  // fluid part

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockf = fluid_field()->block_system_matrix();

  Core::LinAlg::SparseMatrix& fii = blockf->matrix(0, 0);
  Core::LinAlg::SparseMatrix& fig = blockf->matrix(0, 1);
  Core::LinAlg::SparseMatrix& fiG = blockf->matrix(0, 3);
  Core::LinAlg::SparseMatrix& fgi = blockf->matrix(1, 0);
  Core::LinAlg::SparseMatrix& fgg = blockf->matrix(1, 1);
  Core::LinAlg::SparseMatrix& fgG = blockf->matrix(1, 3);
  Core::LinAlg::SparseMatrix& fGi = blockf->matrix(3, 0);
  Core::LinAlg::SparseMatrix& fGg = blockf->matrix(3, 1);
  Core::LinAlg::SparseMatrix& fGG = blockf->matrix(3, 3);

  // mat.Matrix(0,0).UnComplete();
  (*fggtransform_)(fgg, scale * timescale, Core::Adapter::CouplingSlaveConverter(coupsf),
      Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(0, 0), true, true);
  (*fgitransform_)(fgi, scale, Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(0, 1));
  (*fg_gtransform_)(
      fgG, scale, Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(0, 1), true);

  (*figtransform_)(blockf->full_row_map(), blockf->full_col_map(), fig, timescale,
      Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(1, 0));
  (*f_ggtransform_)(blockf->full_row_map(), blockf->full_col_map(), fGg, timescale,
      Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(1, 0), true, true);

  mat.matrix(1, 1).add(fii, false, 1.0, 0.0);
  mat.matrix(1, 1).add(fiG, false, 1.0, 1.0);
  mat.matrix(1, 1).add(fGi, false, 1.0, 1.0);
  mat.matrix(1, 1).add(fGG, false, 1.0, 1.0);

  Teuchos::RCP<Core::LinAlg::SparseMatrix> eye =
      Core::LinAlg::Eye(*fluid_field()->interface()->fsi_cond_map());
  mat.matrix(1, 1).add(*eye, false, 1.0, 1.0);

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> mmm = fluid_field()->shape_derivatives();

  if (mmm != Teuchos::null)
  {
    Core::LinAlg::SparseMatrix& fmii = mmm->matrix(0, 0);
    Core::LinAlg::SparseMatrix& fmig = mmm->matrix(0, 1);
    Core::LinAlg::SparseMatrix& fmiG = mmm->matrix(0, 3);
    Core::LinAlg::SparseMatrix& fmgi = mmm->matrix(1, 0);
    Core::LinAlg::SparseMatrix& fmgg = mmm->matrix(1, 1);
    Core::LinAlg::SparseMatrix& fmgG = mmm->matrix(1, 3);
    Core::LinAlg::SparseMatrix& fmGi = mmm->matrix(3, 0);
    Core::LinAlg::SparseMatrix& fmGg = mmm->matrix(3, 1);
    Core::LinAlg::SparseMatrix& fmGG = mmm->matrix(3, 3);

    // We cannot copy the pressure value. It is not used anyway. So no exact
    // match here.
    (*fmiitransform_)(mmm->full_row_map(), mmm->full_col_map(), fmii, 1.,
        Core::Adapter::CouplingMasterConverter(coupfa), mat.matrix(1, 2), false, false);
    (*fm_gitransform_)(mmm->full_row_map(), mmm->full_col_map(), fmGi, 1.,
        Core::Adapter::CouplingMasterConverter(coupfa), mat.matrix(1, 2), false, true);
    (*fmgitransform_)(fmgi, scale, Core::Adapter::CouplingSlaveConverter(coupsf),
        Core::Adapter::CouplingMasterConverter(coupfa), mat.matrix(0, 2), false, false);

    (*fmigtransform_)(mmm->full_row_map(), mmm->full_col_map(), fmig, 1.,
        Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(1, 0), true, true);
    (*fm_ggtransform_)(mmm->full_row_map(), mmm->full_col_map(), fmGg, 1.,
        Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(1, 0), true, true);

    (*fmggtransform_)(fmgg, scale, Core::Adapter::CouplingSlaveConverter(coupsf),
        Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(0, 0), true, true);

    (*fmi_gtransform_)(*coupfsout_->master_dof_map(), fmiG.col_map(), fmiG, 1.,
        Core::Adapter::CouplingMasterConverter(*coupfsout_), mat.matrix(1, 0), true, true);
    (*fm_g_gtransform_)(*coupfsout_->master_dof_map(), fmGG.col_map(), fmGG, 1.,
        Core::Adapter::CouplingMasterConverter(*coupfsout_), mat.matrix(1, 0), true, true);
    (*fmg_gtransform_)(fmgG, scale, Core::Adapter::CouplingSlaveConverter(coupsf),
        Core::Adapter::CouplingMasterConverter(*coupfsout_), mat.matrix(0, 0), true, true);

    Core::LinAlg::SparseMatrix& addfmGg = AddFluidShapeDerivMatrix_->matrix(3, 1);
    Core::LinAlg::SparseMatrix& addfmGG = AddFluidShapeDerivMatrix_->matrix(3, 3);

    (*addfm_g_gtransform_)(AddFluidShapeDerivMatrix_->full_row_map(),
        AddFluidShapeDerivMatrix_->full_col_map(), addfmGG, 1.,
        Core::Adapter::CouplingMasterConverter(*coupfsout_), mat.matrix(1, 0), true, true);

    (*addfm_ggtransform_)(AddFluidShapeDerivMatrix_->full_row_map(),
        AddFluidShapeDerivMatrix_->full_col_map(), addfmGg, 1.,
        Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(1, 0), true, true);
  }

  /*----------------------------------------------------------------------*/
  // fluid constraint part

  // split in two blocks according to inner and fsi dofs

  Teuchos::RCP<Adapter::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<Adapter::FluidLung>(fluid_field());

  Teuchos::RCP<Epetra_Map> emptymap =
      Teuchos::rcp(new Epetra_Map(-1, 0, nullptr, 0, fluid_field()->discretization()->get_comm()));
  Core::LinAlg::MapExtractor extractor;
  extractor.setup(*ConstrMap_, emptymap, ConstrMap_);

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> fluidconstrblocks =
      FluidConstrMatrix_->split<Core::LinAlg::DefaultBlockMatrixStrategy>(
          extractor, *fluidfield->fsi_interface());

  fluidconstrblocks->complete();

  Core::LinAlg::SparseMatrix& fcii = fluidconstrblocks->matrix(0, 0);
  Core::LinAlg::SparseMatrix& fcgi = fluidconstrblocks->matrix(1, 0);

  // fcii cannot be simply assigned here (in case of fsi amg which is default)
  // due to non-matching maps
  mat.matrix(1, 3).zero();
  mat.matrix(1, 3).add(fcii, false, 1.0, 0.0);

  // add interface part to structure block

  mat.matrix(0, 3).un_complete();
  (*fcgitransform_)(
      fcgi, scale, Core::Adapter::CouplingSlaveConverter(coupsf), mat.matrix(0, 3), true);

  /*----------------------------------------------------------------------*/
  // ale part

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> a = ale_field()->block_system_matrix();

  if (a == Teuchos::null) FOUR_C_THROW("expect ale block matrix");

  a->complete();

  Core::LinAlg::SparseMatrix& aii = a->matrix(0, 0);
  Core::LinAlg::SparseMatrix& aig = a->matrix(0, 1);
  Core::LinAlg::SparseMatrix& aiG = a->matrix(0, 3);

  mat.assign(2, 2, Core::LinAlg::View, aii);

  (*aigtransform_)(a->full_row_map(), a->full_col_map(), aig, 1.,
      Core::Adapter::CouplingSlaveConverter(coupsa), mat.matrix(2, 0));

  (*ai_gtransform_)(a->full_row_map(), a->full_col_map(), aiG, 1.,
      Core::Adapter::CouplingSlaveConverter(*coupsaout_), mat.matrix(2, 0), true, true);

  /*----------------------------------------------------------------------*/
  // constraint part -> structure

  mat.assign(3, 0, Core::LinAlg::View, AddStructConstrMatrix_->matrix(1, 0));

  /*----------------------------------------------------------------------*/
  // constraint part -> fluid

  // split in two blocks according to inner and fsi structure dofs

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> constrfluidblocks =
      ConstrFluidMatrix_->split<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *fluidfield->fsi_interface(), extractor);
  constrfluidblocks->complete();

  Core::LinAlg::SparseMatrix& cfii = constrfluidblocks->matrix(0, 0);

  // cfii cannot be simply assigned here (in case of fsi amg which is default)
  // due to non-matching maps
  mat.matrix(3, 1).zero();
  mat.matrix(3, 1).add(cfii, false, 1.0, 0.0);

  /*----------------------------------------------------------------------*/
  // constraint part -> "ale"

  Core::LinAlg::SparseMatrix& caiG = ConstrAleMatrix_->matrix(0, 3);
  (*cai_gtransform_)(*coupfsout_->master_dof_map(), caiG.col_map(), caiG, 1.0,
      Core::Adapter::CouplingMasterConverter(*coupfsout_), mat.matrix(3, 0), true, true);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.complete();
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

  Teuchos::RCP<Adapter::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<Adapter::FluidLung>(fluid_field());
  Teuchos::RCP<Epetra_Vector> fov = fluidfield->fsi_interface()->extract_other_vector(fv);
  fov = fluidfield->fsi_interface()->insert_other_vector(fov);
  Teuchos::RCP<Epetra_Vector> aov = ale_field()->interface()->extract_other_vector(av);

  if (fluidscale != 0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> fcv = fluid_field()->interface()->extract_fsi_cond_vector(fv);
    Teuchos::RCP<Epetra_Vector> modsv =
        structure_field()->interface()->insert_fsi_cond_vector(fluid_to_struct(fcv));
    modsv->Update(1.0, *sv, fluidscale);

    extractor().insert_vector(*modsv, 0, f);
  }
  else
  {
    extractor().insert_vector(*sv, 0, f);
  }

  extractor().insert_vector(*fov, 1, f);
  extractor().insert_vector(*aov, 2, f);
  extractor().insert_vector(*cv, 3, f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::LungMonolithicFluidSplit::extract_field_vectors");

  Teuchos::RCP<Adapter::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<Adapter::FluidLung>(fluid_field());

  // We have overlap at the interface. Thus we need the interface part of the
  // structure vector and append it to the fluid and ale vector. (With the
  // right translation.)

  sx = extractor().extract_vector(x, 0);
  Teuchos::RCP<const Epetra_Vector> scx =
      structure_field()->interface()->extract_fsi_cond_vector(sx);

  // process fluid unknowns

  Teuchos::RCP<const Epetra_Vector> fox = extractor().extract_vector(x, 1);
  fox = fluidfield->fsi_interface()->extract_other_vector(fox);
  Teuchos::RCP<Epetra_Vector> fcx = struct_to_fluid(scx);

  fluid_field()->displacement_to_velocity(fcx);

  Teuchos::RCP<Epetra_Vector> f = fluidfield->fsi_interface()->insert_other_vector(fox);
  fluid_field()->interface()->insert_fsi_cond_vector(fcx, f);
  fx = f;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = extractor().extract_vector(x, 2);
  Teuchos::RCP<Epetra_Vector> acx = struct_to_ale(scx);
  Teuchos::RCP<Epetra_Vector> a = ale_field()->interface()->insert_other_vector(aox);
  ale_field()->interface()->insert_vector(acx, 1, a);

  Teuchos::RCP<Epetra_Vector> scox =
      structure_field()->interface()->extract_lung_asi_cond_vector(sx);
  Teuchos::RCP<Epetra_Vector> acox = struct_to_ale_outflow(scox);
  ale_field()->interface()->insert_vector(acox, 3, a);

  ax = a;
}

FOUR_C_NAMESPACE_CLOSE
