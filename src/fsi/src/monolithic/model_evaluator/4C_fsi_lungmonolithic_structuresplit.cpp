/*----------------------------------------------------------------------*/
/*! \file
\brief Volume-coupled FSI (structure-split)


\level 3
*/
/*----------------------------------------------------------------------*/
#include "4C_fsi_lungmonolithic_structuresplit.hpp"

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
#include "4C_structure_aux.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

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
  intersectionmaps.push_back(structure_field()->get_dbc_map_extractor()->cond_map());
  intersectionmaps.push_back(structure_field()->interface()->fsi_cond_map());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      Core::LinAlg::MultiMapExtractor::intersect_maps(intersectionmaps);

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
void FSI::LungMonolithicStructureSplit::setup_system()
{
  general_setup();
  // create combined map
  create_combined_dof_row_map();

  fluid_field()->use_block_matrix(false);

  // build ale system matrix in splitted system
  ale_field()->create_system_matrix(ale_field()->interface());

  //-----------------------------------------------------------------------------
  // create block system matrix
  //-----------------------------------------------------------------------------
  create_system_matrix(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::create_combined_dof_row_map()
{
  const Teuchos::RCP<Adapter::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<Adapter::StructureLung>(structure_field());

  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(structfield->fsi_interface()->other_map());
  vecSpaces.push_back(fluid_field()->dof_row_map());
  // remaining (not coupled) dofs of ale field
  vecSpaces.push_back(ale_field()->interface()->Map(0));
  // additional volume constraints
  vecSpaces.push_back(ConstrMap_);

  if (vecSpaces[0]->NumGlobalElements() == 0)
    FOUR_C_THROW("No inner structural equations. Splitting not possible. Panic.");

  set_dof_row_maps(vecSpaces);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::setup_rhs_residual(Epetra_Vector& f)
{
  Teuchos::RCP<Epetra_Vector> structureRHS =
      Teuchos::rcp(new Epetra_Vector(*structure_field()->discretization()->dof_row_map()));
  structureRHS->Update(1.0, *structure_field()->rhs(), 1.0, *AddStructRHS_, 0.0);

  const double fluidscale = fluid_field()->residual_scaling();

  Teuchos::RCP<Epetra_Vector> fluidRHS =
      Teuchos::rcp(new Epetra_Vector(*fluid_field()->discretization()->dof_row_map()));
  fluidRHS->Update(1.0, *fluid_field()->rhs(), 1.0, *AddFluidRHS_, 0.0);

  setup_vector(f, structureRHS, fluidRHS, ale_field()->rhs(), ConstrRHS_, fluidscale);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::setup_rhs_lambda(Epetra_Vector& f)
{
  // ToDo: We still need to implement this.

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::setup_rhs_firstiter(Epetra_Vector& f)
{
  const Teuchos::RCP<Adapter::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<Adapter::StructureLung>(structure_field());
  const Teuchos::RCP<Adapter::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<Adapter::FluidLung>(fluid_field());

  const double fluidscale = fluid_field()->residual_scaling();

  // additional rhs term for ALE equations
  // -dt Aig u(n)
  //
  //    1/dt Delta d(n+1) = theta Delta u(n+1) + u(n)
  //
  // And we are concerned with the u(n) part here.

  //--------------------------------------------------------------------------------
  // ale
  //--------------------------------------------------------------------------------
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> a = ale_field()->block_system_matrix();
  if (a == Teuchos::null) FOUR_C_THROW("expect ale block matrix");

  Core::LinAlg::SparseMatrix& aig = a->matrix(0, 1);

  Teuchos::RCP<Epetra_Vector> fveln = fluid_field()->extract_interface_veln();
  Teuchos::RCP<Epetra_Vector> sveln = fluid_to_struct(fveln);
  Teuchos::RCP<Epetra_Vector> aveln = struct_to_ale(sveln);
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(aig.row_map()));
  aig.Apply(*aveln, *rhs);

  rhs->Scale(-1. * dt());

  extractor().add_vector(*rhs, 2, f);

  //--------------------------------------------------------------------------------
  // structure
  //--------------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> veln = structure_field()->interface()->insert_fsi_cond_vector(sveln);
  rhs = Teuchos::rcp(new Epetra_Vector(veln->Map()));

  Teuchos::RCP<Core::LinAlg::SparseMatrix> s = structure_field()->system_matrix();
  s->Apply(*veln, *rhs);

  Teuchos::RCP<Epetra_Vector> addrhs = Teuchos::rcp(new Epetra_Vector(veln->Map()));
  AddStructConstrMatrix_->matrix(0, 0).Apply(*veln, *addrhs);

  rhs->Update(1.0, *addrhs, 1.0);
  rhs->Scale(-1. * dt());

  veln = structfield->fsi_interface()->extract_other_vector(rhs);
  extractor().add_vector(*veln, 0, f);

  veln = structure_field()->interface()->extract_fsi_cond_vector(rhs);
  veln = fluid_field()->interface()->insert_fsi_cond_vector(struct_to_fluid(veln));

  veln->Scale(1. / fluidscale);
  extractor().add_vector(*veln, 1, f);

  //--------------------------------------------------------------------------------
  // constraint structure
  //--------------------------------------------------------------------------------
  // split in two blocks according to inner and fsi structure dofs
  Teuchos::RCP<Epetra_Map> emptymap = Teuchos::rcp(
      new Epetra_Map(-1, 0, nullptr, 0, structure_field()->discretization()->get_comm()));
  Core::LinAlg::MapExtractor extractor_temp;
  extractor_temp.setup(*ConstrMap_, emptymap, ConstrMap_);

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> constrstructblocks =
      AddStructConstrMatrix_->matrix(1, 0).split<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *structfield->fsi_interface(), extractor_temp);
  constrstructblocks->complete();

  Core::LinAlg::SparseMatrix& csig = constrstructblocks->matrix(0, 1);

  rhs = Teuchos::rcp(new Epetra_Vector(csig.row_map()));
  csig.Apply(*sveln, *rhs);
  rhs->Scale(-1. * dt());
  extractor().add_vector(*rhs, 3, f);

  //--------------------------------------------------------------------------------
  // constraint ale
  //--------------------------------------------------------------------------------
  Core::LinAlg::SparseMatrix& caig = ConstrAleMatrix_->matrix(0, 1);
  caig.Apply(*fveln, *rhs);
  rhs->Scale(-1. * dt());
  extractor().add_vector(*rhs, 3, f);

  //--------------------------------------------------------------------------------
  // shape derivatives
  //--------------------------------------------------------------------------------
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> mmm = fluidfield->shape_derivatives();

  if (mmm != Teuchos::null)
  {
    Core::LinAlg::SparseMatrix& fmig = mmm->matrix(0, 1);
    Core::LinAlg::SparseMatrix& fmgg = mmm->matrix(1, 1);
    Core::LinAlg::SparseMatrix& fmGg = mmm->matrix(3, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmig.row_map()));
    fmig.Apply(*fveln, *rhs);
    veln = fluid_field()->interface()->insert_vector(rhs, 0);

    rhs = Teuchos::rcp(new Epetra_Vector(fmgg.row_map()));
    fmgg.Apply(*fveln, *rhs);
    fluid_field()->interface()->insert_vector(rhs, 1, veln);

    rhs = Teuchos::rcp(new Epetra_Vector(fmGg.row_map()));
    fmGg.Apply(*fveln, *rhs);
    fluid_field()->interface()->insert_vector(rhs, 3, veln);

    veln->Scale(-1. * dt());
    extractor().add_vector(*veln, 1, f);

    Core::LinAlg::SparseMatrix& afmgg = AddFluidShapeDerivMatrix_->matrix(1, 1);
    Core::LinAlg::SparseMatrix& afmGg = AddFluidShapeDerivMatrix_->matrix(3, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(afmgg.row_map()));
    afmgg.Apply(*fveln, *rhs);
    veln = fluid_field()->interface()->insert_vector(rhs, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(afmGg.row_map()));
    afmGg.Apply(*fveln, *rhs);
    fluid_field()->interface()->insert_vector(rhs, 3, veln);

    veln->Scale(-1. * dt());
    extractor().add_vector(*veln, 1, f);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::setup_system_matrix(
    Core::LinAlg::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::setup_system_matrix");

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W

  const Core::Adapter::Coupling& coupsf = structure_fluid_coupling();

  Teuchos::RCP<Core::LinAlg::SparseMatrix> f = fluid_field()->system_matrix();

  double scale = fluid_field()->residual_scaling();
  double timescale = fluid_field()->time_scaling();

  /*----------------------------------------------------------------------*/
  // fluid part
  mat.assign(1, 1, Core::LinAlg::View, *f);

  // fluid linearization with respect to mesh motion block
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> mmm = fluid_field()->shape_derivatives();
  const Core::Adapter::Coupling& coupfa = fluid_ale_coupling();

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

    Core::LinAlg::SparseMatrix& addfmGg = AddFluidShapeDerivMatrix_->matrix(3, 1);
    Core::LinAlg::SparseMatrix& addfmGG = AddFluidShapeDerivMatrix_->matrix(3, 3);

    mat.matrix(1, 1).add(fmig, false, 1. / timescale, 1.0);
    mat.matrix(1, 1).add(fmgg, false, 1. / timescale, 1.0);
    mat.matrix(1, 1).add(fmGg, false, 1. / timescale, 1.0);

    fmiitransform_(mmm->full_row_map(), mmm->full_col_map(), fmii, 1.,
        Core::Adapter::CouplingMasterConverter(coupfa), mat.matrix(1, 2), false, false);
    fmgitransform_(mmm->full_row_map(), mmm->full_col_map(), fmgi, 1.,
        Core::Adapter::CouplingMasterConverter(coupfa), mat.matrix(1, 2), false, true);
    fm_gitransform_(mmm->full_row_map(), mmm->full_col_map(), fmGi, 1.,
        Core::Adapter::CouplingMasterConverter(coupfa), mat.matrix(1, 2), false, true);

    fmg_gtransform_(mmm->full_row_map(), mmm->full_col_map(), fmgG, 1.,
        Core::Adapter::CouplingMasterConverter(*coupfsout_), mat.matrix(1, 0), true, false);
    fmi_gtransform_(mmm->full_row_map(), mmm->full_col_map(), fmiG, 1.,
        Core::Adapter::CouplingMasterConverter(*coupfsout_), mat.matrix(1, 0), true, true);
    fm_g_gtransform_(mmm->full_row_map(), mmm->full_col_map(), fmGG, 1.,
        Core::Adapter::CouplingMasterConverter(*coupfsout_), mat.matrix(1, 0), true, true);

    mat.matrix(1, 1).add(addfmGg, false, 1. / timescale, 1.0);

    addfm_g_gtransform_(AddFluidShapeDerivMatrix_->full_row_map(),
        AddFluidShapeDerivMatrix_->full_col_map(), addfmGG, 1.,
        Core::Adapter::CouplingMasterConverter(*coupfsout_), mat.matrix(1, 0), true, true);
  }

  /*----------------------------------------------------------------------*/
  // structure and additional structure part

  const Teuchos::RCP<Adapter::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<Adapter::StructureLung>(structure_field());

  Core::LinAlg::SparseMatrix& s = *structure_field()->system_matrix();
  s.add(AddStructConstrMatrix_->matrix(0, 0), false, 1.0, 1.0);

  siitransform_(s, *structfield->fsi_interface()->Map(0), *structfield->fsi_interface()->Map(0),
      1.0, nullptr, nullptr, mat.matrix(0, 0), true, false);

  Core::Adapter::CouplingMasterConverter converter(coupsf);
  sigtransform_(s, *structfield->fsi_interface()->Map(0), *structfield->fsi_interface()->Map(1),
      1. / timescale, nullptr, &converter, mat.matrix(0, 1), true, false);

  sgitransform_(s, *structfield->fsi_interface()->Map(1), *structfield->fsi_interface()->Map(0),
      1. / scale, &converter, nullptr, mat.matrix(1, 0), true, true);

  sggtransform_(s, *structfield->fsi_interface()->Map(1), *structfield->fsi_interface()->Map(1),
      1. / (scale * timescale), &converter, &converter, mat.matrix(1, 1), true, true);

  /*----------------------------------------------------------------------*/
  // fluid constraint part
  // add into new matrix instead of view -> avoids matrix UnComplete() operations in the
  // second and subsequent time steps
  mat.matrix(1, 3).add(*FluidConstrMatrix_, false, 1., 0.);

  /*----------------------------------------------------------------------*/
  // structure constraint part

  // split in two blocks according to inner and fsi structure dofs
  Core::LinAlg::MapExtractor extractor;

  sciitransform_(AddStructConstrMatrix_->matrix(0, 1), *structfield->fsi_interface()->Map(0),
      *ConstrMap_, 1., nullptr, nullptr, mat.matrix(0, 3), true, false);

  // add interface part to fluid block
  scgitransform_(AddStructConstrMatrix_->matrix(0, 1), *structfield->fsi_interface()->Map(1),
      *ConstrMap_, 1. / scale, &converter, nullptr, mat.matrix(1, 3), true, true);

  /*----------------------------------------------------------------------*/
  // ale part
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> a = ale_field()->block_system_matrix();

  if (a == Teuchos::null) FOUR_C_THROW("expect ale block matrix");

  a->complete();

  Core::LinAlg::SparseMatrix& aii = a->matrix(0, 0);
  Core::LinAlg::SparseMatrix& aig = a->matrix(0, 1);
  Core::LinAlg::SparseMatrix& aiG = a->matrix(0, 3);

  ai_gtransform_(a->full_row_map(), a->full_col_map(), aiG, 1.,
      Core::Adapter::CouplingSlaveConverter(*coupsaout_), mat.matrix(2, 0));
  aigtransform_(a->full_row_map(), a->full_col_map(), aig, 1. / timescale,
      Core::Adapter::CouplingSlaveConverter(*icoupfa_), mat.matrix(2, 1));
  mat.assign(2, 2, Core::LinAlg::View, aii);

  /*----------------------------------------------------------------------*/
  // constraint part -> fluid
  // add into new matrix instead of view -> avoids matrix UnComplete() operations in the
  // second and subsequent time steps

  mat.matrix(3, 1).add(*ConstrFluidMatrix_, false, 1., 0.);

  /*----------------------------------------------------------------------*/
  // constraint part -> structure
  // split in two blocks according to inner and fsi structure dofs

  csiitransform_(AddStructConstrMatrix_->matrix(1, 0), *ConstrMap_,
      *structfield->fsi_interface()->Map(0), 1., nullptr, nullptr, mat.matrix(3, 0), true, false);
  csigtransform_(AddStructConstrMatrix_->matrix(1, 0), *ConstrMap_,
      *structfield->fsi_interface()->Map(1), 1. / timescale, nullptr, &converter, mat.matrix(3, 1),
      true, true);

  /*----------------------------------------------------------------------*/
  // constraint part -> ale

  Core::LinAlg::SparseMatrix& caiG = ConstrAleMatrix_->matrix(0, 3);
  cai_gtransform_(*coupfsout_->master_dof_map(), caiG.col_map(), caiG, 1.0,
      Core::Adapter::CouplingMasterConverter(*coupfsout_), mat.matrix(3, 0), true, true);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::initial_guess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::initial_guess");

  Teuchos::RCP<Epetra_Vector> ConstraintInitialGuess =
      Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));

  setup_vector(*ig, structure_field()->initial_guess(), fluid_field()->initial_guess(),
      ale_field()->initial_guess(), ConstraintInitialGuess, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::setup_vector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av, Teuchos::RCP<const Epetra_Vector> cv, double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  const Teuchos::RCP<Adapter::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<Adapter::StructureLung>(structure_field());

  Teuchos::RCP<Epetra_Vector> sov = structfield->fsi_interface()->extract_other_vector(sv);

  Teuchos::RCP<Epetra_Vector> aov = ale_field()->interface()->extract_vector(av, 0);

  if (fluidscale != 0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> scv = structure_field()->interface()->extract_fsi_cond_vector(sv);
    Teuchos::RCP<Epetra_Vector> modfv =
        fluid_field()->interface()->insert_fsi_cond_vector(struct_to_fluid(scv));
    modfv->Update(1.0, *fv, 1. / fluidscale);

    extractor().insert_vector(*modfv, 1, f);
  }
  else
  {
    extractor().insert_vector(*fv, 1, f);
  }

  extractor().insert_vector(*sov, 0, f);
  extractor().insert_vector(*aov, 2, f);
  extractor().insert_vector(*cv, 3, f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  const Teuchos::RCP<Adapter::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<Adapter::StructureLung>(structure_field());

  fx = extractor().extract_vector(x, 1);

  // process structure unknowns

  Teuchos::RCP<Epetra_Vector> fcx = fluid_field()->interface()->extract_fsi_cond_vector(fx);
  fluid_field()->velocity_to_displacement(fcx);
  Teuchos::RCP<const Epetra_Vector> sox = extractor().extract_vector(x, 0);
  Teuchos::RCP<Epetra_Vector> scx = fluid_to_struct(fcx);

  Teuchos::RCP<Epetra_Vector> s = structfield->fsi_interface()->insert_other_vector(sox);
  structure_field()->interface()->insert_fsi_cond_vector(scx, s);
  sx = s;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = extractor().extract_vector(x, 2);
  Teuchos::RCP<Epetra_Vector> acx = struct_to_ale(scx);
  Teuchos::RCP<Epetra_Vector> a = ale_field()->interface()->insert_vector(aox, 0);
  ale_field()->interface()->insert_vector(acx, 1, a);

  Teuchos::RCP<Epetra_Vector> scox =
      structure_field()->interface()->extract_lung_asi_cond_vector(sx);
  Teuchos::RCP<Epetra_Vector> acox = struct_to_ale_outflow(scox);
  ale_field()->interface()->insert_vector(acox, 3, a);

  ax = a;
}

FOUR_C_NAMESPACE_CLOSE
