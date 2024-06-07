/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problem with constraints

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_fsi_constrmonolithic_structuresplit.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid.hpp"
#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_adapter_str_structure.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_constraint_manager.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_io_control.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_structure_aux.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

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
  intersectionmaps.push_back(structure_field()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(structure_field()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      Core::LinAlg::MultiMapExtractor::IntersectMaps(intersectionmaps);

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

  sggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
  sgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  sigtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  aigtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  scgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  csigtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::SetupSystem()
{
  GeneralSetup();

  // create combined map
  create_combined_dof_row_map();

  fluid_field()->use_block_matrix(false);

  // Use splitted structure matrix
  structure_field()->use_block_matrix();

  Teuchos::RCP<Epetra_Map> emptymap =
      Teuchos::rcp(new Epetra_Map(-1, 0, nullptr, 0, structure_field()->discretization()->Comm()));
  Teuchos::RCP<Core::LinAlg::MapExtractor> extractor;
  extractor->Setup(*conman_->GetConstraintMap(), emptymap, conman_->GetConstraintMap());
  conman_->use_block_matrix(extractor, structure_field()->Interface());
  scon_t_ =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *structure_field()->Interface(), *extractor, 81, false, true));

  // build ale system matrix in splitted system
  ale_field()->create_system_matrix(ale_field()->Interface());

  create_system_matrix(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::create_combined_dof_row_map()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(structure_field()->Interface()->OtherMap());
  vecSpaces.push_back(fluid_field()->dof_row_map());
  vecSpaces.push_back(ale_field()->Interface()->OtherMap());

  vecSpaces.push_back(conman_->GetConstraintMap());

  if (vecSpaces[0]->NumGlobalElements() == 0)
    FOUR_C_THROW("No inner structural equations. Splitting not possible. Panic.");

  set_dof_row_maps(vecSpaces);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::setup_rhs_residual(Epetra_Vector& f)
{
  const double scale = fluid_field()->residual_scaling();

  setup_vector(f, structure_field()->RHS(), fluid_field()->RHS(), ale_field()->RHS(),
      conman_->GetError(), scale);

  // add additional ale residual
  extractor().AddVector(*aleresidual_, 2, f);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::setup_rhs_lambda(Epetra_Vector& f)
{
  // ToDo: We still need to implement this.

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::setup_rhs_firstiter(Epetra_Vector& f)
{
  // additional rhs term for ALE equations
  // -dt Aig u(n)
  //
  //    1/dt Delta d(n+1) = theta Delta u(n+1) + u(n)
  //
  // And we are concerned with the u(n) part here.

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> a = ale_field()->BlockSystemMatrix();
  if (a == Teuchos::null) FOUR_C_THROW("expect ale block matrix");

  const Core::LinAlg::SparseMatrix& aig = a->Matrix(0, 1);

  Teuchos::RCP<const Epetra_Vector> fveln = fluid_field()->extract_interface_veln();
  Teuchos::RCP<const Epetra_Vector> sveln = fluid_to_struct(fveln);
  Teuchos::RCP<const Epetra_Vector> aveln = struct_to_ale(sveln);
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap()));
  aig.Apply(*aveln, *rhs);

  rhs->Scale(-1. * Dt());

  extractor().AddVector(*rhs, 2, f);

  // structure
  Teuchos::RCP<Epetra_Vector> veln = structure_field()->Interface()->InsertFSICondVector(sveln);
  rhs = Teuchos::rcp(new Epetra_Vector(veln->Map()));

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> s = structure_field()->block_system_matrix();
  s->Apply(*veln, *rhs);

  rhs->Scale(-1. * Dt());

  veln = structure_field()->Interface()->ExtractOtherVector(rhs);
  extractor().AddVector(*veln, 0, f);

  veln = structure_field()->Interface()->ExtractFSICondVector(rhs);
  veln = fluid_field()->Interface()->InsertFSICondVector(struct_to_fluid(veln));

  const double scale = fluid_field()->residual_scaling();

  veln->Scale(1. / scale);

  extractor().AddVector(*veln, 1, f);

  // shape derivatives
  Teuchos::RCP<const Core::LinAlg::BlockSparseMatrixBase> mmm = fluid_field()->ShapeDerivatives();
  if (mmm != Teuchos::null)
  {
    const Core::LinAlg::SparseMatrix& fmig = mmm->Matrix(0, 1);
    const Core::LinAlg::SparseMatrix& fmgg = mmm->Matrix(1, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap()));
    fmig.Apply(*fveln, *rhs);
    veln = fluid_field()->Interface()->InsertOtherVector(rhs);

    rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap()));
    fmgg.Apply(*fveln, *rhs);
    fluid_field()->Interface()->InsertFSICondVector(rhs, veln);

    veln->Scale(-1. * Dt());

    extractor().AddVector(*veln, 1, f);
  }

  //--------------------------------------------------------------------------------
  // constraint
  //--------------------------------------------------------------------------------
  Core::LinAlg::SparseOperator& tmp = *conman_->GetConstrMatrix();
  Core::LinAlg::BlockSparseMatrixBase& scon =
      dynamic_cast<Core::LinAlg::BlockSparseMatrixBase&>(tmp);
  for (int rowblock = 0; rowblock < scon.Rows(); ++rowblock)
  {
    for (int colblock = 0; colblock < scon.Cols(); ++colblock)
    {
      scon_t_->Matrix(colblock, rowblock).Add(scon.Matrix(rowblock, colblock), true, 1.0, 0.0);
    }
  }
  scon_t_->Complete();

  Core::LinAlg::SparseMatrix& csig = scon_t_->Matrix(0, 1);

  rhs = Teuchos::rcp(new Epetra_Vector(csig.RowMap()));
  csig.Apply(*sveln, *rhs);
  rhs->Scale(-1. * Dt());
  extractor().AddVector(*rhs, 3, f);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::setup_system_matrix(
    Core::LinAlg::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::setup_system_matrix");

  // extract Jacobian matrices and put them into composite system
  // matrix W

  const Core::Adapter::Coupling& coupsf = structure_fluid_coupling();
  // const Adapter::Coupling& coupsa = structure_ale_coupling();

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> s = structure_field()->block_system_matrix();
  if (s == Teuchos::null) FOUR_C_THROW("expect structure block matrix");
  Teuchos::RCP<Core::LinAlg::SparseMatrix> f = fluid_field()->SystemMatrix();
  if (f == Teuchos::null) FOUR_C_THROW("expect fluid matrix");
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> a = ale_field()->BlockSystemMatrix();
  if (a == Teuchos::null) FOUR_C_THROW("expect ale block matrix");

  Core::LinAlg::SparseMatrix& aii = a->Matrix(0, 0);
  Core::LinAlg::SparseMatrix& aig = a->Matrix(0, 1);

  /*----------------------------------------------------------------------*/

  double scale = fluid_field()->residual_scaling();
  double timescale = fluid_field()->TimeScaling();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  f->UnComplete();

  mat.Assign(0, 0, Core::LinAlg::View, s->Matrix(0, 0));

  (*sigtransform_)(s->FullRowMap(), s->FullColMap(), s->Matrix(0, 1), 1. / timescale,
      Core::Adapter::CouplingMasterConverter(coupsf), mat.Matrix(0, 1));
  (*sggtransform_)(s->Matrix(1, 1), 1. / (scale * timescale),
      Core::Adapter::CouplingMasterConverter(coupsf),
      Core::Adapter::CouplingMasterConverter(coupsf), *f, true, true);
  (*sgitransform_)(s->Matrix(1, 0), 1. / scale, Core::Adapter::CouplingMasterConverter(coupsf),
      mat.Matrix(1, 0));

  mat.Assign(1, 1, Core::LinAlg::View, *f);

  (*aigtransform_)(a->FullRowMap(), a->FullColMap(), aig, 1. / timescale,
      Core::Adapter::CouplingSlaveConverter(*icoupfa_), mat.Matrix(2, 1));
  mat.Assign(2, 2, Core::LinAlg::View, aii);

  /*----------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> mmm = fluid_field()->ShapeDerivatives();
  if (mmm != Teuchos::null)
  {
    Core::LinAlg::SparseMatrix& fmii = mmm->Matrix(0, 0);
    Core::LinAlg::SparseMatrix& fmig = mmm->Matrix(0, 1);
    Core::LinAlg::SparseMatrix& fmgi = mmm->Matrix(1, 0);
    Core::LinAlg::SparseMatrix& fmgg = mmm->Matrix(1, 1);

    mat.Matrix(1, 1).Add(fmgg, false, 1. / timescale, 1.0);
    mat.Matrix(1, 1).Add(fmig, false, 1. / timescale, 1.0);

    const Core::Adapter::Coupling& coupfa = fluid_ale_coupling();

    (*fmgitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmgi, 1.,
        Core::Adapter::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, false);

    (*fmiitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmii, 1.,
        Core::Adapter::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, true);
  }


  /*----------------------------------------------------------------------*/
  // structure constraint part

  Core::LinAlg::SparseOperator& tmp = *conman_->GetConstrMatrix();
  Core::LinAlg::BlockSparseMatrixBase& scon =
      dynamic_cast<Core::LinAlg::BlockSparseMatrixBase&>(tmp);
  for (int rowblock = 0; rowblock < scon.Rows(); ++rowblock)
  {
    for (int colblock = 0; colblock < scon.Cols(); ++colblock)
    {
      scon_t_->Matrix(colblock, rowblock).Add(scon.Matrix(rowblock, colblock), true, 1.0, 0.0);
    }
  }
  scon_t_->Complete();


  scon.Complete();

  mat.Assign(0, 3, Core::LinAlg::View, scon.Matrix(0, 0));

  (*scgitransform_)(scon.Matrix(1, 0), 1. / scale, Core::Adapter::CouplingMasterConverter(coupsf),
      mat.Matrix(1, 3));

  mat.Assign(3, 0, Core::LinAlg::View, scon_t_->Matrix(0, 0));

  (*csigtransform_)(*coupsf.MasterDofMap(), scon_t_->Matrix(0, 1).ColMap(), scon_t_->Matrix(0, 1),
      1. / timescale, Core::Adapter::CouplingMasterConverter(coupsf), mat.Matrix(3, 1), true);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();

  // Finally, take care of Dirichlet boundary conditions
  mat.ApplyDirichlet(*(dbcmaps_->CondMap()), true);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::initial_guess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::initial_guess");

  Teuchos::RCP<Epetra_Vector> ConstraintInitialGuess =
      Teuchos::rcp(new Epetra_Vector(*(conman_->GetConstraintMap()), true));

  setup_vector(*ig, structure_field()->initial_guess(), fluid_field()->initial_guess(),
      ale_field()->initial_guess(), ConstraintInitialGuess, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::setup_vector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av, Teuchos::RCP<const Epetra_Vector> cv, double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<Epetra_Vector> sov = structure_field()->Interface()->ExtractOtherVector(sv);
  Teuchos::RCP<Epetra_Vector> aov = ale_field()->Interface()->ExtractOtherVector(av);

  if (fluidscale != 0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> scv = structure_field()->Interface()->ExtractFSICondVector(sv);
    Teuchos::RCP<Epetra_Vector> modfv =
        fluid_field()->Interface()->InsertFSICondVector(struct_to_fluid(scv));
    modfv->Update(1.0, *fv, 1. / fluidscale);

    extractor().InsertVector(*modfv, 1, f);
  }
  else
  {
    extractor().InsertVector(*fv, 1, f);
  }

  extractor().InsertVector(*sov, 0, f);
  extractor().InsertVector(*aov, 2, f);

  Epetra_Vector modcv = *cv;
  modcv.Scale(-1.0);
  extractor().InsertVector(modcv, 3, f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::ConstrMonolithicStructureSplit::extract_field_vectors");

  fx = extractor().ExtractVector(x, 1);

  // process structure unknowns

  Teuchos::RCP<Epetra_Vector> fcx = fluid_field()->Interface()->ExtractFSICondVector(fx);
  fluid_field()->velocity_to_displacement(fcx);
  Teuchos::RCP<const Epetra_Vector> sox = extractor().ExtractVector(x, 0);
  Teuchos::RCP<Epetra_Vector> scx = fluid_to_struct(fcx);

  Teuchos::RCP<Epetra_Vector> s = structure_field()->Interface()->InsertOtherVector(sox);
  structure_field()->Interface()->InsertFSICondVector(scx, s);
  sx = s;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = extractor().ExtractVector(x, 2);
  Teuchos::RCP<Epetra_Vector> acx = struct_to_ale(scx);

  Teuchos::RCP<Epetra_Vector> a = ale_field()->Interface()->InsertOtherVector(aox);
  ale_field()->Interface()->InsertFSICondVector(acx, a);

  ax = a;
}

FOUR_C_NAMESPACE_CLOSE
