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
  intersectionmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(StructureField()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      CORE::LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

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

  sggtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowColTransform);
  sgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  sigtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  aigtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  scgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  csigtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::SetupSystem()
{
  GeneralSetup();

  // create combined map
  create_combined_dof_row_map();

  fluid_field()->UseBlockMatrix(false);

  // Use splitted structure matrix
  StructureField()->UseBlockMatrix();

  Teuchos::RCP<Epetra_Map> emptymap =
      Teuchos::rcp(new Epetra_Map(-1, 0, nullptr, 0, StructureField()->Discretization()->Comm()));
  Teuchos::RCP<CORE::LINALG::MapExtractor> extractor;
  extractor->Setup(*conman_->GetConstraintMap(), emptymap, conman_->GetConstraintMap());
  conman_->UseBlockMatrix(extractor, StructureField()->Interface());
  scon_t_ =
      Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *StructureField()->Interface(), *extractor, 81, false, true));

  // build ale system matrix in splitted system
  ale_field()->CreateSystemMatrix(ale_field()->Interface());

  CreateSystemMatrix(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::create_combined_dof_row_map()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(StructureField()->Interface()->OtherMap());
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
  const double scale = fluid_field()->ResidualScaling();

  setup_vector(f, StructureField()->RHS(), fluid_field()->RHS(), ale_field()->RHS(),
      conman_->GetError(), scale);

  // add additional ale residual
  Extractor().AddVector(*aleresidual_, 2, f);

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

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> a = ale_field()->BlockSystemMatrix();
  if (a == Teuchos::null) FOUR_C_THROW("expect ale block matrix");

  const CORE::LINALG::SparseMatrix& aig = a->Matrix(0, 1);

  Teuchos::RCP<const Epetra_Vector> fveln = fluid_field()->extract_interface_veln();
  Teuchos::RCP<const Epetra_Vector> sveln = FluidToStruct(fveln);
  Teuchos::RCP<const Epetra_Vector> aveln = StructToAle(sveln);
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap()));
  aig.Apply(*aveln, *rhs);

  rhs->Scale(-1. * Dt());

  Extractor().AddVector(*rhs, 2, f);

  // structure
  Teuchos::RCP<Epetra_Vector> veln = StructureField()->Interface()->InsertFSICondVector(sveln);
  rhs = Teuchos::rcp(new Epetra_Vector(veln->Map()));

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> s = StructureField()->BlockSystemMatrix();
  s->Apply(*veln, *rhs);

  rhs->Scale(-1. * Dt());

  veln = StructureField()->Interface()->ExtractOtherVector(rhs);
  Extractor().AddVector(*veln, 0, f);

  veln = StructureField()->Interface()->ExtractFSICondVector(rhs);
  veln = fluid_field()->Interface()->InsertFSICondVector(StructToFluid(veln));

  const double scale = fluid_field()->ResidualScaling();

  veln->Scale(1. / scale);

  Extractor().AddVector(*veln, 1, f);

  // shape derivatives
  Teuchos::RCP<const CORE::LINALG::BlockSparseMatrixBase> mmm = fluid_field()->ShapeDerivatives();
  if (mmm != Teuchos::null)
  {
    const CORE::LINALG::SparseMatrix& fmig = mmm->Matrix(0, 1);
    const CORE::LINALG::SparseMatrix& fmgg = mmm->Matrix(1, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap()));
    fmig.Apply(*fveln, *rhs);
    veln = fluid_field()->Interface()->InsertOtherVector(rhs);

    rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap()));
    fmgg.Apply(*fveln, *rhs);
    fluid_field()->Interface()->InsertFSICondVector(rhs, veln);

    veln->Scale(-1. * Dt());

    Extractor().AddVector(*veln, 1, f);
  }

  //--------------------------------------------------------------------------------
  // constraint
  //--------------------------------------------------------------------------------
  CORE::LINALG::SparseOperator& tmp = *conman_->GetConstrMatrix();
  CORE::LINALG::BlockSparseMatrixBase& scon =
      dynamic_cast<CORE::LINALG::BlockSparseMatrixBase&>(tmp);
  for (int rowblock = 0; rowblock < scon.Rows(); ++rowblock)
  {
    for (int colblock = 0; colblock < scon.Cols(); ++colblock)
    {
      scon_t_->Matrix(colblock, rowblock).Add(scon.Matrix(rowblock, colblock), true, 1.0, 0.0);
    }
  }
  scon_t_->Complete();

  CORE::LINALG::SparseMatrix& csig = scon_t_->Matrix(0, 1);

  rhs = Teuchos::rcp(new Epetra_Vector(csig.RowMap()));
  csig.Apply(*sveln, *rhs);
  rhs->Scale(-1. * Dt());
  Extractor().AddVector(*rhs, 3, f);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::setup_system_matrix(
    CORE::LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::setup_system_matrix");

  // extract Jacobian matrices and put them into composite system
  // matrix W

  const CORE::ADAPTER::Coupling& coupsf = structure_fluid_coupling();
  // const ADAPTER::Coupling& coupsa = structure_ale_coupling();

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> s = StructureField()->BlockSystemMatrix();
  if (s == Teuchos::null) FOUR_C_THROW("expect structure block matrix");
  Teuchos::RCP<CORE::LINALG::SparseMatrix> f = fluid_field()->SystemMatrix();
  if (f == Teuchos::null) FOUR_C_THROW("expect fluid matrix");
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> a = ale_field()->BlockSystemMatrix();
  if (a == Teuchos::null) FOUR_C_THROW("expect ale block matrix");

  CORE::LINALG::SparseMatrix& aii = a->Matrix(0, 0);
  CORE::LINALG::SparseMatrix& aig = a->Matrix(0, 1);

  /*----------------------------------------------------------------------*/

  double scale = fluid_field()->ResidualScaling();
  double timescale = fluid_field()->TimeScaling();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  f->UnComplete();

  mat.Assign(0, 0, CORE::LINALG::View, s->Matrix(0, 0));

  (*sigtransform_)(s->FullRowMap(), s->FullColMap(), s->Matrix(0, 1), 1. / timescale,
      CORE::ADAPTER::CouplingMasterConverter(coupsf), mat.Matrix(0, 1));
  (*sggtransform_)(s->Matrix(1, 1), 1. / (scale * timescale),
      CORE::ADAPTER::CouplingMasterConverter(coupsf),
      CORE::ADAPTER::CouplingMasterConverter(coupsf), *f, true, true);
  (*sgitransform_)(s->Matrix(1, 0), 1. / scale, CORE::ADAPTER::CouplingMasterConverter(coupsf),
      mat.Matrix(1, 0));

  mat.Assign(1, 1, CORE::LINALG::View, *f);

  (*aigtransform_)(a->FullRowMap(), a->FullColMap(), aig, 1. / timescale,
      CORE::ADAPTER::CouplingSlaveConverter(*icoupfa_), mat.Matrix(2, 1));
  mat.Assign(2, 2, CORE::LINALG::View, aii);

  /*----------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> mmm = fluid_field()->ShapeDerivatives();
  if (mmm != Teuchos::null)
  {
    CORE::LINALG::SparseMatrix& fmii = mmm->Matrix(0, 0);
    CORE::LINALG::SparseMatrix& fmig = mmm->Matrix(0, 1);
    CORE::LINALG::SparseMatrix& fmgi = mmm->Matrix(1, 0);
    CORE::LINALG::SparseMatrix& fmgg = mmm->Matrix(1, 1);

    mat.Matrix(1, 1).Add(fmgg, false, 1. / timescale, 1.0);
    mat.Matrix(1, 1).Add(fmig, false, 1. / timescale, 1.0);

    const CORE::ADAPTER::Coupling& coupfa = FluidAleCoupling();

    (*fmgitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmgi, 1.,
        CORE::ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, false);

    (*fmiitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmii, 1.,
        CORE::ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, true);
  }


  /*----------------------------------------------------------------------*/
  // structure constraint part

  CORE::LINALG::SparseOperator& tmp = *conman_->GetConstrMatrix();
  CORE::LINALG::BlockSparseMatrixBase& scon =
      dynamic_cast<CORE::LINALG::BlockSparseMatrixBase&>(tmp);
  for (int rowblock = 0; rowblock < scon.Rows(); ++rowblock)
  {
    for (int colblock = 0; colblock < scon.Cols(); ++colblock)
    {
      scon_t_->Matrix(colblock, rowblock).Add(scon.Matrix(rowblock, colblock), true, 1.0, 0.0);
    }
  }
  scon_t_->Complete();


  scon.Complete();

  mat.Assign(0, 3, CORE::LINALG::View, scon.Matrix(0, 0));

  (*scgitransform_)(scon.Matrix(1, 0), 1. / scale, CORE::ADAPTER::CouplingMasterConverter(coupsf),
      mat.Matrix(1, 3));

  mat.Assign(3, 0, CORE::LINALG::View, scon_t_->Matrix(0, 0));

  (*csigtransform_)(*coupsf.MasterDofMap(), scon_t_->Matrix(0, 1).ColMap(), scon_t_->Matrix(0, 1),
      1. / timescale, CORE::ADAPTER::CouplingMasterConverter(coupsf), mat.Matrix(3, 1), true);

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

  setup_vector(*ig, StructureField()->initial_guess(), fluid_field()->initial_guess(),
      ale_field()->initial_guess(), ConstraintInitialGuess, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::ConstrMonolithicStructureSplit::setup_vector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av, Teuchos::RCP<const Epetra_Vector> cv, double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<Epetra_Vector> sov = StructureField()->Interface()->ExtractOtherVector(sv);
  Teuchos::RCP<Epetra_Vector> aov = ale_field()->Interface()->ExtractOtherVector(av);

  if (fluidscale != 0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> scv = StructureField()->Interface()->ExtractFSICondVector(sv);
    Teuchos::RCP<Epetra_Vector> modfv =
        fluid_field()->Interface()->InsertFSICondVector(StructToFluid(scv));
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
void FSI::ConstrMonolithicStructureSplit::extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::ConstrMonolithicStructureSplit::extract_field_vectors");

  fx = Extractor().ExtractVector(x, 1);

  // process structure unknowns

  Teuchos::RCP<Epetra_Vector> fcx = fluid_field()->Interface()->ExtractFSICondVector(fx);
  fluid_field()->velocity_to_displacement(fcx);
  Teuchos::RCP<const Epetra_Vector> sox = Extractor().ExtractVector(x, 0);
  Teuchos::RCP<Epetra_Vector> scx = FluidToStruct(fcx);

  Teuchos::RCP<Epetra_Vector> s = StructureField()->Interface()->InsertOtherVector(sox);
  StructureField()->Interface()->InsertFSICondVector(scx, s);
  sx = s;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x, 2);
  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);

  Teuchos::RCP<Epetra_Vector> a = ale_field()->Interface()->InsertOtherVector(aox);
  ale_field()->Interface()->InsertFSICondVector(acx, a);

  ax = a;
}

FOUR_C_NAMESPACE_CLOSE
