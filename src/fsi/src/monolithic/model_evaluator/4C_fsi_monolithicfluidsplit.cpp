/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problem with matching grids using a monolithic scheme
with condensed fluid interface velocities


\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_fsi_monolithicfluidsplit.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_adapter_str_structure.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_constraint_manager.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_fsi_monolithic_linearsystem.hpp"
#include "4C_fsi_nox_group.hpp"
#include "4C_fsi_overlapprec.hpp"
#include "4C_fsi_overlapprec_fsiamg.hpp"
#include "4C_fsi_statustest.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_structure_aux.hpp"

#include <NOX_Epetra_LinearSystem.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicFluidSplit::MonolithicFluidSplit(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : BlockMonolithic(comm, timeparams),
      lambda_(Teuchos::null),
      lambdaold_(Teuchos::null),
      energysum_(0.0)
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
      Core::LinAlg::MultiMapExtractor::IntersectMaps(intersectionmaps);

  // Check whether the intersection is empty
  if (intersectionmap->NumGlobalElements() != 0)
  {
    //    std::cout << "Slave interface nodes with Dirichlet boundary condition (input file
    //    numbering):" << std::endl; for (int i=0; i <
    //    (int)fluid_field()->discretization()->NumMyRowNodes(); i++)
    //    {
    //      // get all nodes and add them
    //      int gid = fluid_field()->discretization()->NodeRowMap()->GID(i);
    //
    //      // do only nodes that I have in my discretization
    //      if (!fluid_field()->discretization()->NodeColMap()->MyGID(gid)) continue;
    //      Core::Nodes::Node* node = fluid_field()->discretization()->gNode(gid);
    //      if (!node) FOUR_C_THROW("Cannot find node with gid %",gid);
    //
    //      std::vector<int> nodedofs = fluid_field()->discretization()->Dof(node);
    //
    //      for (int j=0; j < (int)nodedofs.size(); j++)
    //      {
    //        for (int k=0; k < (int)intersectionmap->NumGlobalElements(); k++)
    //        {
    //          if (nodedofs[j] == intersectionmap->GID(k))
    //          {
    //            std::cout << gid+1 << std::endl;
    //            k = (int)intersectionmap->GID(k);
    //            j = (int)nodedofs.size();
    //          }
    //        }
    //      }
    //    }


    // It is not allowed, that slave DOFs at the interface hold a Dirichlet
    // boundary condition. Thus --> Error message

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

    FOUR_C_THROW(errormsg.str());
  }
  // ---------------------------------------------------------------------------

  fggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
  fgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  figtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  aigtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
  fmggtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);

  // Recovery of Lagrange multiplier happens on fluid field
  lambda_ = Teuchos::rcp(new Epetra_Vector(*fluid_field()->Interface()->FSICondMap(), true));
  lambdaold_ = Teuchos::rcp(new Epetra_Vector(*fluid_field()->Interface()->FSICondMap(), true));
  fmgiprev_ = Teuchos::null;
  fmgicur_ = Teuchos::null;
  fmggprev_ = Teuchos::null;
  fmggcur_ = Teuchos::null;
  fgiprev_ = Teuchos::null;
  fgicur_ = Teuchos::null;
  fggprev_ = Teuchos::null;
  fggcur_ = Teuchos::null;

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // check whether allocation was successful
  if (fggtransform_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'fggtransform_' failed.");
  }
  if (fgitransform_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'fgitransform_' failed.");
  }
  if (figtransform_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'figtransform_' failed.");
  }
  if (aigtransform_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'aigtransform_' failed.");
  }
  if (fmiitransform_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'fmiitransform_' failed.");
  }
  if (fmgitransform_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'fmgitransform_' failed.");
  }
  if (fmggtransform_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'fmggtransform_' failed.");
  }
  if (lambda_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'lambda_' failed.");
  }
  if (lambdaold_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'lambdaold_' failed.");
  }
#endif

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::SetupSystem()
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  linearsolverstrategy_ =
      Core::UTILS::IntegralValue<Inpar::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

  set_default_parameters(fsidyn, nox_parameter_list());

  // call SetupSystem in base class
  FSI::Monolithic::SetupSystem();

  // create combined map
  create_combined_dof_row_map();

  /*----------------------------------------------------------------------*/
  // Switch fluid to interface split block matrix
  fluid_field()->use_block_matrix(true);

  // build ale system matrix in splitted system
  ale_field()->create_system_matrix(ale_field()->Interface());

  aleresidual_ = Teuchos::rcp(new Epetra_Vector(*ale_field()->Interface()->OtherMap()));

  // ---------------------------------------------------------------------------
  // Build the global Dirichlet map extractor
  setup_dbc_map_extractor();
  // ---------------------------------------------------------------------------

  // enable debugging
  if (Core::UTILS::IntegralValue<int>(fsidyn, "DEBUGOUTPUT") & 2)
  {
    pcdbg_ = Teuchos::rcp(new UTILS::MonolithicDebugWriter(*this));
  }

  // create the system matrix
  create_system_matrix();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::create_system_matrix()
{
  FSI::BlockMonolithic::create_system_matrix(systemmatrix_, false);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::create_combined_dof_row_map()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(structure_field()->dof_row_map());
  vecSpaces.push_back(fluid_field()->dof_row_map());
  vecSpaces.push_back(ale_field()->Interface()->OtherMap());

  if (vecSpaces[1]->NumGlobalElements() == 0)
    FOUR_C_THROW("No inner fluid equations. Splitting not possible. Panic.");

  set_dof_row_maps(vecSpaces);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::setup_dbc_map_extractor()
{
  // Dirichlet maps for structure and fluid do not intersect with interface map.
  // ALE Dirichlet map might intersect with interface map, but ALE interface DOFs
  // are not part of the final system of equations. Hence, we just need the
  // intersection of inner ALE DOFs with Dirichlet ALE DOFs.
  std::vector<Teuchos::RCP<const Epetra_Map>> aleintersectionmaps;
  aleintersectionmaps.push_back(ale_field()->GetDBCMapExtractor()->CondMap());
  aleintersectionmaps.push_back(ale_field()->Interface()->OtherMap());
  Teuchos::RCP<Epetra_Map> aleintersectionmap =
      Core::LinAlg::MultiMapExtractor::IntersectMaps(aleintersectionmaps);

  // Merge Dirichlet maps of structure, fluid and ALE to global FSI Dirichlet map
  std::vector<Teuchos::RCP<const Epetra_Map>> dbcmaps;
  dbcmaps.push_back(structure_field()->GetDBCMapExtractor()->CondMap());
  dbcmaps.push_back(fluid_field()->GetDBCMapExtractor()->CondMap());
  dbcmaps.push_back(aleintersectionmap);
  Teuchos::RCP<const Epetra_Map> dbcmap = Core::LinAlg::MultiMapExtractor::MergeMaps(dbcmaps);

  // Finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new Core::LinAlg::MapExtractor(*dof_row_map(), dbcmap, true));
  if (dbcmaps_ == Teuchos::null)
  {
    FOUR_C_THROW("Creation of FSI Dirichlet map extractor failed.");
  }
  // ---------------------------------------------------------------------------

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> FSI::MonolithicFluidSplit::SystemMatrix() const
{
  return systemmatrix_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::setup_rhs_residual(Epetra_Vector& f)
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = structure_field()->TimIntParam();
  const double ftiparam = fluid_field()->TimIntParam();

  // some scaling factors for fluid
  const double fluidscale = fluid_field()->residual_scaling();

  // get single field residuals
  Teuchos::RCP<const Epetra_Vector> sv = Teuchos::rcp(new Epetra_Vector(*structure_field()->RHS()));
  Teuchos::RCP<const Epetra_Vector> fv = Teuchos::rcp(new Epetra_Vector(*fluid_field()->RHS()));
  Teuchos::RCP<const Epetra_Vector> av = Teuchos::rcp(new Epetra_Vector(*ale_field()->RHS()));

  // extract only inner DOFs from fluid (=slave) and ALE field
  Teuchos::RCP<Epetra_Vector> fov = fluid_field()->Interface()->ExtractOtherVector(fv);
  fov = fluid_field()->Interface()->InsertOtherVector(fov);
  Teuchos::RCP<const Epetra_Vector> aov = ale_field()->Interface()->ExtractOtherVector(av);

  // add fluid interface values to structure vector
  Teuchos::RCP<Epetra_Vector> fcv = fluid_field()->Interface()->ExtractFSICondVector(fv);
  Teuchos::RCP<Epetra_Vector> modsv =
      structure_field()->Interface()->InsertFSICondVector(fluid_to_struct(fcv));
  modsv->Update(1.0, *sv, (1.0 - stiparam) / (1.0 - ftiparam) * fluidscale);

  if (structure_field()->get_stc_algo() == Inpar::STR::stc_currsym)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> stcmat = structure_field()->get_stc_mat();
    stcmat->Multiply(true, *modsv, *modsv);
  }

  // put the single field residuals together
  FSI::Monolithic::combine_field_vectors(f, modsv, fov, aov);

  // add additional ale residual
  extractor().AddVector(*aleresidual_, 2, f);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::setup_rhs_lambda(Epetra_Vector& f)
{
  if (lambdaold_ != Teuchos::null)
  {
    // get time integration parameters of structure and fluid time integrators
    // to enable consistent time integration among the fields
    const double stiparam = structure_field()->TimIntParam();
    const double ftiparam = fluid_field()->TimIntParam();

    // project Lagrange multiplier field onto the master interface DOFs and consider temporal
    // scaling
    Teuchos::RCP<Epetra_Vector> lambdafull =
        structure_field()->Interface()->InsertFSICondVector(fluid_to_struct(lambdaold_));
    lambdafull->Scale(stiparam - (ftiparam * (1.0 - stiparam)) / (1.0 - ftiparam));

    // add Lagrange multiplier
    extractor().AddVector(*lambdafull, 0, f);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::setup_rhs_firstiter(Epetra_Vector& f)
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = structure_field()->TimIntParam();
  const double ftiparam = fluid_field()->TimIntParam();

  // some scaling factors for fluid
  const double timescale = fluid_field()->TimeScaling();
  const double scale = fluid_field()->residual_scaling();

  // old interface velocity of fluid field
  const Teuchos::RCP<const Epetra_Vector> fveln = fluid_field()->extract_interface_veln();

  // get fluid matrix
  const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockf =
      fluid_field()->BlockSystemMatrix();

  // get fluid shape derivatives matrix
  const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> mmm = fluid_field()->ShapeDerivatives();

  // get ale matrix
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blocka = ale_field()->BlockSystemMatrix();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (blockf == Teuchos::null)
  {
    FOUR_C_THROW("Expected Teuchos::rcp to fluid block matrix.");
  }
  if (blocka == Teuchos::null)
  {
    FOUR_C_THROW("Expected Teuchos::rcp to ale block matrix.");
  }
#endif

  // extract fluid and ale submatrices
  const Core::LinAlg::SparseMatrix& fig = blockf->Matrix(0, 1);  // F_{I\Gamma}
  const Core::LinAlg::SparseMatrix& fgg = blockf->Matrix(1, 1);  // F_{\Gamma\Gamma}
  const Core::LinAlg::SparseMatrix& aig = blocka->Matrix(0, 1);  // A_{I\Gamma}

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;

  // Different contributions/terms to the rhs are separated by the following comment line
  // ---------- structural interface DOFs
  /* The following terms are added to the structural interface DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  + (1-stiparam)/(1-ftiparam) * dt / tau * F_{\Gamma\Gamma} * u^{n}_{\Gamma}
   *
   * (2)  - (1-stiparam)/(1-ftiparam) / tau * F_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
   *
   * (3)  - (1-stiparam)/(1-ftiparam) * F^{G}_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
   *
   * Remarks on all terms:
   * +  tau: time scaling factor for interface time integration (tau =
   * 1/fluid_field()->TimeScaling())
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(fgg.RangeMap(), true));

  fgg.Apply(*fveln, *rhs);

  rhs->Scale(scale * (1. - stiparam) / (1. - ftiparam) * Dt() * timescale);

  rhs = fluid_to_struct(rhs);
  rhs = structure_field()->Interface()->InsertFSICondVector(rhs);

  if (structure_field()->get_stc_algo() == Inpar::STR::stc_currsym)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> stcmat = structure_field()->get_stc_mat();
    stcmat->Multiply(true, *rhs, *rhs);
  }

  extractor().AddVector(*rhs, 0, f);
  // ----------end of term 1

  // ----------addressing term 2:
  rhs = Teuchos::rcp(new Epetra_Vector(fgg.RangeMap(), true));

  fgg.Apply(*struct_to_fluid(ddgpred_), *rhs);

  rhs->Scale(-scale * (1. - stiparam) / (1. - ftiparam) * timescale);
  rhs = structure_field()->Interface()->InsertFSICondVector(fluid_to_struct(rhs));

  extractor().AddVector(*rhs, 0, f);
  // ----------end of term 2

  // ----------addressing term 3:
  if (mmm != Teuchos::null)
  {
    // extract F^{G}_{\Gamma\Gamma}
    const Core::LinAlg::SparseMatrix& fmgg = mmm->Matrix(1, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RangeMap(), true));

    fmgg.Apply(*struct_to_fluid(ddgpred_), *rhs);

    rhs->Scale(-(1. - stiparam) / (1. - ftiparam));
    rhs = structure_field()->Interface()->InsertFSICondVector(fluid_to_struct(rhs));

    extractor().AddVector(*rhs, 0, f);
  }
  // ----------end of term 3
  // ----------end of structural interface DOFs

  // ---------- inner fluid DOFs
  /* The following terms are added to the inner fluid DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  + dt / tau * F_{I \Gamma} * u^{n}_{\Gamma}
   *
   * (2)  - 1 / tau F_{I \Gamma} * \Delta d_{\Gamma,p}
   *
   * (3)  - F^{G}_{I \Gamma} * \Delta d_{\Gamma,p}
   *
   * Remarks on all terms:
   * +  tau: time scaling factor for interface time integration (tau =
   * 1/fluid_field()->TimeScaling())
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(fig.RangeMap(), true));

  fig.Apply(*fveln, *rhs);

  rhs->Scale(Dt() * timescale);

  rhs = fluid_field()->Interface()->InsertOtherVector(rhs);

  extractor().AddVector(*rhs, 1, f);
  // ----------end of term 1

  // ----------addressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(fig.RangeMap(), true));

  fig.Apply(*struct_to_fluid(ddgpred_), *rhs);

  rhs->Scale(-timescale);

  rhs = fluid_field()->Interface()->InsertOtherVector(rhs);

  extractor().AddVector(*rhs, 1, f);
  // ----------end of term 2

  // ----------addressing term 3
  if (mmm != Teuchos::null)
  {
    // extract F^{G}_{I \Gamma}
    const Core::LinAlg::SparseMatrix& fmig = mmm->Matrix(0, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmig.RangeMap(), true));

    fmig.Apply(*struct_to_fluid(ddgpred_), *rhs);

    rhs->Scale(-1.);

    rhs = fluid_field()->Interface()->InsertOtherVector(rhs);

    extractor().AddVector(*rhs, 1, f);
  }
  // ----------end of term 3
  // ----------end of inner fluid DOFs

  // ---------- inner ale DOFs
  /* The following terms are added to the inner ale DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - A_{I \Gamma} * \Delta d_{\Gamma,p}
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(aig.RangeMap(), true));

  aig.Apply(*struct_to_ale(ddgpred_), *rhs);
  rhs->Scale(-1.0);

  extractor().AddVector(*rhs, 2, f);
  // ----------end of term 1
  // ---------- end of inner ale DOFs

  // Reset quantities of previous iteration step since they still store values from the last time
  // step
  ddginc_ = Core::LinAlg::CreateVector(*structure_field()->Interface()->FSICondMap(), true);
  duiinc_ = Core::LinAlg::CreateVector(*fluid_field()->Interface()->OtherMap(), true);
  ddialeinc_ = Core::LinAlg::CreateVector(*ale_field()->Interface()->OtherMap(), true);
  soliprev_ = Teuchos::null;
  solgprev_ = Teuchos::null;
  fgicur_ = Teuchos::null;
  fggcur_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::setup_system_matrix");

  const Core::Adapter::Coupling& coupsf = structure_fluid_coupling();
  const Core::Adapter::Coupling& coupsa = structure_ale_coupling();
  const Core::Adapter::Coupling& coupfa = fluid_ale_coupling();

  // get info about STC feature
  Inpar::STR::StcScale stcalgo = structure_field()->get_stc_algo();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> stcmat = Teuchos::null;
  // if STC is to be used, get STC matrix from structure field
  if (stcalgo != Inpar::STR::stc_none) stcmat = structure_field()->get_stc_mat();

  // get single field block matrices
  Teuchos::RCP<Core::LinAlg::SparseMatrix> s =
      structure_field()->system_matrix();  // can't be 'const' --> is modified by STC
  const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> f = fluid_field()->BlockSystemMatrix();
  const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> a = ale_field()->BlockSystemMatrix();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // check whether allocation was successful
  if (s == Teuchos::null)
  {
    FOUR_C_THROW("expect structure block matrix");
  }
  if (f == Teuchos::null)
  {
    FOUR_C_THROW("expect fluid block matrix");
  }
  if (a == Teuchos::null)
  {
    FOUR_C_THROW("expect ale block matrix");
  }
#endif

  // extract submatrices
  Core::LinAlg::SparseMatrix& fii = f->Matrix(0, 0);
  Core::LinAlg::SparseMatrix& fig = f->Matrix(0, 1);
  Core::LinAlg::SparseMatrix& fgi = f->Matrix(1, 0);
  Core::LinAlg::SparseMatrix& fgg = f->Matrix(1, 1);
  Core::LinAlg::SparseMatrix& aii = a->Matrix(0, 0);
  Core::LinAlg::SparseMatrix& aig = a->Matrix(0, 1);

  // scaling factors for fluid
  const double scale = fluid_field()->residual_scaling();
  const double timescale = fluid_field()->TimeScaling();

  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = structure_field()->TimIntParam();
  const double ftiparam = fluid_field()->TimIntParam();

  /*----------------------------------------------------------------------*/

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  s->UnComplete();

  (*fggtransform_)(fgg, (1.0 - stiparam) / (1.0 - ftiparam) * scale * timescale,
      Core::Adapter::CouplingSlaveConverter(coupsf), Core::Adapter::CouplingSlaveConverter(coupsf),
      *s, true, true);

  Teuchos::RCP<Core::LinAlg::SparseMatrix> lfgi =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(s->RowMap(), 81, false));
  (*fgitransform_)(fgi, (1.0 - stiparam) / (1.0 - ftiparam) * scale,
      Core::Adapter::CouplingSlaveConverter(coupsf), *lfgi);

  lfgi->Complete(fgi.DomainMap(), s->RangeMap());

  if (stcalgo == Inpar::STR::stc_currsym)
    lfgi = Core::LinAlg::MLMultiply(*stcmat, true, *lfgi, false, true, true, true);

  mat.Matrix(0, 1).UnComplete();
  mat.Matrix(0, 1).Add(*lfgi, false, 1., 0.0);

  if (stcalgo == Inpar::STR::stc_none)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> lfig =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(fig.RowMap(), 81, false));
    (*figtransform_)(f->FullRowMap(), f->FullColMap(), fig, timescale,
        Core::Adapter::CouplingSlaveConverter(coupsf), mat.Matrix(1, 0));
  }
  else
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> lfig =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(fig.RowMap(), 81, false));
    (*figtransform_)(f->FullRowMap(), f->FullColMap(), fig, timescale,
        Core::Adapter::CouplingSlaveConverter(coupsf), *lfig);

    lfig->Complete(s->DomainMap(), fig.RangeMap());

    lfig = Core::LinAlg::MLMultiply(*lfig, false, *stcmat, false, false, false, true);

    mat.Matrix(1, 0).UnComplete();
    mat.Matrix(1, 0).Add(*lfig, false, 1., 0.0);
  }

  mat.Matrix(1, 1).Add(fii, false, 1., 0.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> eye =
      Core::LinAlg::Eye(*fluid_field()->Interface()->FSICondMap());
  mat.Matrix(1, 1).Add(*eye, false, 1., 1.0);

  if (stcalgo == Inpar::STR::stc_none)
  {
    (*aigtransform_)(a->FullRowMap(), a->FullColMap(), aig, 1.,
        Core::Adapter::CouplingSlaveConverter(coupsa), mat.Matrix(2, 0));
  }
  else
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> laig =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(aii.RowMap(), 81, false));
    (*aigtransform_)(a->FullRowMap(), a->FullColMap(), aig, 1.,
        Core::Adapter::CouplingSlaveConverter(coupsa), *laig);

    laig->Complete(s->DomainMap(), laig->RangeMap());

    if (stcalgo != Inpar::STR::stc_none)
    {
      laig = Core::LinAlg::MLMultiply(*laig, false, *stcmat, false, false, false, true);
    }

    mat.Assign(2, 0, Core::LinAlg::View, *laig);
  }

  mat.Assign(2, 2, Core::LinAlg::View, aii);

  /*----------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> mmm = fluid_field()->ShapeDerivatives();
  if (mmm != Teuchos::null)
  {
    Core::LinAlg::SparseMatrix& fmii = mmm->Matrix(0, 0);
    Core::LinAlg::SparseMatrix& fmgi = mmm->Matrix(1, 0);

    Core::LinAlg::SparseMatrix& fmig = mmm->Matrix(0, 1);
    Core::LinAlg::SparseMatrix& fmgg = mmm->Matrix(1, 1);

    // reuse transform objects to add shape derivative matrices to structural blocks

    if (stcalgo == Inpar::STR::stc_none)
    {
      (*figtransform_)(f->FullRowMap(), f->FullColMap(), fmig, 1.,
          Core::Adapter::CouplingSlaveConverter(coupsf), mat.Matrix(1, 0), false, true);
    }
    else
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> lfmig =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(fmig.RowMap(), 81, false));
      (*figtransform_)(f->FullRowMap(), f->FullColMap(), fmig, 1.,
          Core::Adapter::CouplingSlaveConverter(coupsf), *lfmig, false, true);


      lfmig->Complete(s->DomainMap(), fmig.RangeMap());

      if (stcalgo != Inpar::STR::stc_none)
      {
        lfmig = Core::LinAlg::MLMultiply(*lfmig, false, *stcmat, false, false, false, true);
      }

      mat.Matrix(1, 0).Add(*lfmig, false, 1.0, 1.0);
    }

    (*fmggtransform_)(fmgg, (1.0 - stiparam) / (1.0 - ftiparam) * scale,
        Core::Adapter::CouplingSlaveConverter(coupsf),
        Core::Adapter::CouplingSlaveConverter(coupsf), *s, false, true);

    // We cannot copy the pressure value. It is not used anyway. So no exact
    // match here.
    (*fmiitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmii, 1.,
        Core::Adapter::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false);

    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> lfmgi =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(s->RowMap(), 81, false));
      (*fmgitransform_)(fmgi, (1.0 - stiparam) / (1.0 - ftiparam) * scale,
          Core::Adapter::CouplingSlaveConverter(coupsf),
          Core::Adapter::CouplingMasterConverter(coupfa), *lfmgi, false, false);

      lfmgi->Complete(aii.DomainMap(), s->RangeMap());

      if (stcalgo == Inpar::STR::stc_currsym)
        lfmgi = Core::LinAlg::MLMultiply(*stcmat, true, *lfmgi, false, true, true, true);

      mat.Matrix(0, 2).UnComplete();
      mat.Matrix(0, 2).Add(*lfmgi, false, 1., 0.0);
    }
  }

  s->Complete();

  if (stcalgo == Inpar::STR::stc_none)
  {
    s->UnComplete();
  }
  else  // apply STC matrix on block (0,0) if STC is used
  {
    s = Core::LinAlg::MLMultiply(*s, false, *stcmat, false, true, true, true);
    if (stcalgo == Inpar::STR::stc_currsym)
      s = Core::LinAlg::MLMultiply(*stcmat, true, *s, false, true, true, false);
  }

  // finally assign structure block
  mat.Matrix(0, 0).Assign(Core::LinAlg::View, *s);

  // done. make sure all blocks are filled.
  mat.Complete();

  // Finally, take care of Dirichlet boundary conditions
  mat.ApplyDirichlet(*(dbcmaps_->CondMap()), true);

  // store parts of fluid matrix to know them in the next iteration as previous iteration matrices
  fgiprev_ = fgicur_;
  fggprev_ = fggcur_;
  fgicur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(f->Matrix(1, 0)));
  fggcur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(f->Matrix(1, 1)));

  // store parts of fluid shape derivative matrix to know them in the next iteration as previous
  // iteration matrices
  fmgiprev_ = fmgicur_;
  fmggprev_ = fmggcur_;
  if (mmm != Teuchos::null)
  {
    fmgicur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(mmm->Matrix(1, 0)));
    fmggcur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(mmm->Matrix(1, 1)));
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::scale_system(
    Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)Core::UTILS::IntegralValue<int>(fsimono, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    // do scaling of structure rows
    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if (A->LeftScale(*srowsum_) or A->RightScale(*scolsum_) or
        mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->RightScale(*scolsum_))
      FOUR_C_THROW("structure scaling failed");

    // do scaling of ale rows
    A = mat.Matrix(2, 2).EpetraMatrix();
    arowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    acolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*arowsum_);
    A->InvColSums(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2, 1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1, 2).EpetraMatrix()->RightScale(*acolsum_))
      FOUR_C_THROW("ale scaling failed");

    // do scaling of structure and ale rhs vectors
    Teuchos::RCP<Epetra_Vector> sx = extractor().ExtractVector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = extractor().ExtractVector(b, 2);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ax->Multiply(1.0, *arowsum_, *ax, 0.0)) FOUR_C_THROW("ale scaling failed");

    extractor().InsertVector(*sx, 0, b);
    extractor().InsertVector(*ax, 2, b);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::unscale_solution(
    Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)Core::UTILS::IntegralValue<int>(fsimono, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> sy = extractor().ExtractVector(x, 0);
    Teuchos::RCP<Epetra_Vector> ay = extractor().ExtractVector(x, 2);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ay->Multiply(1.0, *acolsum_, *ay, 0.0)) FOUR_C_THROW("ale scaling failed");

    // get info about STC feature and unscale solution if necessary
    Inpar::STR::StcScale stcalgo = structure_field()->get_stc_algo();
    if (stcalgo != Inpar::STR::stc_none)
    {
      structure_field()->get_stc_mat()->Multiply(false, *sy, *sy);
    }

    extractor().InsertVector(*sy, 0, x);
    extractor().InsertVector(*ay, 2, x);

    Teuchos::RCP<Epetra_Vector> sx = extractor().ExtractVector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = extractor().ExtractVector(b, 2);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ax->ReciprocalMultiply(1.0, *arowsum_, *ax, 0.0)) FOUR_C_THROW("ale scaling failed");

    if (stcalgo != Inpar::STR::stc_none)
    {
      structure_field()->get_stc_mat()->Multiply(false, *sx, *sx);
    }

    extractor().InsertVector(*sx, 0, b);
    extractor().InsertVector(*ax, 2, b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if (A->LeftScale(*srowsum_) or A->RightScale(*scolsum_) or
        mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->RightScale(*scolsum_))
      FOUR_C_THROW("structure scaling failed");

    A = mat.Matrix(2, 2).EpetraMatrix();
    arowsum_->Reciprocal(*arowsum_);
    acolsum_->Reciprocal(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2, 1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1, 2).EpetraMatrix()->RightScale(*acolsum_))
      FOUR_C_THROW("ale scaling failed");
  }

  // very simple hack just to see the linear solution

  Epetra_Vector r(b.Map());
  mat.Apply(x, r);
  r.Update(1., b, 1.);

  Teuchos::RCP<Epetra_Vector> sr = extractor().ExtractVector(r, 0);
  Teuchos::RCP<Epetra_Vector> fr = extractor().ExtractVector(r, 1);
  Teuchos::RCP<Epetra_Vector> ar = extractor().ExtractVector(r, 2);

  // increment additional ale residual
  aleresidual_->Update(-1., *ar, 0.);

  std::ios_base::fmtflags flags = utils()->out().flags();

  double n, ns, nf, na;
  r.Norm2(&n);
  sr->Norm2(&ns);
  fr->Norm2(&nf);
  ar->Norm2(&na);
  if (verbosity_ == Inpar::FSI::verbosity_full)
  {
    utils()->out() << std::scientific << "\nlinear solver quality:\n"
                   << "L_2-norms:\n"
                   << "   |r|=" << n << "   |rs|=" << ns << "   |rf|=" << nf << "   |ra|=" << na
                   << "\n";
  }
  r.NormInf(&n);
  sr->NormInf(&ns);
  fr->NormInf(&nf);
  ar->NormInf(&na);
  if (verbosity_ == Inpar::FSI::verbosity_full)
  {
    utils()->out() << "L_inf-norms:\n"
                   << "   |r|=" << n << "   |rs|=" << ns << "   |rf|=" << nf << "   |ra|=" << na
                   << "\n";
  }
  utils()->out().flags(flags);

  if (structure_field()->get_stc_algo() != Inpar::STR::stc_none)
    structure_field()->system_matrix()->Reset();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Combo> FSI::MonolithicFluidSplit::create_status_test(
    Teuchos::ParameterList& nlParams, Teuchos::RCP<::NOX::Epetra::Group> grp)
{
  // --------------------------------------------------------------------
  // Setup the test framework
  // --------------------------------------------------------------------
  // Create the top-level test combo
  Teuchos::RCP<::NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::OR));

  // Create test combo for convergence of residuals and iterative increments
  Teuchos::RCP<::NOX::StatusTest::Combo> converged =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::AND));

  // Create some other plausibility tests
  Teuchos::RCP<::NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new ::NOX::StatusTest::MaxIters(nlParams.get<int>("Max Iterations")));
  Teuchos::RCP<::NOX::StatusTest::FiniteValue> fv =
      Teuchos::rcp(new ::NOX::StatusTest::FiniteValue);

  // Add single tests to the top-level test combo
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Start filling the 'converged' combo here
  // require one solve
  converged->addStatusTest(Teuchos::rcp(new NOX::FSI::MinIters(1)));


  // --------------------------------------------------------------------
  // setup tests for structural displacement field
  // --------------------------------------------------------------------
  // create ::NOX::StatusTest::Combo for structural displacement field
  Teuchos::RCP<::NOX::StatusTest::Combo> structcombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp_L2 = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "DISPL residual", extractor(), 0, nlParams.get<double>("Tol dis res L2"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "DISPL residual", extractor(), 0, nlParams.get<double>("Tol dis res Inf"),
      ::NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("DISPL update", extractor(), 0,
          nlParams.get<double>("Tol dis inc L2"), ::NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("DISPL update", extractor(), 0,
          nlParams.get<double>("Tol dis inc Inf"), ::NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  add_status_test(structureDisp_L2);

  // add norm-tests to structural displacement ::NOX::StatusTest::Combo
  structcombo->addStatusTest(structureDisp_L2);
  structcombo->addStatusTest(structureDisp_inf);
  structcombo->addStatusTest(structureDispUpdate_L2);
  structcombo->addStatusTest(structureDispUpdate_inf);

  // add structural displacement test combo to top-level test combo
  converged->addStatusTest(structcombo);
  // ---------- end of structural displacement field tests

  // --------------------------------------------------------------------
  // setup tests for interface
  // --------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> interface;
  interface.push_back(structure_field()->Interface()->FSICondMap());
  interface.push_back(Teuchos::null);
  Core::LinAlg::MultiMapExtractor interfaceextract(*dof_row_map(), interface);

  // create ::NOX::StatusTest::Combo for interface
  Teuchos::RCP<::NOX::StatusTest::Combo> interfacecombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> interfaceTest_L2 = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "GAMMA residual", interfaceextract, 0, nlParams.get<double>("Tol fsi res L2"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> interfaceTest_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "GAMMA residual", interfaceextract, 0, nlParams.get<double>("Tol fsi res Inf"),
      ::NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> interfaceTestUpdate_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("GAMMA update", interfaceextract, 0,
          nlParams.get<double>("Tol fsi inc L2"), ::NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> interfaceTestUpdate_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("GAMMA update", interfaceextract, 0,
          nlParams.get<double>("Tol fsi inc Inf"), ::NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  add_status_test(interfaceTest_L2);

  // add norm-tests to interface ::NOX::StatusTest::Combo
  interfacecombo->addStatusTest(interfaceTest_L2);
  interfacecombo->addStatusTest(interfaceTest_inf);
  interfacecombo->addStatusTest(interfaceTestUpdate_L2);
  interfacecombo->addStatusTest(interfaceTestUpdate_inf);

  // add interface test combo to top-level test combo
  converged->addStatusTest(interfacecombo);
  // ---------- end of interface tests

  // --------------------------------------------------------------------
  // setup tests for fluid velocity field
  // --------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidvel;
  fluidvel.push_back(fluid_field()->InnerVelocityRowMap());
  fluidvel.push_back(Teuchos::null);
  Core::LinAlg::MultiMapExtractor fluidvelextract(*dof_row_map(), fluidvel);

  // create ::NOX::StatusTest::Combo for fluid velocity field
  Teuchos::RCP<::NOX::StatusTest::Combo> fluidvelcombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel_L2 = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "VELOC residual", fluidvelextract, 0, nlParams.get<double>("Tol vel res L2"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "VELOC residual", fluidvelextract, 0, nlParams.get<double>("Tol vel res Inf"),
      ::NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("VELOC update", fluidvelextract, 0,
          nlParams.get<double>("Tol vel inc L2"), ::NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("VELOC update", fluidvelextract, 0,
          nlParams.get<double>("Tol vel inc Inf"), ::NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  add_status_test(innerFluidVel_L2);

  // add norm-tests to fluid velocity ::NOX::StatusTest::Combo
  fluidvelcombo->addStatusTest(innerFluidVel_L2);
  fluidvelcombo->addStatusTest(innerFluidVel_inf);
  fluidvelcombo->addStatusTest(innerFluidVelUpdate_L2);
  fluidvelcombo->addStatusTest(innerFluidVelUpdate_inf);

  // add fluid velocity test combo to top-level test combo
  converged->addStatusTest(fluidvelcombo);
  // ---------- end of fluid velocity field tests

  // --------------------------------------------------------------------
  // setup tests for fluid pressure field
  // --------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidpress;
  fluidpress.push_back(fluid_field()->PressureRowMap());
  fluidpress.push_back(Teuchos::null);
  Core::LinAlg::MultiMapExtractor fluidpressextract(*dof_row_map(), fluidpress);

  // create ::NOX::StatusTest::Combo for fluid pressure field
  Teuchos::RCP<::NOX::StatusTest::Combo> fluidpresscombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress_L2 = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "PRESS residual", fluidpressextract, 0, nlParams.get<double>("Tol pre res L2"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "PRESS residual", fluidpressextract, 0, nlParams.get<double>("Tol pre res Inf"),
      ::NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("PRESS update", fluidpressextract, 0,
          nlParams.get<double>("Tol pre inc L2"), ::NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("PRESS update", fluidpressextract, 0,
          nlParams.get<double>("Tol pre inc Inf"), ::NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  add_status_test(fluidPress_L2);

  // add norm-tests to fluid pressure ::NOX::StatusTest::Combo
  fluidpresscombo->addStatusTest(fluidPress_L2);
  fluidpresscombo->addStatusTest(fluidPress_inf);
  fluidpresscombo->addStatusTest(fluidPressUpdate_L2);
  fluidpresscombo->addStatusTest(fluidPressUpdate_inf);

  // add fluid pressure test combo to top-level test combo
  converged->addStatusTest(fluidpresscombo);
  // ---------- end of fluid pressure field tests

  // Finally, return the test combo
  return combo;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::extract_field_vectors");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (ddgpred_ == Teuchos::null)
  {
    FOUR_C_THROW("Vector 'ddgpred_' has not been initialized properly.");
  }
#endif

  // ---------------------------------------------------------------------------
  // process structure unknowns
  // ---------------------------------------------------------------------------
  // extract structure solution increment from NOX increment
  sx = extractor().ExtractVector(x, 0);

  // ---------------------------------------------------------------------------
  // process ale unknowns
  // ---------------------------------------------------------------------------
  // extract inner ALE solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> aox = extractor().ExtractVector(x, 2);

  // convert structure solution increment to ALE solution increment at the interface
  Teuchos::RCP<Epetra_Vector> scx = structure_field()->Interface()->ExtractFSICondVector(sx);
  scx->Update(1.0, *ddgpred_, 1.0);
  Teuchos::RCP<const Epetra_Vector> acx = struct_to_ale(scx);

  // put inner and interface ALE solution increments together
  Teuchos::RCP<Epetra_Vector> a = ale_field()->Interface()->InsertOtherVector(aox);
  ale_field()->Interface()->InsertFSICondVector(acx, a);
  ax = a;

  // ---------------------------------------------------------------------------
  // process fluid unknowns
  // ---------------------------------------------------------------------------
  // extract inner fluid solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> fox = extractor().ExtractVector(x, 1);
  fox = fluid_field()->Interface()->ExtractOtherVector(fox);

  // convert ALE solution increment to fluid solution increment at the interface
  Teuchos::RCP<Epetra_Vector> fcx = ale_to_fluid_interface(acx);
  fluid_field()->displacement_to_velocity(fcx);

  // put inner and interface fluid solution increments together
  Teuchos::RCP<Epetra_Vector> f = fluid_field()->Interface()->InsertOtherVector(fox);
  fluid_field()->Interface()->InsertFSICondVector(fcx, f);
  fx = f;

  // ---------------------------------------------------------------------------

  // Store field vectors to know them later on as previous quantities
  // inner ale displacement increment
  // interface structure displacement increment
  if (disgprev_ != Teuchos::null)
    ddginc_->Update(1.0, *scx, -1.0, *disgprev_, 0.0);  // compute current iteration increment
  else
    ddginc_ = Teuchos::rcp(new Epetra_Vector(*scx));  // first iteration increment

  disgprev_ = scx;  // store current step increment
  // ------------------------------------

  // inner ale displacement increment
  if (solialeprev_ != Teuchos::null)
    ddialeinc_->Update(1.0, *aox, -1.0, *solialeprev_, 0.0);  // compute current iteration increment
  else
    ddialeinc_ = Teuchos::rcp(new Epetra_Vector(*aox));  // first iteration increment

  solialeprev_ = aox;  // store current step increment
  // ------------------------------------

  // fluid solution increment
  if (soliprev_ != Teuchos::null)  // compute current iteration increment
    duiinc_->Update(1.0, *fox, -1.0, *soliprev_, 0.0);
  else  // first iteration increment
    duiinc_ = Teuchos::rcp(new Epetra_Vector(*fox));
  // store current step increment
  soliprev_ = fox;
  // ------------------------------------
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::output()
{
  structure_field()->Output();
  fluid_field()->Output();

  // output Lagrange multiplier
  OutputLambda();

  ale_field()->Output();

  if (structure_field()->get_constraint_manager()->HaveMonitor())
  {
    structure_field()->get_constraint_manager()->compute_monitor_values(
        structure_field()->Dispnp());
    if (Comm().MyPID() == 0) structure_field()->get_constraint_manager()->PrintMonitorValues();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::OutputLambda()
{
  /* 'lambda_' is only defined on the interface. So, insert 'lambda_' into
   * 'lambdafull' that is defined on the entire fluid field. Then, write
   * output or restart data.
   */
  Teuchos::RCP<Epetra_Vector> lambdafull = fluid_field()->Interface()->InsertFSICondVector(lambda_);
  const int uprestart = timeparams_.get<int>("RESTARTEVRY");
  const int upres = timeparams_.get<int>("RESULTSEVRY");
  if ((uprestart != 0 && fluid_field()->Step() % uprestart == 0) or
      (upres != 0 and fluid_field()->Step() % upres == 0))
    fluid_field()->DiscWriter()->WriteVector("fsilambda", lambdafull);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::update()
{
  lambdaold_->Update(1.0, *lambda_, 0.0);

  FSI::BlockMonolithic::update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::read_restart(int step)
{
  structure_field()->read_restart(step);
  fluid_field()->read_restart(step);

  // read Lagrange multiplier
  {
    Teuchos::RCP<Epetra_Vector> lambdafull =
        Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(), true));
    Core::IO::DiscretizationReader reader = Core::IO::DiscretizationReader(
        fluid_field()->discretization(), Global::Problem::Instance()->InputControlFile(), step);
    reader.ReadVector(lambdafull, "fsilambda");
    lambdaold_ = fluid_field()->Interface()->ExtractFSICondVector(lambdafull);
    // Note: the above is normally enough. However, we can use the restart in order to periodically
    // repeat the fsi simulation (see AC-FS3I)
    lambda_ = fluid_field()->Interface()->ExtractFSICondVector(lambdafull);
  }

  ale_field()->read_restart(step);

  SetTimeStep(fluid_field()->Time(), fluid_field()->Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::prepare_time_step()
{
  increment_time_and_step();
  if (verbosity_ >= Inpar::FSI::verbosity_low) print_header();

  prepare_time_step_preconditioner();

  if (structure_field()->get_stc_algo() != Inpar::STR::stc_none)
    structure_field()->system_matrix()->Reset();

  prepare_time_step_fields();

  // Note: it's important to first prepare the single fields and than the fsi problem
  prepare_time_step_fsi();

  return;
}

/*----------------------------------------------------------------------*/
/* Recover the Lagrange multiplier at the interface   mayr.mt (03/2012) */
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::recover_lagrange_multiplier()
{
  // get time integration parameter of fluid time integrator
  // to enable consistent time integration among the fields
  const double ftiparam = fluid_field()->TimIntParam();

  // some scaling factors for fluid
  const double timescale = fluid_field()->TimeScaling();
  const double scale = fluid_field()->residual_scaling();

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> tmpvec =
      Teuchos::null;  // stores intermediate result of terms (3)-(8)
  Teuchos::RCP<Epetra_Vector> auxvec = Teuchos::null;     // just for convenience
  Teuchos::RCP<Epetra_Vector> auxauxvec = Teuchos::null;  // just for convenience

  /* Recovery of Lagrange multiplier \lambda_^{n+1} is done by the following
   * condensation expression:
   *
   * lambda_^{n+1} =
   *
   * (1)  - ftiparam / (1.-ftiparam) * lambda^{n}
   *
   * (2)  - 1. / (1.-ftiparam) * tmpvec
   *
   * with tmpvec =
   *
   * (3)    r_{\Gamma}^{F,n+1}
   *
   * (4)  + 1 / tau * F_{\Gamma\Gamma} * \Delta d_{\Gamma}^{S,n+1}
   *
   * (5)  + F_{\Gamma\Gamma}^{G} * \Delta d_{\Gamma}^{S,n+1}
   *
   * (6)  + F_{\Gamma I} * \Delta u_{I}^{F,n+1}
   *
   * (7)  + F_{\Gamma I}^{G} * \Delta d_{I}^{G,n+1}
   *
   * (8)  + dt / tau * F_{\Gamma\Gamma} * u_{\Gamma}^n]
   *
   * Remark on term (8):
   * Term (8) has to be considered only in the first Newton iteration.
   * Hence, it will usually not be computed since in general we need more
   * than one nonlinear iteration until convergence.
   *
   * Remarks on all terms:
   * +  Division by -(1.0 - ftiparam) will be done in the end
   *    since this is common to all terms
   * +  tau: time scaling factor for interface time integration (tau =
   * 1/fluid_field()->TimeScaling())
   * +  neglecting terms (4)-(8) should not alter the results significantly
   *    since at the end of the time step the solution increments tend to zero.
   *
   *                                                 Matthias Mayr (10/2012)
   */

  // ---------Addressing term (1)
  lambda_->Update(ftiparam, *lambdaold_, 0.0);
  // ---------End of term (1)

  // ---------Addressing term (3)
  Teuchos::RCP<Epetra_Vector> fluidresidual =
      fluid_field()->Interface()->ExtractFSICondVector(fluid_field()->RHS());
  fluidresidual->Scale(-1.0);  // invert sign to obtain residual, not rhs
  tmpvec = Teuchos::rcp(new Epetra_Vector(*fluidresidual));
  // ---------End of term (3)

  // ---------Addressing term (4)
  auxvec = Teuchos::rcp(new Epetra_Vector(fggprev_->RangeMap(), true));
  fggprev_->Apply(*struct_to_fluid(ddginc_), *auxvec);
  tmpvec->Update(timescale, *auxvec, 1.0);
  // ---------End of term (4)

  // ---------Addressing term (5)
  if (fmggprev_ != Teuchos::null)
  {
    auxvec = Teuchos::rcp(new Epetra_Vector(fmggprev_->RangeMap(), true));
    fmggprev_->Apply(*struct_to_fluid(ddginc_), *auxvec);
    tmpvec->Update(1.0, *auxvec, 1.0);
  }
  // ---------End of term (5)

  // ---------Addressing term (6)
  auxvec = Teuchos::rcp(new Epetra_Vector(fgiprev_->RangeMap(), true));
  fgiprev_->Apply(*duiinc_, *auxvec);
  tmpvec->Update(1.0, *auxvec, 1.0);
  // ---------End of term (6)

  // ---------Addressing term (7)
  if (fmgiprev_ != Teuchos::null)
  {
    /* For matrix-vector-product, the DomainMap() of the matrix and the Map() of the vector
     * have to match. DomainMap() contains inner velocity DOFs and all pressure DOFs.
     * The inner ale displacement increment is converted to the fluid map using AleToFluid().
     * This results in a map that contains all velocity but no pressure DOFs.
     *
     * We have to circumvent some trouble with Epetra_BlockMaps since we cannot split
     * an Epetra_BlockMap into inner and interface DOFs.
     *
     * We create a map extractor 'velothermap' in order to extract the inner velocity
     * DOFs after calling AleToFluid(). Afterwards, a second map extractor 'velotherpressuremapext'
     * is used to append pressure DOFs filled with zeros.
     *
     * Finally, maps match and matrix-vector-multiplication can be done.
     */

    // extract inner velocity DOFs after calling AleToFluid()
    Teuchos::RCP<Epetra_Map> velothermap = Core::LinAlg::SplitMap(
        *fluid_field()->VelocityRowMap(), *interface_fluid_ale_coupling().MasterDofMap());
    Core::LinAlg::MapExtractor velothermapext =
        Core::LinAlg::MapExtractor(*fluid_field()->VelocityRowMap(), velothermap, false);
    auxvec = Teuchos::rcp(new Epetra_Vector(*velothermap, true));
    velothermapext.ExtractOtherVector(
        ale_to_fluid(ale_field()->Interface()->InsertOtherVector(ddialeinc_)), auxvec);

    // add pressure DOFs
    Core::LinAlg::MapExtractor velotherpressuremapext =
        Core::LinAlg::MapExtractor(fmgiprev_->DomainMap(), velothermap);
    auxauxvec = Teuchos::rcp(new Epetra_Vector(fmgiprev_->DomainMap(), true));
    velotherpressuremapext.InsertCondVector(auxvec, auxauxvec);

    // prepare vector to store result of matrix-vector-product
    auxvec = Teuchos::rcp(new Epetra_Vector(fmgiprev_->RangeMap(), true));

    // Now, do the actual matrix-vector-product
    fmgiprev_->Apply(*auxauxvec, *auxvec);
    tmpvec->Update(1.0, *auxvec, 1.0);
  }
  // ---------End of term (7)

  // ---------Addressing term (8)
  if (firstcall_)
  {
    auxvec = Teuchos::rcp(new Epetra_Vector(fggprev_->RangeMap(), true));
    fggprev_->Apply(*fluid_field()->extract_interface_veln(), *auxvec);
    tmpvec->Update(Dt() * timescale, *auxvec, 1.0);
  }
  // ---------End of term (8)

  // ---------Addressing term (2)
  lambda_->Update(scale, *tmpvec, 1.0);  // scale with residual_scaling() to get [N/m^2]
  // ---------End of term (2)

  // Finally, divide by (1.0-ftiparam) which is common to all terms
  lambda_->Scale(-1.0 / (1.0 - ftiparam));

  // Finally, the Lagrange multiplier 'lambda_' is recovered here.
  // It represents nodal forces acting onto the structure.

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::calculate_interface_energy_increment()
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = structure_field()->TimIntParam();
  const double ftiparam = fluid_field()->TimIntParam();

  // interface traction weighted by time integration factors
  Teuchos::RCP<Epetra_Vector> tractionfluid = Teuchos::rcp(new Epetra_Vector(lambda_->Map(), true));
  tractionfluid->Update(stiparam - ftiparam, *lambdaold_, ftiparam - stiparam, *lambda_, 0.0);
  Teuchos::RCP<Epetra_Vector> tractionstructure = fluid_to_struct(tractionfluid);

  // displacement increment of this time step
  Teuchos::RCP<Epetra_Vector> deltad =
      Teuchos::rcp(new Epetra_Vector(*structure_field()->dof_row_map(), true));
  deltad->Update(1.0, *structure_field()->Dispnp(), -1.0, *structure_field()->Dispn(), 0.0);

  // calculate the energy increment
  double energy = 0.0;
  tractionstructure->Dot(*structure_field()->Interface()->ExtractFSICondVector(deltad), &energy);

  energysum_ += energy;

  write_interface_energy_file(energy, energysum_);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFluidSplit::combine_field_vectors(Epetra_Vector& v,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av, const bool slave_vectors_contain_interface_dofs)
{
  if (slave_vectors_contain_interface_dofs)
  {
    // extract inner DOFs from slave vectors
    Teuchos::RCP<Epetra_Vector> fov = fluid_field()->Interface()->ExtractOtherVector(fv);
    fov = fluid_field()->Interface()->InsertOtherVector(fov);
    Teuchos::RCP<Epetra_Vector> aov = ale_field()->Interface()->ExtractOtherVector(av);

    // put them together
    FSI::Monolithic::combine_field_vectors(v, sv, fov, aov);
  }
  else
    FSI::Monolithic::combine_field_vectors(v, sv, fv, av);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::MonolithicFluidSplit::select_dt_error_based() const
{
  // get time step size suggestions based on some error norms
  const double dtstr = get_ada_str_dt();           // based on all structure DOFs
  const double dtstrfsi = get_ada_str_fsi_dt();    // based on structure FSI DOFs
  const double dtflinner = get_ada_fl_inner_dt();  // based on inner fluid DOFs

  double dt = Dt();

  // select time step size based on error estimation
  if (is_ada_structure() and is_ada_fluid())
    dt = std::min(std::min(dtstr, dtstrfsi), dtflinner);
  else if (is_ada_structure() and (not is_ada_fluid()))
    dt = std::min(dtstr, dtstrfsi);
  else if ((not is_ada_structure()) and is_ada_fluid())
    dt = dtflinner;
  else
  {
    // no change in time step size based on structure or fluid field error estimation
  }

  return dt;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicFluidSplit::set_accepted() const
{
  // get error norms
  const double strnorm = get_ada_strnorm();            // based on all structure DOFs
  const double strfsinorm = get_ada_str_fs_inorm();    // based on structure FSI DOFs
  const double flinnernorm = get_ada_fl_inner_norm();  // based on inner fluid DOFs

  bool accepted = std::max(strnorm, strfsinorm) < errtolstr_ and flinnernorm < errtolfl_;

  // in case error estimation in the fluid field is turned off:
  if (not is_ada_fluid()) accepted = std::max(strnorm, strfsinorm) < errtolstr_;

  // in case error estimation in the structure field is turned off:
  if (not is_ada_structure()) accepted = flinnernorm < errtolfl_;

  // no error based time adaptivity
  if ((not is_ada_structure()) and (not is_ada_fluid())) accepted = true;

  return accepted;
}

FOUR_C_NAMESPACE_CLOSE
