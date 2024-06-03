/*--------------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problem with sliding grids using a monolithic scheme
with condensed structure interface displacements

\level 2

*/
/*--------------------------------------------------------------------------*/



#include "4C_fsi_slidingmonolithic_structuresplit.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_constraint_manager.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_fsi_monolithic_linearsystem.hpp"
#include "4C_fsi_overlapprec.hpp"
#include "4C_fsi_overlapprec_fsiamg.hpp"
#include "4C_fsi_statustest.hpp"
#include "4C_fsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_structure_aux.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::SlidingMonolithicStructureSplit::SlidingMonolithicStructureSplit(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : BlockMonolithic(comm, timeparams),
      comm_(comm),
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
  intersectionmaps.push_back(structure_field()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(structure_field()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      CORE::LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  // Check whether the intersection is empty
  if (intersectionmap->NumGlobalElements() != 0)
  {
    //    std::cout << "Slave interface nodes with Dirichlet boundary condition "
    //                 "(input file numbering):" << std::endl;
    //    for (int i=0; i < (int)fluid_field()->discretization()->NumMyRowNodes(); i++)
    //    {
    //      // get all nodes and add them
    //      int gid = structure_field()->discretization()->NodeRowMap()->GID(i);
    //
    //      // do only nodes that I have in my discretization
    //      if (!structure_field()->discretization()->NodeColMap()->MyGID(gid)) continue;
    //      CORE::Nodes::Node* node = structure_field()->discretization()->gNode(gid);
    //      if (!node) FOUR_C_THROW("Cannot find node with gid %",gid);
    //
    //      std::vector<int> nodedofs = structure_field()->discretization()->Dof(node);
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

    FOUR_C_THROW(errormsg.str());
  }
  // ---------------------------------------------------------------------------

  notsetup_ = true;

  coupsfm_ = Teuchos::rcp(new CORE::ADAPTER::CouplingMortar(GLOBAL::Problem::Instance()->NDim(),
      GLOBAL::Problem::Instance()->mortar_coupling_params(),
      GLOBAL::Problem::Instance()->contact_dynamic_params(),
      GLOBAL::Problem::Instance()->spatial_approximation_type()));
  fscoupfa_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());

  aigtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fsaigtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fsmgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);

  SetLambda();
  ddiinc_ = Teuchos::null;
  disiprev_ = Teuchos::null;
  disgprev_ = Teuchos::null;
  sgiprev_ = Teuchos::null;
  sggprev_ = Teuchos::null;

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (coupsfm_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'coupsfm_' failed.");
  }
  if (fscoupfa_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'fscoupfa_' failed.");
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
  if (fsaigtransform_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'fsaigtransform_' failed.");
  }
  if (fsmgitransform_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'fsmgitransform_' failed.");
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

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::SetLambda()
{
  lambda_ = Teuchos::rcp(new Epetra_Vector(*structure_field()->Interface()->FSICondMap(), true));
  lambdaold_ = Teuchos::rcp(new Epetra_Vector(*structure_field()->Interface()->FSICondMap(), true));

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::SetupSystem()
{
  if (notsetup_)
  {
    const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
    linearsolverstrategy_ =
        CORE::UTILS::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

    aleproj_ = CORE::UTILS::IntegralValue<INPAR::FSI::SlideALEProj>(fsidyn, "SLIDEALEPROJ");

    set_default_parameters(fsidyn, nox_parameter_list());

    // we use non-matching meshes at the interface
    // mortar with: structure = slave, fluid = master

    const int ndim = GLOBAL::Problem::Instance()->NDim();

    // get coupling objects
    CORE::ADAPTER::Coupling& icoupfa = interface_fluid_ale_coupling();

    /* structure to fluid
     * coupling condition at the fsi interface:
     * displacements (=number spatial dimensions) are coupled
     * e.g.: 3D: coupleddof = [1, 1, 1]
     */
    std::vector<int> coupleddof(ndim, 1);

    coupsfm_->Setup(fluid_field()->discretization(), structure_field()->discretization(),
        ale_field()->write_access_discretization(), coupleddof, "FSICoupling", comm_,
        GLOBAL::Problem::Instance()->FunctionManager(), false);

    // fluid to ale at the interface
    icoupfa.setup_condition_coupling(*fluid_field()->discretization(),
        fluid_field()->Interface()->FSICondMap(), *ale_field()->discretization(),
        ale_field()->Interface()->FSICondMap(), "FSICoupling", ndim);

    // we might have a free surface
    if (fluid_field()->Interface()->FSCondRelevant())
    {
      fscoupfa_->setup_condition_coupling(*fluid_field()->discretization(),
          fluid_field()->Interface()->FSCondMap(), *ale_field()->discretization(),
          ale_field()->Interface()->FSCondMap(), "FREESURFCoupling", ndim);
    }

    // the fluid-ale coupling always matches
    const Epetra_Map* fluidnodemap = fluid_field()->discretization()->NodeRowMap();
    const Epetra_Map* alenodemap = ale_field()->discretization()->NodeRowMap();

    CORE::ADAPTER::Coupling& coupfa = fluid_ale_coupling();

    coupfa.setup_coupling(*fluid_field()->discretization(), *ale_field()->discretization(),
        *fluidnodemap, *alenodemap, ndim);

    fluid_field()->SetMeshMap(coupfa.MasterDofMap());

    // create combined map
    create_combined_dof_row_map();

    // Use normal matrix for fluid equations but build (splitted) mesh movement
    // linearization (if requested in the input file)
    fluid_field()->use_block_matrix(false);

    // Use splitted structure matrix
    structure_field()->use_block_matrix();

    // build ale system matrix in splitted system
    ale_field()->create_system_matrix(ale_field()->Interface());

    aleresidual_ = Teuchos::rcp(new Epetra_Vector(*FsiAleField()->FsiInterface()->OtherMap()));

    // -------------------------------------------------------------------------
    // Build the global Dirichlet map extractor
    setup_dbc_map_extractor();
    // -------------------------------------------------------------------------

    // enable debugging
    if (CORE::UTILS::IntegralValue<int>(fsidyn, "DEBUGOUTPUT") & 2)
    {
      pcdbg_ = Teuchos::rcp(new UTILS::MonolithicDebugWriter(*this));
    }

    create_system_matrix();

    // set up sliding ale if necessary
    if (aleproj_ != INPAR::FSI::ALEprojection_none)
    {
      // mesh_init possibly modifies reference configuration of slave side --> recompute element
      // volume in initialize_elements()
      structure_field()->discretization()->fill_complete(false, true, true);
      // set up sliding ale utils
      slideale_ = Teuchos::rcp(new FSI::UTILS::SlideAleUtils(structure_field()->discretization(),
          fluid_field()->discretization(), *coupsfm_, false, aleproj_));

      iprojdisp_ = Teuchos::rcp(new Epetra_Vector(*coupsfm_->MasterDofMap(), true));
      iprojdispinc_ = Teuchos::rcp(new Epetra_Vector(*coupsfm_->MasterDofMap(), true));
    }
    notsetup_ = false;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::create_combined_dof_row_map()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(structure_field()->Interface()->OtherMap());
  vecSpaces.push_back(fluid_field()->dof_row_map());
  vecSpaces.push_back(FsiAleField()->FsiInterface()->OtherMap());

  if (vecSpaces[0]->NumGlobalElements() == 0)
    FOUR_C_THROW("No inner structural equations. Splitting not possible. Panic.");

  set_dof_row_maps(vecSpaces);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::setup_dbc_map_extractor()
{
  /* Dirichlet maps for structure and fluid do not intersect with interface map.
   * ALE Dirichlet map might intersect with interface map, but ALE interface
   * DOFs are not part of the final system of equations. Hence, we just need the
   * intersection of inner ALE DOFs with Dirichlet ALE DOFs.
   */
  std::vector<Teuchos::RCP<const Epetra_Map>> aleintersectionmaps;
  aleintersectionmaps.push_back(ale_field()->GetDBCMapExtractor()->CondMap());
  aleintersectionmaps.push_back(FsiAleField()->FsiInterface()->OtherMap());
  Teuchos::RCP<Epetra_Map> aleintersectionmap =
      CORE::LINALG::MultiMapExtractor::IntersectMaps(aleintersectionmaps);

  // Merge Dirichlet maps of structure, fluid and ALE to global FSI Dirichlet map
  std::vector<Teuchos::RCP<const Epetra_Map>> dbcmaps;
  dbcmaps.push_back(structure_field()->GetDBCMapExtractor()->CondMap());
  dbcmaps.push_back(fluid_field()->GetDBCMapExtractor()->CondMap());
  dbcmaps.push_back(aleintersectionmap);
  Teuchos::RCP<const Epetra_Map> dbcmap = CORE::LINALG::MultiMapExtractor::MergeMaps(dbcmaps);

  // Finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new CORE::LINALG::MapExtractor(*dof_row_map(), dbcmap, true));
  if (dbcmaps_ == Teuchos::null) FOUR_C_THROW("Creation of FSI Dirichlet map extractor failed.");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>
FSI::SlidingMonolithicStructureSplit::SystemMatrix() const
{
  return systemmatrix_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::setup_rhs_residual(Epetra_Vector& f)
{
  /* get time integration parameters of structure and fluid time integrators
   * to enable consistent time integration among the fields
   */
  const double stiparam = structure_field()->TimIntParam();
  const double ftiparam = fluid_field()->TimIntParam();

  // some scaling factors for fluid
  const double fluidscale = fluid_field()->residual_scaling();

  // get the Mortar matrix M
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarMatrixP();

  // get single field residuals
  Teuchos::RCP<const Epetra_Vector> sv = Teuchos::rcp(new Epetra_Vector(*structure_field()->RHS()));
  Teuchos::RCP<const Epetra_Vector> fv = Teuchos::rcp(new Epetra_Vector(*fluid_field()->RHS()));
  Teuchos::RCP<const Epetra_Vector> av = Teuchos::rcp(new Epetra_Vector(*ale_field()->RHS()));

  // extract only inner DOFs from structure (=slave) and ALE field
  Teuchos::RCP<const Epetra_Vector> sov = structure_field()->Interface()->ExtractOtherVector(sv);
  Teuchos::RCP<const Epetra_Vector> aov = FsiAleField()->FsiInterface()->ExtractOtherVector(av);

  // add structure interface residual to fluid interface residual considering temporal scaling
  Teuchos::RCP<const Epetra_Vector> scv = structure_field()->Interface()->ExtractFSICondVector(sv);
  Teuchos::RCP<Epetra_Vector> fcv =
      CORE::LINALG::CreateVector(*fluid_field()->Interface()->FSICondMap(), true);
  mortarp->Multiply(true, *scv, *fcv);
  Teuchos::RCP<Epetra_Vector> modfv = fluid_field()->Interface()->InsertFSICondVector(fcv);
  modfv->Update(1.0, *fv, (1.0 - ftiparam) / ((1.0 - stiparam) * fluidscale));

  // put the single field residuals together
  FSI::Monolithic::combine_field_vectors(f, sov, modfv, aov);

  // add additional ale residual
  extractor().AddVector(*aleresidual_, 2, f);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::setup_rhs_lambda(Epetra_Vector& f)
{
  if (lambdaold_ != Teuchos::null)
  {
    // get time integration parameters of structure and fluid time integrators
    // to enable consistent time integration among the fields
    const double stiparam = structure_field()->TimIntParam();
    const double ftiparam = fluid_field()->TimIntParam();

    // some scaling factors for fluid
    const double fluidscale = fluid_field()->residual_scaling();

    // get the Mortar matrix M
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> mortarm = coupsfm_->GetMortarMatrixM();

    // project Lagrange multiplier field onto the master interface DOFs and consider temporal
    // scaling
    Teuchos::RCP<Epetra_Vector> lambda =
        Teuchos::rcp(new Epetra_Vector(mortarm->DomainMap(), true));
    mortarm->Multiply(true, *lambdaold_, *lambda);
    Teuchos::RCP<Epetra_Vector> lambdafull =
        fluid_field()->Interface()->InsertFSICondVector(lambda);
    lambdafull->Scale((-ftiparam + (stiparam * (1.0 - ftiparam)) / (1.0 - stiparam)) / fluidscale);

    // add Lagrange multiplier
    extractor().AddVector(*lambdafull, 1, f);
  }

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::setup_rhs_firstiter(Epetra_Vector& f)
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = structure_field()->TimIntParam();
  const double ftiparam = fluid_field()->TimIntParam();

  // some scaling factors for fluid
  const double scale = fluid_field()->residual_scaling();

  // old interface velocity of fluid field
  const Teuchos::RCP<const Epetra_Vector> fveln = fluid_field()->extract_interface_veln();

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<const CORE::LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarMatrixP();

  // get fluid shape derivatives matrix
  const Teuchos::RCP<const CORE::LINALG::BlockSparseMatrixBase> mmm =
      fluid_field()->ShapeDerivatives();

  // get structure matrix
  const Teuchos::RCP<const CORE::LINALG::BlockSparseMatrixBase> blocks =
      structure_field()->BlockSystemMatrix();

  // get ale matrix
  const Teuchos::RCP<const CORE::LINALG::BlockSparseMatrixBase> blocka =
      ale_field()->BlockSystemMatrix();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (mortarp == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix P.");
  if (blocks == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to structure block matrix.");
  if (blocka == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to ALE block matrix.");
#endif

  // extract submatrices
  const CORE::LINALG::SparseMatrix& sig = blocks->Matrix(0, 1);  // S_{I\Gamma}
  const CORE::LINALG::SparseMatrix& sgg = blocks->Matrix(1, 1);  // S_{\Gamma\Gamma}
  const CORE::LINALG::SparseMatrix& aig = blocka->Matrix(0, 1);  // A_{I\Gamma}

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;     // right hand side of single set of DOFs
  Teuchos::RCP<Epetra_Vector> auxvec = Teuchos::null;  // just for convenience
  Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::null;  // just for convenience

  /* Different contributions/terms to the rhs are separated by the following
   * comment line */
  // ---------- inner structure DOFs
  /* The following terms are added to the inner structure DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - dt * S_{I \Gamma} * P * u^{n}_{\Gamma}
   *
   * (2)  + S_{I \Gamma} * \Delta d_{\Gamma,p}
   *
   * Remarks on all terms:
   * +  tau: time scaling factor for interface time integration
   *         (tau = 1/fluid_field()->TimeScaling())
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(sig.RangeMap(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(), true));

  mortarp->Apply(*fveln, *auxvec);
  sig.Apply(*auxvec, *rhs);

  rhs->Scale(-Dt());

  extractor().AddVector(*rhs, 0, f);
  // ----------end of term 1

  // ----------addressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(sig.RangeMap(), true));

  sig.Apply(*ddgpred_, *rhs);

  extractor().AddVector(*rhs, 0, f);
  // ----------end of term 2
  // ----------end of inner structure DOFs

  // ---------- inner fluid DOFs
  /* The following terms are added to the inner fluid DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - dt * F^{G}_{I \Gamma} * u^{n}_{\Gamma}
   *
   */
  // ----------addressing term 1
  if (mmm != Teuchos::null)
  {
    // extract F^{G}_{I \Gamma}
    const CORE::LINALG::SparseMatrix& fmig = mmm->Matrix(0, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmig.RangeMap(), true));

    fmig.Apply(*fveln, *rhs);

    rhs->Scale(-Dt());
    rhs = FsiFluidField()->FsiInterface()->InsertOtherVector(rhs);

    extractor().AddVector(*rhs, 1, f);
  }
  // ----------end of term 1
  // ----------end of inner fluid DOFs

  // ---------- interface fluid DOFs
  /* The following terms are added to the interface fluid DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - dt * F^{G}_{\Gamma\Gamma} * u^{n}_{\Gamma}
   *
   * (2)  - (1-ftiparam) / (1-stiparam) * dt * P^{T} * S_{\Gamma\Gamma} * P * u^{n}_{\Gamma}
   *
   * (3)  + (1-ftiparam) / (1-stiparam) * P^{T} * S_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
   *
   * Remarks on all terms:
   * +  tau: time scaling factor for interface time integration (tau =
   * 1/fluid_field()->TimeScaling())
   *
   */
  // ----------addressing term 1
  if (mmm != Teuchos::null)
  {
    // extract F^{G}_{\Gamma\Gamma}
    const CORE::LINALG::SparseMatrix& fmgg = mmm->Matrix(1, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RangeMap(), true));

    fmgg.Apply(*fveln, *rhs);

    rhs->Scale(-Dt());
    rhs = fluid_field()->Interface()->InsertFSICondVector(rhs);

    extractor().AddVector(*rhs, 1, f);
  }
  // ----------end of term 1

  // ----------addressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(mortarp->DomainMap(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(sgg.RangeMap(), true));
  tmpvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(), true));

  mortarp->Apply(*fveln, *tmpvec);
  sgg.Apply(*tmpvec, *auxvec);
  mortarp->Multiply(true, *auxvec, *rhs);

  rhs->Scale(-(1. - ftiparam) / (1. - stiparam) * Dt() / scale);
  rhs = fluid_field()->Interface()->InsertFSICondVector(rhs);

  extractor().AddVector(*rhs, 1, f);
  // ----------end of term 2

  // ----------addressing term 3
  rhs = Teuchos::rcp(new Epetra_Vector(mortarp->DomainMap(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(sgg.RangeMap(), true));

  sgg.Apply(*ddgpred_, *auxvec);
  mortarp->Multiply(true, *auxvec, *rhs);

  rhs->Scale((1. - ftiparam) / (1. - stiparam) / scale);
  rhs = fluid_field()->Interface()->InsertFSICondVector(rhs);

  extractor().AddVector(*rhs, 1, f);
  // ----------end of term 3
  // ----------end of interface fluid DOFs

  // ---------- inner ALE DOFs
  /* The following terms are added to the inner ALE DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - dt * A_{I \Gamma} * u^{n}_{\Gamma}
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(aig.RangeMap(), true));

  aig.Apply(*fluid_to_ale_interface(fveln), *rhs);

  rhs->Scale(-Dt());

  extractor().AddVector(*rhs, 2, f);
  // ----------end of term 1
  // ----------end of inner ALE DOFs

  // only if relative movement between ale and structure is possible
  if (aleproj_ != INPAR::FSI::ALEprojection_none)
  {
    rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap(), true));

    aig.Apply(*fluid_to_ale_interface(iprojdispinc_), *rhs);

    extractor().AddVector(*rhs, 2, f);
  }

  // if there is a free surface
  if (fluid_field()->Interface()->FSCondRelevant())
  {
    // here we extract the free surface submatrices from position 2
    const CORE::LINALG::SparseMatrix& afsig = blocka->Matrix(0, 2);

    // extract fluid free surface velocities.
    Teuchos::RCP<Epetra_Vector> fveln = fluid_field()->extract_free_surface_veln();
    Teuchos::RCP<Epetra_Vector> aveln = fluid_to_ale_interface(fveln);

    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(afsig.RowMap(), true));
    afsig.Apply(*aveln, *rhs);

    rhs->Scale(-1. * Dt());

    extractor().AddVector(*rhs, 1, f);

    // shape derivatives
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> mmm = fluid_field()->ShapeDerivatives();
    if (mmm != Teuchos::null)
    {
      // here we extract the free surface submatrices from position 2
      CORE::LINALG::SparseMatrix& fmfsig = mmm->Matrix(0, 2);
      CORE::LINALG::SparseMatrix& fmfsgg = mmm->Matrix(2, 2);

      rhs = Teuchos::rcp(new Epetra_Vector(fmfsig.RowMap(), true));
      fmfsig.Apply(*fveln, *rhs);
      Teuchos::RCP<Epetra_Vector> veln = FsiFluidField()->FsiInterface()->InsertOtherVector(rhs);

      rhs = Teuchos::rcp(new Epetra_Vector(fmfsgg.RowMap(), true));
      fmfsgg.Apply(*fveln, *rhs);
      fluid_field()->Interface()->InsertFSCondVector(rhs, veln);

      veln->Scale(-1. * Dt());

      extractor().AddVector(*veln, 0, f);
    }
  }

  /* Reset quantities for previous iteration step since they still store values
   * from the last time step */
  ddiinc_ = CORE::LINALG::CreateVector(*structure_field()->Interface()->OtherMap(), true);
  disiprev_ = Teuchos::null;
  disgprev_ = Teuchos::null;
  sgicur_ = Teuchos::null;
  sggcur_ = Teuchos::null;

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::setup_system_matrix(
    CORE::LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::SlidingMonolithicStructureSplit::setup_system_matrix");

  // get the Mortar projection matrix P = D^{-1} * M
  Teuchos::RCP<CORE::LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarMatrixP();

  // get single field block matrices
  const Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> s =
      structure_field()->BlockSystemMatrix();
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> f = fluid_field()->SystemMatrix();
  const Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> a = ale_field()->BlockSystemMatrix();


#ifdef FOUR_C_ENABLE_ASSERTIONS
  // check whether allocation was successful
  if (mortarp == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix P.");
  if (s == Teuchos::null)
  {
    FOUR_C_THROW("expect structure block matrix");
  }
  if (f == Teuchos::null)
  {
    FOUR_C_THROW("expect fluid matrix");
  }
  if (a == Teuchos::null)
  {
    FOUR_C_THROW("expect ale block matrix");
  }

  // some checks whether maps for matrix-matrix-multiplication do really match
  if (!s->Matrix(0, 1).DomainMap().PointSameAs(mortarp->RangeMap()))
    FOUR_C_THROW("Maps do not match.");
  if (!s->Matrix(1, 0).RangeMap().PointSameAs(mortarp->RangeMap()))
    FOUR_C_THROW("Maps do not match.");
  if (!s->Matrix(1, 1).DomainMap().PointSameAs(mortarp->RangeMap()))
    FOUR_C_THROW("Maps do not match.");
#endif

  // extract submatrices
  CORE::LINALG::SparseMatrix& aii = a->Matrix(0, 0);
  const CORE::LINALG::SparseMatrix& aig = a->Matrix(0, 1);

  // scaling factors for fluid
  const double scale = fluid_field()->residual_scaling();
  const double timescale = fluid_field()->TimeScaling();

  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = structure_field()->TimIntParam();
  const double ftiparam = fluid_field()->TimIntParam();

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  f->UnComplete();

  // ---------------------------------------------------------------------------
  // BEGIN building the global 4x4 system matrix
  // ---------------------------------------------------------------------------
  // Contributions to blocks in system matrix are listed separately.
  // Block numbering in comments ranges from (1,1) to (4,4).

  mat.Assign(0, 0, CORE::LINALG::View, s->Matrix(0, 0));

  // ----------Addressing contribution to block (1,3)
  Teuchos::RCP<CORE::LINALG::SparseMatrix> sig =
      MLMultiply(s->Matrix(0, 1), false, *mortarp, false, false, false, true);
  Teuchos::RCP<CORE::LINALG::SparseMatrix> lsig =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(sig->RowMap(), 81, false));

  lsig->Add(*sig, false, 1. / timescale, 0.0);
  lsig->Complete(f->DomainMap(), sig->RangeMap());

  mat.Assign(0, 1, CORE::LINALG::View, *lsig);

  // ----------Addressing contribution to block (3,1)
  Teuchos::RCP<CORE::LINALG::SparseMatrix> sgi =
      MLMultiply(*mortarp, true, s->Matrix(1, 0), false, false, false, true);
  Teuchos::RCP<CORE::LINALG::SparseMatrix> lsgi =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(f->RowMap(), 81, false));

  lsgi->Add(*sgi, false, (1. - ftiparam) / ((1. - stiparam) * scale), 0.0);
  lsgi->Complete(sgi->DomainMap(), f->RangeMap());

  mat.Assign(1, 0, CORE::LINALG::View, *lsgi);

  // ----------Addressing contribution to block (3,3)
  Teuchos::RCP<CORE::LINALG::SparseMatrix> sgg =
      MLMultiply(s->Matrix(1, 1), false, *mortarp, false, false, false, true);
  sgg = MLMultiply(*mortarp, true, *sgg, false, false, false, true);

  f->Add(*sgg, false, (1. - ftiparam) / ((1. - stiparam) * scale * timescale), 1.0);
  mat.Assign(1, 1, CORE::LINALG::View, *f);

  (*aigtransform_)(a->FullRowMap(), a->FullColMap(), aig, 1. / timescale,
      CORE::ADAPTER::CouplingSlaveConverter(interface_fluid_ale_coupling()), mat.Matrix(2, 1));
  mat.Assign(2, 2, CORE::LINALG::View, aii);

  /*--------------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> mmm = fluid_field()->ShapeDerivatives();
  if (mmm != Teuchos::null)
  {
    // extract submatrices
    const CORE::LINALG::SparseMatrix& fmii = mmm->Matrix(0, 0);
    const CORE::LINALG::SparseMatrix& fmig = mmm->Matrix(0, 1);
    const CORE::LINALG::SparseMatrix& fmgi = mmm->Matrix(1, 0);
    const CORE::LINALG::SparseMatrix& fmgg = mmm->Matrix(1, 1);

    // ----------Addressing contribution to block (3,3)
    mat.Matrix(1, 1).Add(fmgg, false, 1. / timescale, 1.0);

    // ----------Addressing contribution to block (2,3)
    mat.Matrix(1, 1).Add(fmig, false, 1. / timescale, 1.0);

    const CORE::ADAPTER::Coupling& coupfa = fluid_ale_coupling();

    (*fmgitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmgi, 1.,
        CORE::ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, false);

    (*fmiitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmii, 1.,
        CORE::ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, true);
  }

  // if there is a free surface
  if (fluid_field()->Interface()->FSCondRelevant())
  {
    // here we extract the free surface submatrices from position 2
    CORE::LINALG::SparseMatrix& aig = a->Matrix(0, 2);

    (*fsaigtransform_)(a->FullRowMap(), a->FullColMap(), aig, 1. / timescale,
        CORE::ADAPTER::CouplingSlaveConverter(*fscoupfa_), mat.Matrix(2, 1));

    if (mmm != Teuchos::null)
    {
      // We assume there is some space between fsi interface and free
      // surface. Thus the matrices mmm->Matrix(1,2) and mmm->Matrix(2,1) are
      // zero.

      // here we extract the free surface submatrices from position 2
      CORE::LINALG::SparseMatrix& fmig = mmm->Matrix(0, 2);
      CORE::LINALG::SparseMatrix& fmgi = mmm->Matrix(2, 0);
      CORE::LINALG::SparseMatrix& fmgg = mmm->Matrix(2, 2);

      mat.Matrix(1, 1).Add(fmgg, false, 1. / timescale, 1.0);
      mat.Matrix(1, 1).Add(fmig, false, 1. / timescale, 1.0);

      const CORE::ADAPTER::Coupling& coupfa = fluid_ale_coupling();

      (*fsmgitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmgi, 1.,
          CORE::ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(1, 2), false, false);
    }
  }

  // done. make sure all blocks are filled.
  mat.Complete();

  // Finally, take care of Dirichlet boundary conditions
  mat.ApplyDirichlet(*(dbcmaps_->CondMap()), true);
  //
  // ---------------------------------------------------------------------------
  // END building the global system matrix
  // ---------------------------------------------------------------------------

  // store parts of structural matrix to know them in the next iteration as previous iteration
  // matrices
  sgiprev_ = sgicur_;
  sggprev_ = sggcur_;
  sgicur_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(s->Matrix(1, 0)));
  sggcur_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(s->Matrix(1, 1)));
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::update()
{
  lambdaold_->Update(1.0, *lambda_, 0.0);

  // update history variables for sliding ale
  if (aleproj_ != INPAR::FSI::ALEprojection_none)
  {
    iprojdisp_ = Teuchos::rcp(new Epetra_Vector(*coupsfm_->MasterDofMap(), true));
    Teuchos::RCP<Epetra_Vector> idispale = ale_to_fluid_interface(
        ale_field()->Interface()->ExtractFSICondVector(ale_field()->Dispnp()));

    slideale_->Remeshing(*structure_field(), fluid_field()->discretization(), idispale, iprojdisp_,
        *coupsfm_, Comm());

    iprojdispinc_->Update(-1.0, *iprojdisp_, 1.0, *idispale, 0.0);

    slideale_->EvaluateMortar(structure_field()->extract_interface_dispnp(), iprojdisp_, *coupsfm_);
    slideale_->EvaluateFluidMortar(idispale, iprojdisp_);

    Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*iprojdisp_));
    temp->ReplaceMap(idispale->Map());
    Teuchos::RCP<Epetra_Vector> acx = fluid_to_ale_interface(temp);
    ale_field()->apply_interface_displacements(acx);
    fluid_field()->apply_mesh_displacement(ale_to_fluid(ale_field()->Dispnp()));

    Teuchos::RCP<Epetra_Vector> unew =
        slideale_->InterpolateFluid(fluid_field()->extract_interface_velnp());
    fluid_field()->apply_interface_velocities(unew);
  }

  // call Update()-routine in base class to handle the single fields
  FSI::BlockMonolithic::update();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::scale_system(
    CORE::LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)CORE::UTILS::IntegralValue<int>(fsimono, "INFNORMSCALING");

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


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::unscale_solution(
    CORE::LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)CORE::UTILS::IntegralValue<int>(fsimono, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> sy = extractor().ExtractVector(x, 0);
    Teuchos::RCP<Epetra_Vector> ay = extractor().ExtractVector(x, 2);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ay->Multiply(1.0, *acolsum_, *ay, 0.0)) FOUR_C_THROW("ale scaling failed");

    extractor().InsertVector(*sy, 0, x);
    extractor().InsertVector(*ay, 2, x);

    Teuchos::RCP<Epetra_Vector> sx = extractor().ExtractVector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = extractor().ExtractVector(b, 2);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ax->ReciprocalMultiply(1.0, *arowsum_, *ax, 0.0)) FOUR_C_THROW("ale scaling failed");

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
  utils()->out() << std::scientific << "\nlinear solver quality:\n"
                 << "L_2-norms:\n"
                 << "   |r|=" << n << "   |rs|=" << ns << "   |rf|=" << nf << "   |ra|=" << na
                 << "\n";
  r.NormInf(&n);
  sr->NormInf(&ns);
  fr->NormInf(&nf);
  ar->NormInf(&na);
  utils()->out() << "L_inf-norms:\n"
                 << "   |r|=" << n << "   |rs|=" << ns << "   |rf|=" << nf << "   |ra|=" << na
                 << "\n";

  utils()->out().flags(flags);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Epetra::LinearSystem>
FSI::SlidingMonolithicStructureSplit::create_linear_system(Teuchos::ParameterList& nlParams,
    ::NOX::Epetra::Vector& noxSoln, Teuchos::RCP<::NOX::Utils> utils)
{
  Teuchos::RCP<::NOX::Epetra::LinearSystem> linSys;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  //  Teuchos::ParameterList* lsParams = nullptr;
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  //  // in case of nonlinCG the linear solver list is somewhere else
  //  if (dirParams.get("Method","User Defined")=="User Defined")
  //    lsParams = &(newtonParams.sublist("Linear Solver"));
  //  else if (dirParams.get("Method","User Defined")=="NonlinearCG")
  //    lsParams = &(dirParams.sublist("Nonlinear CG").sublist("Linear Solver"));
  //  else FOUR_C_THROW("Unknown nonlinear method");

  ::NOX::Epetra::Interface::Jacobian* iJac = this;
  ::NOX::Epetra::Interface::Preconditioner* iPrec = this;
  const Teuchos::RCP<Epetra_Operator> J = systemmatrix_;
  const Teuchos::RCP<Epetra_Operator> M = systemmatrix_;

  switch (linearsolverstrategy_)
  {
    case INPAR::FSI::PreconditionedKrylov:
      linSys = Teuchos::rcp(new ::NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
          Teuchos::rcp(iJac, false), J, Teuchos::rcp(iPrec, false), M, noxSoln));
      break;
    default:
      FOUR_C_THROW("unsupported linear block solver strategy: %d", linearsolverstrategy_);
      break;
  }

  return linSys;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Combo> FSI::SlidingMonolithicStructureSplit::create_status_test(
    Teuchos::ParameterList& nlParams, Teuchos::RCP<::NOX::Epetra::Group> grp)
{
  // ---------------------------------------------------------------------------
  // Setup the test framework
  // ---------------------------------------------------------------------------
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


  // ---------------------------------------------------------------------------
  // setup tests for structural displacement field
  // ---------------------------------------------------------------------------
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

  // ---------------------------------------------------------------------------
  // setup tests for interface
  // ---------------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> interface;
  interface.push_back(fluid_field()->Interface()->FSICondMap());
  interface.push_back(Teuchos::null);
  CORE::LINALG::MultiMapExtractor interfaceextract(*dof_row_map(), interface);

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

  // ---------------------------------------------------------------------------
  // setup tests for fluid velocity field
  // ---------------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidvel;
  fluidvel.push_back(fluid_field()->InnerVelocityRowMap());
  fluidvel.push_back(Teuchos::null);
  CORE::LINALG::MultiMapExtractor fluidvelextract(*dof_row_map(), fluidvel);

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

  // ---------------------------------------------------------------------------
  // setup tests for fluid pressure field
  // ---------------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidpress;
  fluidpress.push_back(fluid_field()->PressureRowMap());
  fluidpress.push_back(Teuchos::null);
  CORE::LINALG::MultiMapExtractor fluidpressextract(*dof_row_map(), fluidpress);

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


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::extract_field_vectors(
    Teuchos::RCP<const Epetra_Vector> x, Teuchos::RCP<const Epetra_Vector>& sx,
    Teuchos::RCP<const Epetra_Vector>& fx, Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::SlidingMonolithicStructureSplit::extract_field_vectors");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (ddgpred_ == Teuchos::null)
    FOUR_C_THROW("Vector 'ddgpred_' has not been initialized properly.");
#endif

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<const CORE::LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarMatrixP();

  // ---------------------------------------------------------------------------
  // process fluid unknowns
  // ---------------------------------------------------------------------------
  // extract fluid solution increment from NOX increment
  Teuchos::RCP<Epetra_Vector> f = extractor().ExtractVector(x, 1);
  fluid_field()->UpdateSlaveDOF(f);
  fx = f;

  // ---------------------------------------------------------------------------
  // process ale unknowns
  // ---------------------------------------------------------------------------
  // extract inner ALE solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> aox = extractor().ExtractVector(x, 2);

  // convert fluid interface velocities into ALE interface displacements
  Teuchos::RCP<Epetra_Vector> fcx = fluid_field()->Interface()->ExtractFSICondVector(fx);
  fluid_field()->velocity_to_displacement(fcx);
  Teuchos::RCP<Epetra_Vector> acx = fluid_to_ale_interface(fcx);

  // put inner and interface ALE solution increments together
  Teuchos::RCP<Epetra_Vector> a = FsiAleField()->FsiInterface()->InsertOtherVector(aox);
  ale_field()->Interface()->InsertFSICondVector(acx, a);

  // if there is a free surface
  if (fluid_field()->Interface()->FSCondRelevant())
  {
    Teuchos::RCP<Epetra_Vector> fcx = fluid_field()->Interface()->ExtractFSCondVector(fx);
    fluid_field()->free_surf_velocity_to_displacement(fcx);

    Teuchos::RCP<Epetra_Vector> acx = fscoupfa_->MasterToSlave(fcx);
    ale_field()->Interface()->InsertFSCondVector(acx, a);
  }

  ale_field()->UpdateSlaveDOF(a);
  ax = a;

  // ---------------------------------------------------------------------------
  // process structure unknowns
  // ---------------------------------------------------------------------------
  // extract inner structure solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> sox = extractor().ExtractVector(x, 0);

  // convert ALE interface displacements to structure interface displacements
  Teuchos::RCP<Epetra_Vector> scx =
      CORE::LINALG::CreateVector(*structure_field()->Interface()->FSICondMap());
  acx = ale_to_fluid_interface(acx);
  mortarp->Apply(*acx, *scx);
  scx->Update(-1.0, *ddgpred_, 1.0);

  // put inner and interface structure solution increments together
  Teuchos::RCP<Epetra_Vector> s = structure_field()->Interface()->InsertOtherVector(sox);
  structure_field()->Interface()->InsertFSICondVector(scx, s);
  sx = s;

  // ---------------------------------------------------------------------------

  // Store field vectors to know them later on as previous quantities
  if (disiprev_ != Teuchos::null)
    ddiinc_->Update(1.0, *sox, -1.0, *disiprev_, 0.0);  // compute current iteration increment
  else
    ddiinc_ = Teuchos::rcp(new Epetra_Vector(*sox));  // first iteration increment

  disiprev_ = sox;  // store current step increment

  if (velgprev_ != Teuchos::null)
    duginc_->Update(1.0, *fcx, -1.0, *velgprev_, 0.0);  // compute current iteration increment
  else
    duginc_ = Teuchos::rcp(new Epetra_Vector(*fcx));  // first iteration increment

  velgprev_ = fcx;  // store current step increment
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::output()
{
  structure_field()->Output();

  // output Lagrange multiplier
  OutputLambda();

  fluid_field()->Output();

  if (aleproj_ != INPAR::FSI::ALEprojection_none)
  {
    int uprestart = timeparams_.get<int>("RESTARTEVRY");
    if (uprestart != 0 && fluid_field()->Step() % uprestart == 0)
    {
      fluid_field()->DiscWriter()->WriteVector("slideALE", iprojdisp_);
      fluid_field()->DiscWriter()->WriteVector("slideALEincr", iprojdispinc_);
      slideale_->output_restart(*fluid_field()->DiscWriter());
    }
  }
  ale_field()->Output();

  if (structure_field()->get_constraint_manager()->HaveMonitor())
  {
    structure_field()->get_constraint_manager()->compute_monitor_values(
        structure_field()->Dispnp());
    if (comm_.MyPID() == 0) structure_field()->get_constraint_manager()->PrintMonitorValues();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::OutputLambda()
{
  /* 'lambda_' is only defined on the interface. So, insert 'lambda_' into

   * output or restart data.
   */
  Teuchos::RCP<Epetra_Vector> lambdafull =
      structure_field()->Interface()->InsertFSICondVector(lambda_);
  const int uprestart = timeparams_.get<int>("RESTARTEVRY");
  const int upres = timeparams_.get<int>("RESULTSEVRY");
  if ((uprestart != 0 && fluid_field()->Step() % uprestart == 0) ||
      fluid_field()->Step() % upres == 0)
    structure_field()->DiscWriter()->WriteVector("fsilambda", lambdafull);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::read_restart(int step)
{
  auto input_control_file = GLOBAL::Problem::Instance()->InputControlFile();

  // read Lagrange multiplier
  {
    Teuchos::RCP<Epetra_Vector> lambdafull =
        Teuchos::rcp(new Epetra_Vector(*structure_field()->dof_row_map(), true));
    IO::DiscretizationReader reader =
        IO::DiscretizationReader(structure_field()->discretization(), input_control_file, step);
    reader.ReadVector(lambdafull, "fsilambda");
    lambdaold_ = structure_field()->Interface()->ExtractFSICondVector(lambdafull);
    // Note: the above is normally enough. However, we can use the restart in order to periodically
    // repeat the fsi simulation (see AC-FS3I)
    lambda_ = structure_field()->Interface()->ExtractFSICondVector(lambdafull);
  }

  structure_field()->read_restart(step);
  fluid_field()->read_restart(step);

  SetupSystem();

  if (aleproj_ != INPAR::FSI::ALEprojection_none)
  {
    IO::DiscretizationReader reader =
        IO::DiscretizationReader(fluid_field()->discretization(), input_control_file, step);
    reader.ReadVector(iprojdisp_, "slideALE");
    reader.ReadVector(iprojdispinc_, "slideALEincr");
    slideale_->read_restart(reader);
  }

  ale_field()->read_restart(step);

  SetTimeStep(fluid_field()->Time(), fluid_field()->Step());

  if (aleproj_ != INPAR::FSI::ALEprojection_none)
    slideale_->EvaluateMortar(structure_field()->extract_interface_dispn(), iprojdisp_, *coupsfm_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::recover_lagrange_multiplier()
{
  // get time integration parameter of structural time integrator
  // to enable consistent time integration among the fields
  const double stiparam = structure_field()->TimIntParam();

  // some scaling factors for fluid
  //  const double timescale = fluid_field()->TimeScaling();

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortarp = coupsfm_->GetMortarMatrixP();

  // get the inverted Mortar matrix D^{-1}
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortardinv = coupsfm_->GetMortarMatrixDinv();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (mortarp == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix P.");
  if (mortardinv == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix D^{-1}.");
#endif

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> tmpvec =
      Teuchos::null;  // stores intermediate result of terms (3)-(8)
  Teuchos::RCP<Epetra_Vector> auxvec = Teuchos::null;     // just for convenience
  Teuchos::RCP<Epetra_Vector> auxauxvec = Teuchos::null;  // just for convenience

  /* Recovery of Lagrange multiplier lambda^{n+1} is done by the following
   * condensation expression:
   *
   * lambda^{n+1} =
   *
   * (1)  - stiparam / (1.-stiparam) * lambda^{n}
   *
   * (2)  + 1. / (1.-stiparam) * D^{-T} * tmpvec
   *
   * with tmpvec =
   *
   * (3)    r_{\Gamma}^{S,n+1}
   *
   * (4)  + S_{\Gamma I} * \Delta d_{I}^{S,n+1}
   *
   * (5)  + tau * S_{\Gamma\Gamma} * P * \Delta u_{\Gamma}^{F,n+1}
   *
   * (6)  + dt * S_{\Gamma\Gamma} * P * u_{\Gamma}^n]
   *
   * Remark on term (6):
   * Term (6) has to be considered only in the first Newton iteration.
   * Hence, it will usually not be computed since in general we need more
   * than one nonlinear iteration until convergence.
   *
   * Remarks on all terms:
   * +  Division by -(1.0 - stiparam) will be done in the end
   *    since this is common to all terms
   * +  tau: time scaling factor for interface time integration (tau =
   * 1/fluid_field()->TimeScaling())
   * +  neglecting terms (4)-(6) should not alter the results significantly
   *    since at the end of the time step the solution increments tend to zero.
   *
   *                                                 Matthias Mayr (10/2012)
   */

  // ---------Addressing term (1)
  lambda_->Update(-stiparam, *lambdaold_, 0.0);
  // ---------End of term (1)

  // ---------Addressing term (3)
  Teuchos::RCP<Epetra_Vector> structureresidual = Teuchos::rcp(new Epetra_Vector(
      *structure_field()->Interface()->ExtractFSICondVector(structure_field()->RHS())));
  structureresidual->Scale(-1.0);  // invert sign to obtain residual, not rhs
  tmpvec = Teuchos::rcp(new Epetra_Vector(*structureresidual));
  // ---------End of term (3)

  /* You might want to comment out terms (4) to (6) since they tend to
   * introduce oscillations in the Lagrange multiplier field for certain
   * material properties of the structure.
   *                                                    Matthias Mayr 11/2012
  // ---------Addressing term (4)
  auxvec = Teuchos::rcp(new Epetra_Vector(sgiprev_->RangeMap(),true));
  sgiprev_->Apply(*ddiinc_,*auxvec);
  tmpvec->Update(1.0,*auxvec,1.0);
  // ---------End of term (4)

  // ---------Addressing term (5)
  auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(),true));
  mortarp->Apply(*duginc_,*auxvec);
  auxauxvec = Teuchos::rcp(new Epetra_Vector(sggprev_->RangeMap(),true));
  sggprev_->Apply(*auxvec,*auxauxvec);
  tmpvec->Update(1.0/timescale,*auxauxvec,1.0);
  // ---------End of term (5)

  // ---------Addressing term (6)
  if (firstcall_)
  {
    auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->RangeMap(),true));
    mortarp->Apply(*fluid_field()->extract_interface_veln(),*auxvec);
    auxauxvec = Teuchos::rcp(new Epetra_Vector(sggprev_->RangeMap(),true));
    sggprev_->Apply(*auxvec,*auxauxvec);
    tmpvec->Update(Dt(),*auxauxvec,1.0);
  }
  // ---------End of term (6)
   *
   */

  // ---------Addressing term (2)
  auxvec = Teuchos::rcp(new Epetra_Vector(mortardinv->DomainMap(), true));
  mortardinv->Multiply(true, *tmpvec, *auxvec);
  lambda_->Update(1.0, *auxvec, 1.0);
  // ---------End of term (2)

  // finally, divide by -(1.-stiparam) which is common to all terms
  lambda_->Scale(1. / (1.0 - stiparam));

  /* Finally, the Lagrange multiplier lambda_ is recovered here. It has the
   * unit [N/m^2]. Actual nodal forces are obtained by multiplication with
   * mortar matrices M or D later on.
   */

  //  check_kinematic_constraint();
  //  check_dynamic_equilibrium();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::calculate_interface_energy_increment()
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = structure_field()->TimIntParam();
  const double ftiparam = fluid_field()->TimIntParam();

  // interface traction weighted by time integration factors
  Teuchos::RCP<Epetra_Vector> tractionstructure =
      Teuchos::rcp(new Epetra_Vector(lambda_->Map(), true));
  tractionstructure->Update(stiparam - ftiparam, *lambdaold_, ftiparam - stiparam, *lambda_, 0.0);

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

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::check_kinematic_constraint()
{
  // some scaling factors for fluid
  const double timescale = fluid_field()->TimeScaling();

  // get the Mortar matrices D and M
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortard = coupsfm_->GetMortarMatrixD();
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortarm = coupsfm_->GetMortarMatrixM();

  // get interface displacements and velocities
  Teuchos::RCP<Epetra_Vector> disnp = structure_field()->extract_interface_dispnp();
  Teuchos::RCP<Epetra_Vector> disn = structure_field()->extract_interface_dispn();
  Teuchos::RCP<Epetra_Vector> velnp = fluid_field()->extract_interface_velnp();
  Teuchos::RCP<Epetra_Vector> veln = fluid_field()->extract_interface_veln();

  // prepare vectors for projected interface quantities
  Teuchos::RCP<Epetra_Vector> disnpproj =
      Teuchos::rcp(new Epetra_Vector(mortard->RangeMap(), true));
  Teuchos::RCP<Epetra_Vector> disnproj = Teuchos::rcp(new Epetra_Vector(mortard->RangeMap(), true));
  Teuchos::RCP<Epetra_Vector> velnpproj =
      Teuchos::rcp(new Epetra_Vector(mortarm->RangeMap(), true));
  Teuchos::RCP<Epetra_Vector> velnproj = Teuchos::rcp(new Epetra_Vector(mortarm->RangeMap(), true));

  // projection of interface displacements
  mortard->Apply(*disnp, *disnpproj);
  mortard->Apply(*disn, *disnproj);

  // projection of interface velocities
  mortarm->Apply(*velnp, *velnpproj);
  mortarm->Apply(*veln, *velnproj);

  // calculate violation of kinematic interface constraint
  Teuchos::RCP<Epetra_Vector> violation = Teuchos::rcp(new Epetra_Vector(*disnpproj));
  violation->Update(-1.0, *disnproj, 1.0);
  violation->Update(-1.0 / timescale, *velnpproj, 1.0 / timescale, *velnproj, 1.0);
  violation->Update(-Dt(), *velnproj, 1.0);

  // calculate some norms
  double violationl2 = 0.0;
  double violationinf = 0.0;
  violation->Norm2(&violationl2);
  violation->NormInf(&violationinf);

  // scale L2-Norm with length of vector
  violationl2 /= sqrt(violation->MyLength());

  // output to screen
  std::ios_base::fmtflags flags = utils()->out().flags();

  utils()->out() << std::scientific << "\nViolation of kinematic interface constraint:\n"
                 << "L_2-norm: " << violationl2 << "        L_inf-norm: " << violationinf << "\n";
  utils()->out().flags(flags);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::check_dynamic_equilibrium()
{
  // get the Mortar matrices D and M
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortard = coupsfm_->GetMortarMatrixD();
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortarm = coupsfm_->GetMortarMatrixM();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (mortarm == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix M.");
  if (mortard == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix D.");
#endif

  // auxiliary vectors
  Teuchos::RCP<Epetra_Vector> tractionmaster =
      Teuchos::rcp(new Epetra_Vector(mortarm->DomainMap(), true));
  Teuchos::RCP<Epetra_Vector> tractionslave =
      Teuchos::rcp(new Epetra_Vector(mortard->DomainMap(), true));

  // calculate forces on master and slave side
  mortarm->Multiply(true, *lambda_, *tractionmaster);
  mortard->Multiply(true, *lambda_, *tractionslave);

  // calculate violation of dynamic equilibrium
  Teuchos::RCP<Epetra_Vector> violation = Teuchos::rcp(new Epetra_Vector(*tractionmaster));
  violation->Update(-1.0, *tractionslave, 1.0);

  // calculate some norms
  double violationl2 = 0.0;
  double violationinf = 0.0;
  violation->Norm2(&violationl2);
  violation->NormInf(&violationinf);

  // scale L2-Norm with sqrt of length of interface vector
  violationl2 /= sqrt(structure_field()->Interface()->FSICondMap()->NumGlobalElements());

  // output to screen
  std::ios_base::fmtflags flags = utils()->out().flags();

  utils()->out() << std::scientific << "\nViolation of dynamic interface equilibrium:\n"
                 << "L_2-norm: " << violationl2 << "        L_inf-norm: " << violationinf << "\n";
  utils()->out().flags(flags);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::combine_field_vectors(Epetra_Vector& v,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av, const bool slave_vectors_contain_interface_dofs)
{
  if (slave_vectors_contain_interface_dofs)
  {
    // extract inner DOFs from slave vectors
    Teuchos::RCP<Epetra_Vector> sov = structure_field()->Interface()->ExtractOtherVector(sv);
    Teuchos::RCP<Epetra_Vector> aov = FsiAleField()->FsiInterface()->ExtractOtherVector(av);

    // put them together
    FSI::Monolithic::combine_field_vectors(v, sov, fv, aov);
  }
  else
    FSI::Monolithic::combine_field_vectors(v, sv, fv, av);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double FSI::SlidingMonolithicStructureSplit::select_dt_error_based() const
{
  // get time step size suggestions based on some error norms
  const double dtfl = get_ada_fl_dt();               // based on all fluid DOFs
  const double dtflfsi = get_ada_fl_fsi_dt();        // based on fluid FSI DOFs
  const double dtstrinner = get_ada_str_inner_dt();  // based on inner structural DOFs

  double dt = Dt();

  // select time step size based on error estimation
  if (is_ada_structure() and is_ada_fluid())
    dt = std::min(std::min(dtfl, dtflfsi), dtstrinner);
  else if (is_ada_structure() and (not is_ada_fluid()))
    dt = dtstrinner;
  else if ((not is_ada_structure()) and is_ada_fluid())
    dt = std::min(dtfl, dtflfsi);
  else
  {
    // no change in time step size based on structure or fluid field error estimation
  }

  return dt;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool FSI::SlidingMonolithicStructureSplit::set_accepted() const
{
  // get error norms
  const double flnorm = get_ada_flnorm();               // based on all fluid DOFs
  const double flfsinorm = get_ada_fl_fs_inorm();       // based on fluid FSI DOFs
  const double strinnernorm = get_ada_str_innernorm();  // based on inner structural DOFs

  bool accepted = std::max(flnorm, flfsinorm) < errtolfl_ && strinnernorm < errtolstr_;

  // in case error estimation in the fluid field is turned off:
  if (not is_ada_fluid()) accepted = strinnernorm < errtolstr_;

  // in case error estimation in the structure field is turned off:
  if (not is_ada_structure()) accepted = std::max(flnorm, flfsinorm) < errtolfl_;

  // no error based time adaptivity
  if ((not is_ada_structure()) and (not is_ada_fluid())) accepted = true;

  return accepted;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::create_system_matrix()
{
  FSI::BlockMonolithic::create_system_matrix(systemmatrix_, true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::SlidingMonolithicStructureSplit::create_node_owner_relationship(
    std::map<int, int>* nodeOwner, std::map<int, std::list<int>>* inverseNodeOwner,
    std::map<int, CORE::Nodes::Node*>* structurenodesPtr,
    std::map<int, CORE::Nodes::Node*>* fluidgnodesPtr,
    Teuchos::RCP<DRT::Discretization> structuredis, Teuchos::RCP<DRT::Discretization> fluiddis,
    const INPAR::FSI::Redistribute domain)
{
  /*******************************************/
  /* distribute masternodes to future owners */
  /*******************************************/
  /* Idea:
   * - the P matrix maps dofs from fluid to structure (fluidsplit) and from structure to fluid
   * (structuresplit)
   * - find "neighboring dofs" via entries in row of P matrix
   * - get owner of the node of the corresponding slave dof and global node
   *   id of the node of the corresponding master dof
   * - save this pair in map nodeOwner: global id of master node -> desired owner of this node
   */

  int numproc = comm_.NumProc();
  int myrank = comm_.MyPID();

  // get P matrix
  Teuchos::RCP<CORE::LINALG::SparseMatrix> P =
      coupsfm_->GetMortarMatrixP();  // P matrix that couples structure dofs with fluid dofs
  Teuchos::RCP<Epetra_CrsMatrix> P_;
  P_ = P->EpetraMatrix();
  const Epetra_Map& P_Map = P_->RowMap();  // maps fluid dofs to procs

  int NumMyElements = P_Map.NumMyElements();
  int numFluidDofs = P_->NumGlobalCols();  // number of related fluid dofs on interface

  int NumMaxElements;
  comm_.MaxAll(&NumMyElements, &NumMaxElements, 1);

  int skip = GLOBAL::Problem::Instance()->NDim();  // Only evaluate one dof per node, skip the other
                                                   // dofs related to the node. It is nsd_=skip.

  int re[2];   // return: re[0] = global node id, re[1] = owner of node
  int flnode;  // fluid
  int flowner;
  int stnode;  // structure
  int stowner = -1;


  std::map<int, int>::iterator nodeOwnerIt;

  // loop over fluid dofs, i.e. rows in the P matrix
  /* We loop over NumMaxElements because all procs have to stay in the loop until
   * the last proc checked each of its rows.
   */


  for (int stdofLID = 0; stdofLID < NumMaxElements; stdofLID = stdofLID + skip)
  {
    flnode = -1;
    flowner = -1;
    stnode = -1;
    stowner = -1;

    int NumEntries = 0;
    std::vector<double> Values(numFluidDofs);
    std::vector<int> fldofGID(numFluidDofs);

    if (NumMyElements > 0)
    {
      int stdofGID = P_Map.GID(stdofLID);  // gid of structure dof
      // find related node and owner to stdofGID
      find_node_related_to_dof(structurenodesPtr, stdofGID, structuredis, re);  // fluid
      stnode = re[0];
      stowner = re[1];

      P_->ExtractGlobalRowCopy(stdofGID, numFluidDofs, NumEntries, Values.data(), fldofGID.data());
    }

    // Loop over related structure dofs and get related nodes.

    int maxNumEntries;
    comm_.MaxAll(&NumEntries, &maxNumEntries, 1);

    for (int j = 0; j < maxNumEntries; ++j)
    {
      if (j < NumEntries)
      {
        find_node_related_to_dof(fluidgnodesPtr, fldofGID[j], fluiddis, re);
        flnode = re[0];
        flowner = re[1];
      }

      // A processor can only find a neighbored node if it belongs to its domain (which in general
      // is not the case). Therefore, we have to search in the other processors' domains.

      int copy;
      int foundNode = -2;
      int sendNode = -2;
      int saveNode = -2;
      int foundOwner;
      int sendOwner;
      int saveOwner = -2;
      int dofid;

      for (int proc = 0; proc < numproc; ++proc)
      {
        copy = flnode;
        comm_.Broadcast(&copy, 1, proc);  // send nodeid to check if the processor needs help to
                                          // find its neighbor, code: copy=-2
        if (copy == -2)
        {
          dofid = fldofGID[j];
          comm_.Broadcast(&dofid, 1, proc);
          find_node_related_to_dof(fluidgnodesPtr, dofid, fluiddis,
              re);  // let each processor look for the node related to gstid
          foundNode = re[0];
          foundOwner = re[1];
          for (int j = 0; j < numproc; ++j)
          {
            sendNode = foundNode;
            sendOwner = foundOwner;
            comm_.Broadcast(&sendNode, 1, j);
            comm_.Broadcast(&sendOwner, 1, j);
            if (sendNode != -2)
            {  // check which processor found the node
              saveNode = sendNode;
              saveOwner = sendOwner;
              break;
            }
          }
          if (myrank == proc && saveNode != -2)
          {  // save the nodegid on the respective processor
            stnode = saveNode;
            stowner = saveOwner;
          }
        }
      }

      if (domain == INPAR::FSI::Redistribute_structure && stnode != -1)
      {  // map structure nodes to fluid owners
        (*nodeOwner)[stnode] = flowner;
      }
      else if (domain == INPAR::FSI::Redistribute_fluid)  // map fluid nodes to structure owners
        (*nodeOwner)[flnode] = stowner;

    }  // for (int j = 0; j < maxNumEntries; ++j){

    if (NumMyElements != 0)
    {
      NumMyElements -= skip;
    }

  }  // end loop over structure dofs


  //    std::map<int,int>::iterator nodeOwnerPrint;
  //    for (nodeOwnerPrint = nodeOwner.begin(); nodeOwnerPrint != nodeOwner.end();
  //    ++nodeOwnerPrint){
  //      std::cout<<"\nNode: "<<nodeOwnerPrint->first<<" Owner: "<<nodeOwnerPrint->second<<" I am
  //      proc "<<myrank;
  //    }


  // If the structure is redistributed, it might occur that one node is contained several times in
  // the list because several procs might have introduced it.
  for (int proc = 0; proc < numproc; ++proc)
  {
    // std::map<int,int>::iterator nodeOwnerIt;   already declared
    nodeOwnerIt = nodeOwner->begin();
    int nodeSize = (int)nodeOwner->size();
    comm_.Broadcast(&nodeSize, 1, proc);
    int node;
    for (int i = 0; i < nodeSize; ++i)
    {
      node = nodeOwnerIt->first;
      comm_.Broadcast(&node, 1, proc);
      if (myrank > proc)
      {
        try
        {
          nodeOwner->at(node);
          nodeOwner->erase(node);
        }
        catch (std::exception& exc)
        {
        }
      }
      if (myrank == proc)
      {
        nodeOwnerIt++;
      }
    }
  }

  //  std::map<int,int>::iterator nodeOwnerPrint;
  //  for (nodeOwnerPrint = nodeOwner.begin(); nodeOwnerPrint != nodeOwner.end(); ++nodeOwnerPrint){
  //    std::cout<<"\nNode: "<<nodeOwnerPrint->first<<" Owner: "<<nodeOwnerPrint->second<<" I am
  //    proc "<<comm_.MyPID();
  //  }

  std::list<int>::iterator listIt;
  int sendPair[2];
  nodeOwnerIt = nodeOwner->begin();  // already declared

  for (int proc = 0; proc < numproc; ++proc)
  {
    int numNodes = nodeOwner->size();
    comm_.Broadcast(&numNodes, 1, proc);

    for (int i = 0; i < numNodes; ++i)
    {
      if (myrank == proc)
      {
        sendPair[0] = nodeOwnerIt->first;
        sendPair[1] = nodeOwnerIt->second;
      }
      comm_.Broadcast(sendPair, 2, proc);
      listIt = (*inverseNodeOwner)[sendPair[1]].begin();
      (*inverseNodeOwner)[sendPair[1]].insert(listIt, sendPair[0]);

      if (myrank == proc) nodeOwnerIt++;
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
