/*--------------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problem with sliding grids using a monolithic scheme
with condensed fluid interface velocities

\level 2

*/
/*--------------------------------------------------------------------------*/


#include "4C_fsi_slidingmonolithic_fluidsplit.hpp"

#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_constraint_manager.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_fsi_statustest.hpp"
#include "4C_fsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_structure_aux.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
FSI::SlidingMonolithicFluidSplit::SlidingMonolithicFluidSplit(
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
  intersectionmaps.push_back(fluid_field()->get_dbc_map_extractor()->cond_map());
  intersectionmaps.push_back(fluid_field()->interface()->fsi_cond_map());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      Core::LinAlg::MultiMapExtractor::intersect_maps(intersectionmaps);

  // Check whether the intersection is empty
  if (intersectionmap->NumGlobalElements() != 0)
  {
    //    std::cout << "Slave interface nodes with Dirichlet boundary condition "
    //              "(input file numbering):" << std::endl;
    //    for (int i=0; i < (int)fluid_field()->discretization()->NumMyRowNodes(); i++)
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

  notsetup_ = true;

  coupsfm_ = Teuchos::rcp(new Coupling::Adapter::CouplingMortar(
      Global::Problem::instance()->n_dim(), Global::Problem::instance()->mortar_coupling_params(),
      Global::Problem::instance()->contact_dynamic_params(),
      Global::Problem::instance()->spatial_approximation_type()));
  fscoupfa_ = Teuchos::rcp(new Coupling::Adapter::Coupling());

  aigtransform_ = Teuchos::rcp(new Coupling::Adapter::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new Coupling::Adapter::MatrixColTransform);

  coupsfm_ = Teuchos::rcp(new Coupling::Adapter::CouplingMortar(
      Global::Problem::instance()->n_dim(), Global::Problem::instance()->mortar_coupling_params(),
      Global::Problem::instance()->contact_dynamic_params(),
      Global::Problem::instance()->spatial_approximation_type()));
  fscoupfa_ = Teuchos::rcp(new Coupling::Adapter::Coupling());

  aigtransform_ = Teuchos::rcp(new Coupling::Adapter::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new Coupling::Adapter::MatrixColTransform);

  // Recovery of Lagrange multiplier happens on fluid field
  set_lambda();

  fmgiprev_ = Teuchos::null;
  fmgicur_ = Teuchos::null;
  fmggprev_ = Teuchos::null;
  fmggcur_ = Teuchos::null;
  fgiprev_ = Teuchos::null;
  fgicur_ = Teuchos::null;
  fggprev_ = Teuchos::null;
  fggcur_ = Teuchos::null;

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
void FSI::SlidingMonolithicFluidSplit::set_lambda()
{
  lambda_ = Teuchos::rcp(new Epetra_Vector(*fluid_field()->interface()->fsi_cond_map(), true));
  lambdaold_ = Teuchos::rcp(new Epetra_Vector(*fluid_field()->interface()->fsi_cond_map(), true));

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::setup_system()
{
  if (notsetup_)
  {
    const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
    const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
    linearsolverstrategy_ =
        Core::UTILS::integral_value<Inpar::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

    aleproj_ = Core::UTILS::integral_value<Inpar::FSI::SlideALEProj>(fsidyn, "SLIDEALEPROJ");

    set_default_parameters(fsidyn, nox_parameter_list());

    // we use non-matching meshes at the interface
    // mortar with: structure = master, fluid = slave

    const int ndim = Global::Problem::instance()->n_dim();

    // get coupling objects
    Coupling::Adapter::Coupling& icoupfa = interface_fluid_ale_coupling();

    /* structure to fluid
     * coupling condition at the fsi interface:
     * displacements (=number spacial dimensions) are coupled
     * e.g.: 3D: coupleddof = [1, 1, 1]
     */
    std::vector<int> coupleddof(ndim, 1);

    coupsfm_->setup(structure_field()->discretization(), fluid_field()->discretization(),
        ale_field()->write_access_discretization(), coupleddof, "FSICoupling", comm_,
        Global::Problem::instance()->function_manager(),
        Global::Problem::instance()->binning_strategy_params(),
        Global::Problem::instance()->discretization_map(),
        Global::Problem::instance()->output_control_file(),
        Global::Problem::instance()->spatial_approximation_type(), true);

    // fluid to ale at the interface

    icoupfa.setup_condition_coupling(*fluid_field()->discretization(),
        fluid_field()->interface()->fsi_cond_map(), *ale_field()->discretization(),
        ale_field()->interface()->fsi_cond_map(), "FSICoupling", ndim);

    // we might have a free surface
    if (fluid_field()->interface()->fs_cond_relevant())
    {
      fscoupfa_->setup_condition_coupling(*fluid_field()->discretization(),
          fluid_field()->interface()->fs_cond_map(), *ale_field()->discretization(),
          ale_field()->interface()->fs_cond_map(), "FREESURFCoupling", ndim);
    }

    Coupling::Adapter::Coupling& coupfa = fluid_ale_coupling();

    // the fluid-ale coupling always matches
    const Epetra_Map* fluidnodemap = fluid_field()->discretization()->node_row_map();
    const Epetra_Map* alenodemap = ale_field()->discretization()->node_row_map();

    coupfa.setup_coupling(*fluid_field()->discretization(), *ale_field()->discretization(),
        *fluidnodemap, *alenodemap, ndim);

    fluid_field()->set_mesh_map(coupfa.master_dof_map());

    // create combined map
    create_combined_dof_row_map();

    /*------------------------------------------------------------------------*/
    // Switch fluid to interface split block matrix
    fluid_field()->use_block_matrix(true);

    // build ale system matrix in splitted system
    ale_field()->create_system_matrix(ale_field()->interface());

    aleresidual_ = Teuchos::rcp(new Epetra_Vector(*fsi_ale_field()->fsi_interface()->other_map()));

    // -------------------------------------------------------------------------
    // Build the global Dirichlet map extractor
    setup_dbc_map_extractor();
    // -------------------------------------------------------------------------

    // enable debugging
    if (Core::UTILS::integral_value<int>(fsidyn, "DEBUGOUTPUT") & 2)
    {
      pcdbg_ = Teuchos::rcp(new UTILS::MonolithicDebugWriter(*this));
    }

    create_system_matrix();

    if (aleproj_ != Inpar::FSI::ALEprojection_none)
    {
      // set up sliding ale utils
      slideale_ = Teuchos::rcp(new FSI::UTILS::SlideAleUtils(structure_field()->discretization(),
          fluid_field()->discretization(), *coupsfm_, true, aleproj_));

      iprojdispinc_ = Teuchos::rcp(new Epetra_Vector(*coupsfm_->slave_dof_map(), true));
      iprojdisp_ = Teuchos::rcp(new Epetra_Vector(*coupsfm_->slave_dof_map(), true));
    }
    notsetup_ = false;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::create_combined_dof_row_map()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(structure_field()->dof_row_map());
  vecSpaces.push_back(fluid_field()->dof_row_map());
  vecSpaces.push_back(fsi_ale_field()->fsi_interface()->other_map());

  if (vecSpaces[1]->NumGlobalElements() == 0)
    FOUR_C_THROW("No inner fluid equations. Splitting not possible.");

  set_dof_row_maps(vecSpaces);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::setup_dbc_map_extractor()
{
  /* Dirichlet maps for structure and fluid do not intersect with interface map.
   * ALE Dirichlet map might intersect with interface map, but ALE interface
   * DOFs are not part of the final system of equations. Hence, we just need the
   * intersection of inner ALE DOFs with Dirichlet ALE DOFs.
   */
  std::vector<Teuchos::RCP<const Epetra_Map>> aleintersectionmaps;
  aleintersectionmaps.push_back(ale_field()->get_dbc_map_extractor()->cond_map());
  aleintersectionmaps.push_back(fsi_ale_field()->fsi_interface()->other_map());
  Teuchos::RCP<Epetra_Map> aleintersectionmap =
      Core::LinAlg::MultiMapExtractor::intersect_maps(aleintersectionmaps);

  // Merge Dirichlet maps of structure, fluid and ALE to global FSI Dirichlet map
  std::vector<Teuchos::RCP<const Epetra_Map>> dbcmaps;
  dbcmaps.push_back(structure_field()->get_dbc_map_extractor()->cond_map());
  dbcmaps.push_back(fluid_field()->get_dbc_map_extractor()->cond_map());
  dbcmaps.push_back(aleintersectionmap);
  Teuchos::RCP<const Epetra_Map> dbcmap = Core::LinAlg::MultiMapExtractor::merge_maps(dbcmaps);

  // Finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new Core::LinAlg::MapExtractor(*dof_row_map(), dbcmap, true));
  if (dbcmaps_ == Teuchos::null) FOUR_C_THROW("Creation of FSI Dirichlet map extractor failed.");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> FSI::SlidingMonolithicFluidSplit::system_matrix()
    const
{
  return systemmatrix_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::setup_rhs_residual(Epetra_Vector& f)
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = structure_field()->tim_int_param();
  const double ftiparam = fluid_field()->tim_int_param();

  // some scaling factors for fluid
  const double fluidscale = fluid_field()->residual_scaling();

  // get the Mortar matrix M
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortarp = coupsfm_->get_mortar_matrix_p();

  // get single field residuals
  Teuchos::RCP<const Epetra_Vector> sv = Teuchos::rcp(new Epetra_Vector(*structure_field()->rhs()));
  Teuchos::RCP<const Epetra_Vector> fv = Teuchos::rcp(new Epetra_Vector(*fluid_field()->rhs()));
  Teuchos::RCP<const Epetra_Vector> av = Teuchos::rcp(new Epetra_Vector(*ale_field()->rhs()));

  // extract only inner DOFs from fluid (=slave) and ALE field
  Teuchos::RCP<Epetra_Vector> fov = fsi_fluid_field()->fsi_interface()->extract_other_vector(fv);
  fov = fsi_fluid_field()->fsi_interface()->insert_other_vector(fov);
  Teuchos::RCP<const Epetra_Vector> aov =
      fsi_ale_field()->fsi_interface()->extract_other_vector(av);

  /* add fluid interface residual to structure interface residual considering
   * temporal scaling
   */
  Teuchos::RCP<Epetra_Vector> fcv = fluid_field()->interface()->extract_fsi_cond_vector(fv);
  Teuchos::RCP<Epetra_Vector> scv =
      Core::LinAlg::create_vector(*structure_field()->interface()->fsi_cond_map(), true);
  mortarp->multiply(true, *fcv, *scv);
  Teuchos::RCP<Epetra_Vector> modsv = structure_field()->interface()->insert_fsi_cond_vector(scv);
  modsv->Update(1.0, *sv, (1.0 - stiparam) / (1.0 - ftiparam) * fluidscale);

  if (structure_field()->get_stc_algo() == Inpar::Solid::stc_currsym)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> stcmat = structure_field()->get_stc_mat();
    stcmat->multiply(true, *modsv, *modsv);
  }

  // put the single field residuals together
  FSI::Monolithic::combine_field_vectors(f, modsv, fov, aov);

  // add additional ale residual
  extractor().add_vector(*aleresidual_, 2, f);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::setup_rhs_lambda(Epetra_Vector& f)
{
  if (lambdaold_ != Teuchos::null)
  {
    // get time integration parameters of structure and fluid time integrators
    // to enable consistent time integration among the fields
    const double stiparam = structure_field()->tim_int_param();
    const double ftiparam = fluid_field()->tim_int_param();

    // get the Mortar matrix M
    const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortarm = coupsfm_->get_mortar_matrix_m();

    /* project Lagrange multiplier field onto the master interface DOFs and
     * consider temporal scaling */
    Teuchos::RCP<Epetra_Vector> lambda =
        Teuchos::rcp(new Epetra_Vector(mortarm->domain_map(), true));
    mortarm->multiply(true, *lambdaold_, *lambda);
    Teuchos::RCP<Epetra_Vector> lambdafull =
        structure_field()->interface()->insert_fsi_cond_vector(lambda);
    lambdafull->Scale(stiparam - (ftiparam * (1.0 - stiparam)) / (1.0 - ftiparam));

    // add Lagrange multiplier
    extractor().add_vector(*lambdafull, 0, f);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::setup_rhs_firstiter(Epetra_Vector& f)
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = structure_field()->tim_int_param();
  const double ftiparam = fluid_field()->tim_int_param();

  // some scaling factors for fluid
  const double timescale = fluid_field()->time_scaling();
  const double scale = fluid_field()->residual_scaling();

  // old interface velocity of fluid field
  const Teuchos::RCP<const Epetra_Vector> fveln = fluid_field()->extract_interface_veln();

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<const Core::LinAlg::SparseMatrix> mortarp = coupsfm_->get_mortar_matrix_p();

  // get fluid matrix
  const Teuchos::RCP<const Core::LinAlg::BlockSparseMatrixBase> blockf =
      fluid_field()->block_system_matrix();

  // get fluid shape derivatives matrix
  const Teuchos::RCP<const Core::LinAlg::BlockSparseMatrixBase> mmm =
      fluid_field()->shape_derivatives();

  // get ale matrix
  const Teuchos::RCP<const Core::LinAlg::BlockSparseMatrixBase> blocka =
      ale_field()->block_system_matrix();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (mortarp == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix P.");
  if (blockf == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to fluid block matrix.");
  if (blocka == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to ale block matrix.");
#endif

  // extract fluid and ale submatrices
  const Core::LinAlg::SparseMatrix& fig = blockf->matrix(0, 1);  // F_{I\Gamma}
  const Core::LinAlg::SparseMatrix& fgg = blockf->matrix(1, 1);  // F_{\Gamma\Gamma}
  const Core::LinAlg::SparseMatrix& aig = blocka->matrix(0, 1);  // A_{I\Gamma}

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;     // right hand side of single set of DOFs
  Teuchos::RCP<Epetra_Vector> auxvec = Teuchos::null;  // just for convenience
  Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::null;  // just for convenience

  /* Different contributions/terms to the rhs are separated by the following
   * comment line */
  // ---------- interface structure DOFs
  /* The following terms are added to the interface structure DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  + (1-stiparam)/(1-ftiparam) * dt / tau * P^{T} * F_{\Gamma\Gamma} * u^{n}_{\Gamma}
   *
   * (2)  - (1-stiparam)/(1-ftiparam) / tau * P^{T} * F_{\Gamma\Gamma} * P * \Delta d_{\Gamma,p}
   *
   * (3)  - (1-stiparam)/(1-ftiparam) * P^{T} * F^{G}_{\Gamma\Gamma} * P * \Delta d_{\Gamma,p}
   *
   * Remarks on all terms:
   * +  tau: time scaling factor for interface time integration (tau =
   * 1/fluid_field()->TimeScaling())
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(mortarp->domain_map(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(fgg.row_map(), true));

  fgg.Apply(*fveln, *auxvec);
  mortarp->multiply(true, *auxvec, *rhs);

  if (structure_field()->get_stc_algo() == Inpar::Solid::stc_currsym)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> stcmat = structure_field()->get_stc_mat();
    stcmat->multiply(true, *rhs, *rhs);
  }

  rhs->Scale(scale * (1. - stiparam) / (1. - ftiparam) * dt() * timescale);
  rhs = structure_field()->interface()->insert_fsi_cond_vector(rhs);

  extractor().add_vector(*rhs, 0, f);
  // ----------end of term 1

  // ----------addressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(mortarp->domain_map(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(fgg.range_map(), true));
  tmpvec = Teuchos::rcp(new Epetra_Vector(mortarp->range_map(), true));

  mortarp->Apply(*ddgpred_, *tmpvec);
  fgg.Apply(*tmpvec, *auxvec);
  mortarp->multiply(true, *auxvec, *rhs);

  rhs->Scale(-scale * (1. - stiparam) / (1. - ftiparam) * timescale);
  rhs = structure_field()->interface()->insert_fsi_cond_vector(rhs);

  extractor().add_vector(*rhs, 0, f);
  // ----------end of term 2

  // ----------addressing term 3
  if (mmm != Teuchos::null)
  {
    // extract F^{G}_{\Gamma\Gamma}
    const Core::LinAlg::SparseMatrix& fmgg = mmm->matrix(1, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(mortarp->domain_map(), true));
    auxvec = Teuchos::rcp(new Epetra_Vector(fmgg.range_map(), true));
    tmpvec = Teuchos::rcp(new Epetra_Vector(mortarp->range_map(), true));

    mortarp->Apply(*ddgpred_, *tmpvec);
    fmgg.Apply(*tmpvec, *auxvec);
    mortarp->multiply(true, *auxvec, *rhs);

    rhs->Scale(-(1. - stiparam) / (1. - ftiparam));
    rhs = structure_field()->interface()->insert_fsi_cond_vector(rhs);

    extractor().add_vector(*rhs, 0, f);
  }
  // ----------end of term 3
  // ----------end of interface structure DOFs

  // ---------- inner fluid DOFs
  /* The following terms are added to the inner fluid DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  + dt / tau * F_{I \Gamma} * u^{n}_{\Gamma}
   *
   * (2)  - 1 / tau * F_{I \Gamma} * P * \Delta d_{\Gamma,p}
   *
   * (3)  - F^{G}_{I \Gamma} * P * \Delta d_{\Gamma,p}
   *
   * Remarks on all terms:
   * +  tau: time scaling factor for interface time integration
   *         (tau = 1/fluid_field()->TimeScaling())
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(fig.row_map(), true));

  fig.Apply(*fveln, *rhs);

  rhs->Scale(dt() * timescale);

  rhs = fsi_fluid_field()->fsi_interface()->insert_other_vector(rhs);

  extractor().add_vector(*rhs, 1, f);
  // ----------end of term 1

  // ----------addressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(fig.range_map(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->range_map(), true));

  mortarp->Apply(*ddgpred_, *auxvec);
  fig.Apply(*auxvec, *rhs);

  rhs->Scale(-timescale);

  rhs = fsi_fluid_field()->fsi_interface()->insert_other_vector(rhs);

  extractor().add_vector(*rhs, 1, f);
  // ----------end of term 2

  // ----------addressing term 3
  if (mmm != Teuchos::null)
  {
    // extract F^{G}_{I \Gamma}
    const Core::LinAlg::SparseMatrix& fmig = mmm->matrix(0, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fmig.range_map(), true));
    auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->range_map(), true));

    mortarp->Apply(*ddgpred_, *auxvec);
    fmig.Apply(*auxvec, *rhs);

    rhs->Scale(-1.);

    rhs = fsi_fluid_field()->fsi_interface()->insert_other_vector(rhs);

    extractor().add_vector(*rhs, 1, f);
  }
  // ----------end of term 3
  // ----------end of inner fluid DOFs

  // ---------- inner ALE DOFs
  /* The following terms are added to the inner ALE DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  - A_{I \Gamma} * P * \Delta d_{\Gamma,p}
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(aig.range_map(), true));
  auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->range_map(), true));

  mortarp->Apply(*ddgpred_, *auxvec);
  aig.Apply(*fluid_to_ale_interface(auxvec), *rhs);

  rhs->Scale(-1.0);

  extractor().add_vector(*rhs, 2, f);
  // ----------end of term 1
  // ----------end of inner ALE DOFs

  // only if relative movement between ale and structure is possible
  if (aleproj_ != Inpar::FSI::ALEprojection_none)
  {
    // get block ale matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> a = ale_field()->block_system_matrix();
    if (a == Teuchos::null)
    {
      FOUR_C_THROW("expect ale block matrix");
    }

    rhs = Teuchos::rcp(new Epetra_Vector(a->matrix(0, 1).row_map()));

    a->matrix(0, 1).Apply(*fluid_to_ale_interface(iprojdispinc_), *rhs);

    extractor().add_vector(*rhs, 2, f);

    // get fluid shape derivative matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> mmm = fluid_field()->shape_derivatives();
    if (mmm != Teuchos::null)
    {
      // extract submatrices
      Core::LinAlg::SparseMatrix& fmig = mmm->matrix(0, 1);
      Core::LinAlg::SparseMatrix& fmgg = mmm->matrix(1, 1);

      rhs = Teuchos::rcp(new Epetra_Vector(fmgg.row_map()));

      fmgg.Apply(*iprojdispinc_, *rhs);

      Teuchos::RCP<Epetra_Vector> tmprhs = Teuchos::rcp(new Epetra_Vector(mortarp->domain_map()));
      mortarp->multiply(true, *rhs, *tmprhs);

      rhs = structure_field()->interface()->insert_fsi_cond_vector(tmprhs);

      auto zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(), true));
      Core::LinAlg::apply_dirichlet_to_system(
          *rhs, *zeros, *(structure_field()->get_dbc_map_extractor()->cond_map()));

      if (structure_field()->get_stc_algo() == Inpar::Solid::stc_currsym)
      {
        Teuchos::RCP<Core::LinAlg::SparseMatrix> stcmat = structure_field()->get_stc_mat();
        stcmat->multiply(true, *rhs, *rhs);
      }

      rhs->Scale(dt() * timescale * (1. - stiparam) / (1. - ftiparam));
      extractor().add_vector(*rhs, 0, f);

      rhs = Teuchos::rcp(new Epetra_Vector(fmig.row_map()));

      fmig.Apply(*iprojdispinc_, *rhs);

      rhs = fsi_fluid_field()->fsi_interface()->insert_other_vector(rhs);

      zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(), true));
      Core::LinAlg::apply_dirichlet_to_system(
          *rhs, *zeros, *(structure_field()->get_dbc_map_extractor()->cond_map()));

      rhs->Scale(-timescale * dt());

      extractor().add_vector(*rhs, 1, f);
    }
  }

  // Reset quantities of previous iteration step since they still store values from the last time
  // step
  ddginc_ = Core::LinAlg::create_vector(*structure_field()->interface()->fsi_cond_map(), true);
  duiinc_ = Core::LinAlg::create_vector(*fsi_fluid_field()->fsi_interface()->other_map(), true);
  veliprev_ = Teuchos::null;
  velgprev_ = Teuchos::null;
  fgicur_ = Teuchos::null;
  fggcur_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::setup_system_matrix");

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortarp = coupsfm_->get_mortar_matrix_p();

  // get info about STC feature
  Inpar::Solid::StcScale stcalgo = structure_field()->get_stc_algo();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> stcmat = Teuchos::null;
  if (stcalgo != Inpar::Solid::stc_none) stcmat = structure_field()->get_stc_mat();

  const Coupling::Adapter::Coupling& coupfa = fluid_ale_coupling();

  // get single field block matrices
  Teuchos::RCP<Core::LinAlg::SparseMatrix> s =
      structure_field()->system_matrix();  // can't be 'const' --> is modified by STC
  const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> f = fluid_field()->block_system_matrix();
  const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> a = ale_field()->block_system_matrix();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // check whether allocation was successful
  if (mortarp == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix P.");
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

  // some checks whether maps for matrix-matrix-multiplication do really match
  if (!f->matrix(0, 1).domain_map().PointSameAs(mortarp->range_map()))
    FOUR_C_THROW("Maps do not match.");
  if (!f->matrix(1, 0).range_map().PointSameAs(mortarp->range_map()))
    FOUR_C_THROW("Maps do not match.");
  if (!f->matrix(1, 1).domain_map().PointSameAs(mortarp->range_map()))
    FOUR_C_THROW("Maps do not match.");
#endif

  // extract submatrices
  Core::LinAlg::SparseMatrix& aii = a->matrix(0, 0);  // A_{II}
  Core::LinAlg::SparseMatrix& aig = a->matrix(0, 1);  // A_{I\Gamma}
  Core::LinAlg::SparseMatrix& fii = f->matrix(0, 0);  // F_{II}

  // scaling factors for fluid
  const double scale = fluid_field()->residual_scaling();
  const double timescale = fluid_field()->time_scaling();

  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = structure_field()->tim_int_param();
  const double ftiparam = fluid_field()->tim_int_param();

  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  s->un_complete();

  // ---------------------------------------------------------------------------
  // BEGIN building the global 4x4 system matrix
  // ---------------------------------------------------------------------------
  // Contributions to blocks in system matrix are listed separately.
  // Block numbering in comments ranges from (1,1) to (4,4).

  // ---------Addressing contribution to block (2,2)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> fgg =
      matrix_multiply(f->matrix(1, 1), false, *mortarp, false, false, false, true);
  fgg = matrix_multiply(*mortarp, true, *fgg, false, false, false, true);

  s->add(*fgg, false, scale * timescale * (1. - stiparam) / (1. - ftiparam), 1.0);

  // ---------Addressing contribution to block (2,3)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> fgi =
      matrix_multiply(*mortarp, true, f->matrix(1, 0), false, false, false, true);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> lfgi =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(s->row_map(), 81, false));

  lfgi->add(*fgi, false, scale, 0.0);
  lfgi->complete(fgi->domain_map(), s->range_map());

  if (stcalgo == Inpar::Solid::stc_currsym)
    lfgi = Core::LinAlg::matrix_multiply(*stcmat, true, *lfgi, false, true, true, true);

  mat.matrix(0, 1).un_complete();
  mat.matrix(0, 1).add(*lfgi, false, (1. - stiparam) / (1. - ftiparam), 0.0);

  // ---------Addressing contribution to block (3,2)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> fig =
      matrix_multiply(f->matrix(0, 1), false, *mortarp, false, false, false, true);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> lfig =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(fig->row_map(), 81, false));

  lfig->add(*fig, false, timescale, 0.0);
  lfig->complete(s->domain_map(), fig->range_map());

  if (stcalgo != Inpar::Solid::stc_none)
  {
    lfig = Core::LinAlg::matrix_multiply(*lfig, false, *stcmat, false, false, false, true);
  }

  mat.matrix(1, 0).un_complete();
  mat.matrix(1, 0).add(*lfig, false, 1., 0.0);

  // ---------Addressing contribution to block (3,3)
  mat.matrix(1, 1).un_complete();
  mat.matrix(1, 1).add(fii, false, 1., 0.0);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> eye =
      Core::LinAlg::create_identity_matrix(*fluid_field()->interface()->fsi_cond_map());
  mat.matrix(1, 1).add(*eye, false, 1., 1.0);

  // ---------Addressing contribution to block (4,2)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> laig =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(aii.row_map(), 81, false));
  (*aigtransform_)(a->full_row_map(), a->full_col_map(), aig, 1.,
      Coupling::Adapter::CouplingSlaveConverter(interface_fluid_ale_coupling()), *laig);

  laig->complete(f->matrix(1, 1).domain_map(), aii.range_map());
  Teuchos::RCP<Core::LinAlg::SparseMatrix> llaig =
      matrix_multiply(*laig, false, *mortarp, false, false, false, true);
  laig = Teuchos::rcp(new Core::LinAlg::SparseMatrix(llaig->row_map(), 81, false));

  laig->add(*llaig, false, 1.0, 0.0);
  laig->complete(s->domain_map(), llaig->range_map());

  if (stcalgo != Inpar::Solid::stc_none)
  {
    laig = Core::LinAlg::matrix_multiply(*laig, false, *stcmat, false, false, false, true);
  }

  mat.assign(2, 0, Core::LinAlg::View, *laig);

  // ---------Addressing contribution to block (4,4)
  mat.assign(2, 2, Core::LinAlg::View, aii);

  /*--------------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> mmm = fluid_field()->shape_derivatives();
  if (mmm != Teuchos::null)
  {
    // extract submatrices
    Core::LinAlg::SparseMatrix& fmii = mmm->matrix(0, 0);
    Core::LinAlg::SparseMatrix& fmgi = mmm->matrix(1, 0);

    // reuse transform objects to add shape derivative matrices to structural blocks

    // ---------Addressing contribution to block (2,2)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> fmgg =
        matrix_multiply(mmm->matrix(1, 1), false, *mortarp, false, false, false, true);
    fmgg = matrix_multiply(*mortarp, true, *fmgg, false, false, false, true);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> lfmgg =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(fmgg->row_map(), 81, false));
    lfmgg->add(*fmgg, false, 1.0, 0.0);
    lfmgg->complete(s->domain_map(), fmgg->range_map());

    s->add(*lfmgg, false, scale * (1. - stiparam) / (1. - ftiparam), 1.0);

    // ---------Addressing contribution to block (3,2)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> fmig =
        matrix_multiply(mmm->matrix(0, 1), false, *mortarp, false, false, false, true);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> lfmig =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(fmig->row_map(), 81, false));

    lfmig->add(*fmig, false, 1.0, 0.0);
    lfmig->complete(s->domain_map(), fmig->range_map());

    if (stcalgo != Inpar::Solid::stc_none)
    {
      lfmig = Core::LinAlg::matrix_multiply(*lfmig, false, *stcmat, false, false, false, true);
    }

    mat.matrix(1, 0).add(*lfmig, false, 1.0, 1.0);

    // We cannot copy the pressure value. It is not used anyway. So no exact
    // match here.
    (*fmiitransform_)(mmm->full_row_map(), mmm->full_col_map(), fmii, 1.,
        Coupling::Adapter::CouplingMasterConverter(coupfa), mat.matrix(1, 2), false);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> lfmgi =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(fmgi.row_map(), 81, false));
    (*fmiitransform_)(mmm->full_row_map(), mmm->full_col_map(), fmgi, 1.0,
        Coupling::Adapter::CouplingMasterConverter(coupfa), *lfmgi, false);

    // ---------Addressing contribution to block (2,4)
    lfmgi->complete(aii.domain_map(), mortarp->range_map());
    Teuchos::RCP<Core::LinAlg::SparseMatrix> llfmgi =
        matrix_multiply(*mortarp, true, *lfmgi, false, false, false, true);
    lfmgi = Teuchos::rcp(new Core::LinAlg::SparseMatrix(s->row_map(), 81, false));

    lfmgi->add(*llfmgi, false, scale, 0.0);
    lfmgi->complete(aii.domain_map(), s->range_map());

    if (stcalgo == Inpar::Solid::stc_currsym)
      lfmgi = Core::LinAlg::matrix_multiply(*stcmat, true, *lfmgi, false, true, true, false);
    lfmgi->scale((1. - stiparam) / (1. - ftiparam));
    mat.assign(0, 2, Core::LinAlg::View, *lfmgi);
  }

  s->complete();

  if (stcalgo != Inpar::Solid::stc_none)
  {
    s = Core::LinAlg::matrix_multiply(*s, false, *stcmat, false, true, true, true);

    if (stcalgo == Inpar::Solid::stc_currsym)
      s = Core::LinAlg::matrix_multiply(*stcmat, true, *s, false, true, true, false);
  }
  else
  {
    s->un_complete();
  }
  // finally assign structure matrix to block (0,0)
  mat.assign(0, 0, Core::LinAlg::View, *s);

  // done. make sure all blocks are filled.
  mat.complete();

  // Finally, take care of Dirichlet boundary conditions
  mat.apply_dirichlet(*(dbcmaps_->cond_map()), true);
  //
  // ---------------------------------------------------------------------------
  // END building the global system matrix
  // ---------------------------------------------------------------------------

  //  Teuchos::RCP<Epetra_CrsMatrix> matrix = mat.Matrix(0,0).EpetraMatrix();
  //  Core::LinAlg::PrintMatrixInMatlabFormat("mat.dat",*matrix,true);

  //  Core::LinAlg::PrintBlockMatrixInMatlabFormat("mat.dat",mat);
  //  std::cout<<"\nWROTE MATRIX!!";

  // ---------------------------------------------------------------------------
  // NOX related stuff needed for recovery of Lagrange multiplier
  // ---------------------------------------------------------------------------
  // store parts of fluid matrix to know them in the next iteration as previous
  // iteration matrices
  fgiprev_ = fgicur_;
  fggprev_ = fggcur_;
  fgicur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(f->matrix(1, 0)));
  fggcur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(f->matrix(1, 1)));

  // store parts of fluid shape derivative matrix to know them in the next
  // iteration as previous iteration matrices
  fmgiprev_ = fmgicur_;
  fmggprev_ = fmggcur_;
  if (mmm != Teuchos::null)
  {
    fmgicur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(mmm->matrix(1, 0)));
    fmggcur_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(mmm->matrix(1, 1)));
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::scale_system(
    Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)Core::UTILS::integral_value<int>(fsimono, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    // do scaling of structure rows
    Teuchos::RCP<Epetra_CrsMatrix> A = mat.matrix(0, 0).epetra_matrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if (A->LeftScale(*srowsum_) or A->RightScale(*scolsum_) or
        mat.matrix(0, 1).epetra_matrix()->LeftScale(*srowsum_) or
        mat.matrix(0, 2).epetra_matrix()->LeftScale(*srowsum_) or
        mat.matrix(1, 0).epetra_matrix()->RightScale(*scolsum_) or
        mat.matrix(2, 0).epetra_matrix()->RightScale(*scolsum_))
      FOUR_C_THROW("structure scaling failed");

    // do scaling of ale rows
    A = mat.matrix(2, 2).epetra_matrix();
    arowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    acolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*arowsum_);
    A->InvColSums(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_) or
        mat.matrix(2, 0).epetra_matrix()->LeftScale(*arowsum_) or
        mat.matrix(2, 1).epetra_matrix()->LeftScale(*arowsum_) or
        mat.matrix(0, 2).epetra_matrix()->RightScale(*acolsum_) or
        mat.matrix(1, 2).epetra_matrix()->RightScale(*acolsum_))
      FOUR_C_THROW("ale scaling failed");

    // do scaling of structure and ale rhs vectors
    Teuchos::RCP<Epetra_Vector> sx = extractor().extract_vector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = extractor().extract_vector(b, 2);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ax->Multiply(1.0, *arowsum_, *ax, 0.0)) FOUR_C_THROW("ale scaling failed");

    extractor().insert_vector(*sx, 0, b);
    extractor().insert_vector(*ax, 2, b);
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::unscale_solution(
    Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)Core::UTILS::integral_value<int>(fsimono, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> sy = extractor().extract_vector(x, 0);
    Teuchos::RCP<Epetra_Vector> ay = extractor().extract_vector(x, 2);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ay->Multiply(1.0, *acolsum_, *ay, 0.0)) FOUR_C_THROW("ale scaling failed");

    // get info about STC feature and unscale solution if necessary
    Inpar::Solid::StcScale stcalgo = structure_field()->get_stc_algo();
    if (stcalgo != Inpar::Solid::stc_none)
    {
      structure_field()->get_stc_mat()->multiply(false, *sy, *sy);
    }

    extractor().insert_vector(*sy, 0, x);
    extractor().insert_vector(*ay, 2, x);

    Teuchos::RCP<Epetra_Vector> sx = extractor().extract_vector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = extractor().extract_vector(b, 2);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ax->ReciprocalMultiply(1.0, *arowsum_, *ax, 0.0)) FOUR_C_THROW("ale scaling failed");

    // get info about STC feature
    if (stcalgo != Inpar::Solid::stc_none)
    {
      structure_field()->get_stc_mat()->multiply(false, *sx, *sx);
    }

    extractor().insert_vector(*sx, 0, b);
    extractor().insert_vector(*ax, 2, b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.matrix(0, 0).epetra_matrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if (A->LeftScale(*srowsum_) or A->RightScale(*scolsum_) or
        mat.matrix(0, 1).epetra_matrix()->LeftScale(*srowsum_) or
        mat.matrix(0, 2).epetra_matrix()->LeftScale(*srowsum_) or
        mat.matrix(1, 0).epetra_matrix()->RightScale(*scolsum_) or
        mat.matrix(2, 0).epetra_matrix()->RightScale(*scolsum_))
      FOUR_C_THROW("structure scaling failed");

    A = mat.matrix(2, 2).epetra_matrix();
    arowsum_->Reciprocal(*arowsum_);
    acolsum_->Reciprocal(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_) or
        mat.matrix(2, 0).epetra_matrix()->LeftScale(*arowsum_) or
        mat.matrix(2, 1).epetra_matrix()->LeftScale(*arowsum_) or
        mat.matrix(0, 2).epetra_matrix()->RightScale(*acolsum_) or
        mat.matrix(1, 2).epetra_matrix()->RightScale(*acolsum_))
      FOUR_C_THROW("ale scaling failed");
  }

  // very simple hack just to see the linear solution

  Epetra_Vector r(b.Map());
  mat.Apply(x, r);
  r.Update(1., b, 1.);

  Teuchos::RCP<Epetra_Vector> sr = extractor().extract_vector(r, 0);
  Teuchos::RCP<Epetra_Vector> fr = extractor().extract_vector(r, 1);
  Teuchos::RCP<Epetra_Vector> ar = extractor().extract_vector(r, 2);

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

  if (structure_field()->get_stc_algo() != Inpar::Solid::stc_none)
    structure_field()->system_matrix()->reset();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Combo> FSI::SlidingMonolithicFluidSplit::create_status_test(
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
  interface.push_back(structure_field()->interface()->fsi_cond_map());
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

  // ---------------------------------------------------------------------------
  // setup tests for fluid velocity field
  // ---------------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidvel;
  fluidvel.push_back(fluid_field()->inner_velocity_row_map());
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

  // ---------------------------------------------------------------------------
  // setup tests for fluid pressure field
  // ---------------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidpress;
  fluidpress.push_back(fluid_field()->pressure_row_map());
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

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
    Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::extract_field_vectors");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (ddgpred_ == Teuchos::null)
    FOUR_C_THROW("Vector 'ddgpred_' has not been initialized properly.");
#endif

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<const Core::LinAlg::SparseMatrix> mortarp = coupsfm_->get_mortar_matrix_p();

  // ---------------------------------------------------------------------------
  // process structure unknowns
  // ---------------------------------------------------------------------------
  // extract structure solution increment from NOX increment
  sx = extractor().extract_vector(x, 0);

  // ---------------------------------------------------------------------------
  // process ale unknowns
  // ---------------------------------------------------------------------------
  // extract inner ALE solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> aox = extractor().extract_vector(x, 2);

  // convert structure solution increment to ALE solution increment at the interface
  Teuchos::RCP<Epetra_Vector> scx = structure_field()->interface()->extract_fsi_cond_vector(sx);
  scx->Update(1.0, *ddgpred_, 1.0);
  Teuchos::RCP<Epetra_Vector> acx =
      Core::LinAlg::create_vector(*fluid_field()->interface()->fsi_cond_map());
  mortarp->Apply(*scx, *acx);
  acx = fluid_to_ale_interface(acx);

  // put inner and interface ALE solution increments together
  Teuchos::RCP<Epetra_Vector> a = fsi_ale_field()->fsi_interface()->insert_other_vector(aox);
  ale_field()->interface()->insert_fsi_cond_vector(acx, a);
  ale_field()->update_slave_dof(a);
  ax = a;

  // ---------------------------------------------------------------------------
  // process fluid unknowns
  // ---------------------------------------------------------------------------
  // extract inner fluid solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> fox = extractor().extract_vector(x, 1);
  fox = fsi_fluid_field()->fsi_interface()->extract_other_vector(fox);

  // convert ALE solution increment to fluid solution increment at the interface
  Teuchos::RCP<Epetra_Vector> fcx = ale_to_fluid_interface(acx);
  fluid_field()->displacement_to_velocity(fcx);

  // put inner and interface fluid solution increments together
  Teuchos::RCP<Epetra_Vector> f = fsi_fluid_field()->fsi_interface()->insert_other_vector(fox);
  fluid_field()->interface()->insert_fsi_cond_vector(fcx, f);
  fluid_field()->update_slave_dof(f);
  fx = f;

  // ---------------------------------------------------------------------------

  // Store field vectors to know them later on as previous quantities
  // interface structure displacement increment
  if (disgprev_ != Teuchos::null)
    ddginc_->Update(1.0, *scx, -1.0, *disgprev_, 0.0);  // compute current iteration increment
  else
    ddginc_ = Teuchos::rcp(new Epetra_Vector(*scx));  // first iteration increment

  disgprev_ = scx;  // store current step increment
  // ------------------------------------

  // inner ale displacement increment
  if (aleiprev_ != Teuchos::null)
    ddialeinc_->Update(1.0, *aox, -1.0, *aleiprev_, 0.0);  // compute current iteration increment
  else
    ddialeinc_ = Teuchos::rcp(new Epetra_Vector(*aox));  // first iteration increment

  aleiprev_ = aox;  // store current step increment
  // ------------------------------------

  // inner fluid solution increment
  if (veliprev_ != Teuchos::null)  // compute current iteration increment
    duiinc_->Update(1.0, *fox, -1.0, *veliprev_, 0.0);
  else
    // first iteration increment
    duiinc_ = Teuchos::rcp(new Epetra_Vector(*fox));
  // store current step increment
  veliprev_ = fox;
  // ------------------------------------
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::update()
{
  lambdaold_->Update(1.0, *lambda_, 0.0);

  // update history variabels for sliding ale
  if (aleproj_ != Inpar::FSI::ALEprojection_none)
  {
    iprojdisp_ = Teuchos::rcp(new Epetra_Vector(*coupsfm_->slave_dof_map(), true));
    Teuchos::RCP<Epetra_Vector> idispale = ale_to_fluid_interface(
        ale_field()->interface()->extract_fsi_cond_vector(ale_field()->dispnp()));

    slideale_->remeshing(*structure_field(), fluid_field()->discretization(), idispale, iprojdisp_,
        *coupsfm_, get_comm());

    iprojdispinc_->Update(-1.0, *iprojdisp_, 1.0, *idispale, 0.0);

    slideale_->evaluate_mortar(
        structure_field()->extract_interface_dispnp(), iprojdisp_, *coupsfm_);
    slideale_->evaluate_fluid_mortar(idispale, iprojdisp_);

    Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*iprojdisp_));
    temp->ReplaceMap(idispale->Map());
    Teuchos::RCP<Epetra_Vector> acx = fluid_to_ale_interface(temp);
    ale_field()->apply_interface_displacements(acx);
    fluid_field()->apply_mesh_displacement(ale_to_fluid(ale_field()->dispnp()));

    Teuchos::RCP<Epetra_Vector> unew =
        slideale_->interpolate_fluid(fluid_field()->extract_interface_velnp());
    fluid_field()->apply_interface_velocities(unew);
  }

  // call update()-routine in base class to handle the single fields
  FSI::BlockMonolithic::update();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::output()
{
  structure_field()->output();
  fluid_field()->output();

  if (aleproj_ != Inpar::FSI::ALEprojection_none)
  {
    int uprestart = timeparams_.get<int>("RESTARTEVRY");
    if (uprestart != 0 && fluid_field()->step() % uprestart == 0)
    {
      fluid_field()->disc_writer()->write_vector("slideALE", iprojdisp_);
      fluid_field()->disc_writer()->write_vector("slideALEincr", iprojdispinc_);
      slideale_->output_restart(*fluid_field()->disc_writer());
    }
  }

  // output Lagrange multiplier
  output_lambda();

  ale_field()->output();

  if (structure_field()->get_constraint_manager()->have_monitor())
  {
    structure_field()->get_constraint_manager()->compute_monitor_values(
        structure_field()->dispnp());
    if (comm_.MyPID() == 0) structure_field()->get_constraint_manager()->print_monitor_values();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::output_lambda()
{
  /* 'lambda_' is only defined on the interface. So, insert 'lambda_' into
   * 'lambdafull' that is defined on the entire fluid field. Then, write
   * output or restart data.
   */
  Teuchos::RCP<Epetra_Vector> lambdafull =
      fluid_field()->interface()->insert_fsi_cond_vector(lambda_);
  const int uprestart = timeparams_.get<int>("RESTARTEVRY");
  const int upres = timeparams_.get<int>("RESULTSEVRY");
  if ((uprestart != 0 && fluid_field()->step() % uprestart == 0) ||
      fluid_field()->step() % upres == 0)
    fluid_field()->disc_writer()->write_vector("fsilambda", lambdafull);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::read_restart(int step)
{
  structure_field()->read_restart(step);
  fluid_field()->read_restart(step);

  auto input_control_file = Global::Problem::instance()->input_control_file();

  // read Lagrange multiplier
  {
    Teuchos::RCP<Epetra_Vector> lambdafull =
        Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(), true));
    Core::IO::DiscretizationReader reader =
        Core::IO::DiscretizationReader(fluid_field()->discretization(), input_control_file, step);
    reader.read_vector(lambdafull, "fsilambda");
    lambdaold_ = fluid_field()->interface()->extract_fsi_cond_vector(lambdafull);
    // Note: the above is normally enough. However, we can use the restart in order to periodically
    // repeat the fsi simulation (see AC-FS3I)
    lambda_ = fluid_field()->interface()->extract_fsi_cond_vector(lambdafull);
  }

  setup_system();

  if (aleproj_ != Inpar::FSI::ALEprojection_none)
  {
    Core::IO::DiscretizationReader reader =
        Core::IO::DiscretizationReader(fluid_field()->discretization(), input_control_file, step);
    reader.read_vector(iprojdisp_, "slideALE");
    reader.read_vector(iprojdispinc_, "slideALEincr");
    slideale_->read_restart(reader);
  }
  ale_field()->read_restart(step);

  set_time_step(fluid_field()->time(), fluid_field()->step());

  if (aleproj_ != Inpar::FSI::ALEprojection_none)
    slideale_->evaluate_mortar(structure_field()->extract_interface_dispn(), iprojdisp_, *coupsfm_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::prepare_time_step()
{
  precondreusecount_ = 0;

  increment_time_and_step();
  print_header();

  prepare_time_step_preconditioner();

  if (structure_field()->get_stc_algo() != Inpar::Solid::stc_none)
    structure_field()->system_matrix()->reset();

  prepare_time_step_fields();

  // Note: it's important to first prepare the single fields and than the fsi problem
  prepare_time_step_fsi();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::recover_lagrange_multiplier()
{
  // get time integration parameter of fluid time integrator
  // to enable consistent time integration among the fields
  const double ftiparam = fluid_field()->tim_int_param();

  // some scaling factors for fluid
  const double timescale = fluid_field()->time_scaling();
  const double scale = fluid_field()->residual_scaling();

  // get the Mortar projection matrix P = D^{-1} * M
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortarp = coupsfm_->get_mortar_matrix_p();

  // get the inverted Mortar matrix D^{-1}
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortardinv = coupsfm_->get_mortar_matrix_dinv();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (mortarp == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix P.");
  if (mortardinv == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix D^{-1}.");
#endif

  // get fluid shape derivative matrix
  const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> mmm = fluid_field()->shape_derivatives();

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
   * (2)  - 1. / (1.-ftiparam) * D^{-T} * tmpvec
   *
   * with tmpvec =
   *
   * (3)    r_{\Gamma}^{F,n+1}
   *
   * (4)  + 1 / tau * F_{\Gamma\Gamma} * P * \Delta d_{\Gamma}^{S,n+1}
   *
   * (5)  + F_{\Gamma\Gamma}^{G} * P * \Delta d_{\Gamma}^{S,n+1}
   *
   * (6)  + F_{\Gamma I} * \Delta u_{I}^{F,n+1}
   *
   * (7)  + F_{\Gamma I}^{G} * \Delta d_{I}^{G,n+1}
   *
   * (8)  - dt / tau * F_{\Gamma\Gamma} * u_{\Gamma}^n]
   *
   * Remark on term (8):
   * Term (8) has to be considered only in the first Newton iteration.
   * Hence, it will usually not be computed since in general we need more
   * than one nonlinear iteration until convergence.
   *
   * Remarks on all terms:
   * +  Division by (1.0 - ftiparam) will be done in the end
   *    since this is common to all terms
   * +  tau: time scaling factor for interface time integration (tau =
   * 1/fluid_field()->TimeScaling())
   * +  neglecting terms (4)-(8) should not alter the results significantly
   *    since at the end of the time step the solution increments tend to zero.
   *
   * Terms arising from field specific predictors have to be considered only in
   * the first Newton iteration. Since we usually need more than on iteration,
   * these terms are not implemented, yet.
   *
   *                                                 Matthias Mayr (10/2012)
   */

  // ---------Addressing term (1)
  lambda_->Update(ftiparam, *lambdaold_, 0.0);
  // ---------End of term (1)

  // ---------Addressing term (3)
  Teuchos::RCP<Epetra_Vector> fluidresidual =
      fluid_field()->interface()->extract_fsi_cond_vector(fluid_field()->rhs());
  fluidresidual->Scale(-1.0);
  tmpvec = Teuchos::rcp(new Epetra_Vector(*fluidresidual));
  // ---------End of term (3)

  // ---------Addressing term (4)
  auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->range_map(), true));
  mortarp->Apply(*ddginc_, *auxvec);
  auxauxvec = Teuchos::rcp(new Epetra_Vector(fggprev_->range_map(), true));
  fggprev_->Apply(*auxvec, *auxauxvec);
  tmpvec->Update(timescale, *auxauxvec, 1.0);
  // ---------End of term (4)

  // ---------Addressing term (5)
  if (fmggprev_ != Teuchos::null)
  {
    auxvec = Teuchos::rcp(new Epetra_Vector(mortarp->range_map(), true));
    mortarp->Apply(*ddginc_, *auxvec);
    fmggprev_->Apply(*auxvec, *auxauxvec);
    tmpvec->Update(1.0, *auxauxvec, 1.0);
  }
  // ---------End of term (5)

  // ---------Addressing term (6)
  auxvec = Teuchos::rcp(new Epetra_Vector(fgiprev_->range_map(), true));
  fgiprev_->Apply(*duiinc_, *auxvec);
  tmpvec->Update(1.0, *auxvec, 1.0);
  // ---------End of term (6)

  // ---------Addressing term (7)
  if (fmgiprev_ != Teuchos::null)
  {
    /* For matrix-vector-product, the DomainMap() of the matrix and the Map()
     * of the vector have to match. DomaintMap() contains inner velocity DOFs
     * and all pressure DOFs. The inner ale displacement increment is converted
     * to the fluid map using AleToFluid(). This results in a map that contains
     * all velocity but no pressure DOFs.
     *
     * We have to circumvent some trouble with Epetra_BlockMaps since we cannot
     * split an Epetra_BlockMap into inner and interface DOFs.
     *
     * We create a map extractor 'velothermap' in order to extract the inner
     * velocity DOFs after calling AleToFluid(). Afterwards, a second map
     * extractor 'velotherpressuremapext' is used to append pressure DOFs filled
     * with zeros.
     *
     * Finally, maps match and matrix-vector-multiplication can be done.
     */

    // extract inner velocity DOFs after calling AleToFluid()
    Teuchos::RCP<Epetra_Map> velothermap = Core::LinAlg::split_map(
        *fluid_field()->velocity_row_map(), *interface_fluid_ale_coupling().master_dof_map());
    Core::LinAlg::MapExtractor velothermapext =
        Core::LinAlg::MapExtractor(*fluid_field()->velocity_row_map(), velothermap, false);
    auxvec = Teuchos::rcp(new Epetra_Vector(*velothermap, true));
    velothermapext.extract_other_vector(
        ale_to_fluid(fsi_ale_field()->fsi_interface()->insert_other_vector(ddialeinc_)), auxvec);

    // add pressure DOFs
    Core::LinAlg::MapExtractor velotherpressuremapext =
        Core::LinAlg::MapExtractor(fmgiprev_->domain_map(), velothermap);
    auxauxvec = Teuchos::rcp(new Epetra_Vector(fmgiprev_->domain_map(), true));
    velotherpressuremapext.insert_cond_vector(auxvec, auxauxvec);

    // prepare vector to store result of matrix-vector-product
    auxvec = Teuchos::rcp(new Epetra_Vector(fmgiprev_->range_map(), true));

    // Now, do the actual matrix-vector-product
    fmgiprev_->Apply(*auxauxvec, *auxvec);
    tmpvec->Update(1.0, *auxvec, 1.0);
  }
  // ---------End of term (7)

  // ---------Addressing term (8)
  if (firstcall_)
  {
    auxvec = Teuchos::rcp(new Epetra_Vector(fggprev_->range_map(), true));
    fggprev_->Apply(*fluid_field()->extract_interface_veln(), *auxvec);
    tmpvec->Update(dt() * timescale, *auxvec, 1.0);
  }
  // ---------End of term (8)

  // ---------Addressing term (2)
  auxvec = Teuchos::rcp(new Epetra_Vector(mortardinv->domain_map(), true));
  mortardinv->multiply(true, *tmpvec, *auxvec);
  lambda_->Update(scale, *auxvec, 1.0);  // scale with residual_scaling() to get [N/m^2]
  // ---------End of term (2)

  // Finally, divide by (1.0-ftiparam) which is common to all terms
  lambda_->Scale(-1.0 / (1.0 - ftiparam));

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
void FSI::SlidingMonolithicFluidSplit::calculate_interface_energy_increment()
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = structure_field()->tim_int_param();
  const double ftiparam = fluid_field()->tim_int_param();

  // get the Mortar matrix M
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortarm = coupsfm_->get_mortar_matrix_m();

  // interface traction weighted by time integration factors
  Teuchos::RCP<Epetra_Vector> tractionfluid = Teuchos::rcp(new Epetra_Vector(lambda_->Map(), true));
  Teuchos::RCP<Epetra_Vector> tractionstructure =
      Teuchos::rcp(new Epetra_Vector(*structure_field()->interface()->fsi_cond_map(), true));
  tractionfluid->Update(stiparam - ftiparam, *lambdaold_, ftiparam - stiparam, *lambda_, 0.0);
  mortarm->multiply(true, *tractionfluid, *tractionstructure);

  // displacement increment of this time step
  Teuchos::RCP<Epetra_Vector> deltad =
      Teuchos::rcp(new Epetra_Vector(*structure_field()->dof_row_map(), true));
  deltad->Update(1.0, *structure_field()->dispnp(), -1.0, *structure_field()->dispn(), 0.0);

  // calculate the energy increment
  double energy = 0.0;
  tractionstructure->Dot(*structure_field()->interface()->extract_fsi_cond_vector(deltad), &energy);

  energysum_ += energy;

  write_interface_energy_file(energy, energysum_);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::check_kinematic_constraint()
{
  // some scaling factors for fluid
  const double timescale = fluid_field()->time_scaling();

  // get the Mortar matrices D and M
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortard = coupsfm_->get_mortar_matrix_d();
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortarm = coupsfm_->get_mortar_matrix_m();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (mortarm == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix M.");
  if (mortard == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix D.");
#endif

  // get interface displacements and velocities
  const Teuchos::RCP<Epetra_Vector> disnp = structure_field()->extract_interface_dispnp();
  const Teuchos::RCP<Epetra_Vector> disn = structure_field()->extract_interface_dispn();
  const Teuchos::RCP<Epetra_Vector> velnp = fluid_field()->extract_interface_velnp();
  const Teuchos::RCP<Epetra_Vector> veln = fluid_field()->extract_interface_veln();

  // prepare vectors for projected interface quantities
  Teuchos::RCP<Epetra_Vector> disnpproj =
      Teuchos::rcp(new Epetra_Vector(mortarm->range_map(), true));
  Teuchos::RCP<Epetra_Vector> disnproj =
      Teuchos::rcp(new Epetra_Vector(mortarm->range_map(), true));
  Teuchos::RCP<Epetra_Vector> velnpproj =
      Teuchos::rcp(new Epetra_Vector(mortard->range_map(), true));
  Teuchos::RCP<Epetra_Vector> velnproj =
      Teuchos::rcp(new Epetra_Vector(mortard->range_map(), true));

  // projection of interface displacements
  mortarm->Apply(*disnp, *disnpproj);
  mortarm->Apply(*disn, *disnproj);

  // projection of interface velocities
  mortard->Apply(*velnp, *velnpproj);
  mortard->Apply(*veln, *velnproj);

  // calculate violation of kinematic interface constraint
  Teuchos::RCP<Epetra_Vector> violation = Teuchos::rcp(new Epetra_Vector(*disnpproj));
  violation->Update(-1.0, *disnproj, 1.0);
  violation->Update(-1.0 / timescale, *velnpproj, 1.0 / timescale, *velnproj, 1.0);
  violation->Update(-dt(), *velnproj, 1.0);

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
void FSI::SlidingMonolithicFluidSplit::check_dynamic_equilibrium()
{
  // get the Mortar matrices D and M
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortard = coupsfm_->get_mortar_matrix_d();
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortarm = coupsfm_->get_mortar_matrix_m();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (mortarm == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix M.");
  if (mortard == Teuchos::null) FOUR_C_THROW("Expected Teuchos::rcp to mortar matrix D.");
#endif

  // auxiliary vectors
  Teuchos::RCP<Epetra_Vector> tractionmaster =
      Teuchos::rcp(new Epetra_Vector(mortarm->domain_map(), true));
  Teuchos::RCP<Epetra_Vector> tractionslave =
      Teuchos::rcp(new Epetra_Vector(mortard->domain_map(), true));

  // calculate tractions on master and slave side
  mortarm->multiply(true, *lambda_, *tractionmaster);
  mortard->multiply(true, *lambda_, *tractionslave);

  // calculate violation of dynamic equilibrium
  Teuchos::RCP<Epetra_Vector> violation = Teuchos::rcp(new Epetra_Vector(*tractionmaster));
  violation->Update(-1.0, *tractionslave, 1.0);

  // calculate some norms
  double violationl2 = 0.0;
  double violationinf = 0.0;
  violation->Norm2(&violationl2);
  violation->NormInf(&violationinf);

  // scale L2-Norm with sqrt of length of interface vector
  violationl2 /= sqrt(fluid_field()->interface()->fsi_cond_map()->NumGlobalElements());

  // output to screen
  std::ios_base::fmtflags flags = utils()->out().flags();

  utils()->out() << std::scientific << "\nViolation of dynamic interface equilibrium:\n"
                 << "L_2-norm: " << violationl2 << "        L_inf-norm: " << violationinf << "\n";
  utils()->out().flags(flags);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::combine_field_vectors(Epetra_Vector& v,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv,
    Teuchos::RCP<const Epetra_Vector> av, const bool slave_vectors_contain_interface_dofs)
{
  if (slave_vectors_contain_interface_dofs)
  {
    // extract inner DOFs from slave vectors
    Teuchos::RCP<Epetra_Vector> fov = fsi_fluid_field()->fsi_interface()->extract_other_vector(fv);
    fov = fsi_fluid_field()->fsi_interface()->insert_other_vector(fov);
    Teuchos::RCP<Epetra_Vector> aov = fsi_ale_field()->fsi_interface()->extract_other_vector(av);

    // put them together
    FSI::Monolithic::combine_field_vectors(v, sv, fov, aov);
  }
  else
    FSI::Monolithic::combine_field_vectors(v, sv, fv, av);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double FSI::SlidingMonolithicFluidSplit::select_dt_error_based() const
{
  // get time step size suggestions
  const double dtstr = get_ada_str_dt();           // based on all structure DOFs
  const double dtstrfsi = get_ada_str_fsi_dt();    // based on structure FSI DOFs
  const double dtflinner = get_ada_fl_inner_dt();  // based on inner fluid DOFs

  double dt = SlidingMonolithicFluidSplit::dt();

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

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool FSI::SlidingMonolithicFluidSplit::set_accepted() const
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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::SlidingMonolithicFluidSplit::create_system_matrix()
{
  FSI::BlockMonolithic::create_system_matrix(systemmatrix_, false);
}

FOUR_C_NAMESPACE_CLOSE
