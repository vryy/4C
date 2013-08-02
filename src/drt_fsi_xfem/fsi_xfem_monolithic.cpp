/*!----------------------------------------------------------------------
\file fsi_monolithic_xfem.cpp
\brief Control routine for monolithic XFSI using XFEM

<pre>
Maintainer:  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>
*----------------------------------------------------------------------*/

#include "fsi_xfem_monolithic.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_fsi/fsi_matrixtransform.H"
#include "../drt_adapter/ad_fld_xfluid_fsi.H"

//#include "fsi_overlapprec_fsiamg.H"
#include "../drt_fsi/fsi_debugwriter.H"
#include "../drt_fsi/fsi_statustest.H"
#include "../drt_fsi/fsi_monolithic_linearsystem.H"

#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_solver.H"

#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_structure/stru_aux.H"
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_ale/ale.H"
#include "../drt_inpar/inpar_xfem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_io/io_pstream.H"

#include "../drt_inpar/inpar_ale.H"

#include <Teuchos_TimeMonitor.hpp>



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicXFEM::MonolithicXFEM(const Epetra_Comm& comm,
                                    const Teuchos::ParameterList& timeparams)
  : AlgorithmXFEM(comm, timeparams),
    fsidyn_(DRT::Problem::Instance()->FSIDynamicParams())
{

  blockrowdofmap_ = Teuchos::rcp(new LINALG::MultiMapExtractor);

  firstcall_ = true;

  itermin_ = 0;
  iter_ = 0;
  iter_outer_ = 0;
  itermax_outer_ = 5;

  tolrhs_ = 1e-12;

//  // Throw an error if there are DBCs on structural interface DOFs.
//  std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmaps;
//  intersectionmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
//  intersectionmaps.push_back(StructureField()->Interface()->FSICondMap());
//  Teuchos::RCP<Epetra_Map> intersectionmap = LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);
//
//  if (intersectionmap->NumGlobalElements() != 0)
//  {
//    // remove interface DOFs from structural DBC map
//    StructureField()->RemoveDirichCond(intersectionmap);
//
//    // give a warning to the user that Dirichlet boundary conditions might not be correct
//
//   if (comm.MyPID() == 0)
//    {
//      IO::cout << "  +---------------------------------------------------------------------------------------------+" << IO::endl;
//      IO::cout << "  |                                        PLEASE NOTE:                                         |" << IO::endl;
//      IO::cout << "  +---------------------------------------------------------------------------------------------+" << IO::endl;
//      IO::cout << "  | You run a monolithic structure split scheme. Hence, there are no structural interface DOFs. |" << IO::endl;
//      IO::cout << "  | Structure Dirichlet boundary conditions on the interface will be neglected.                 |" << IO::endl;
//      IO::cout << "  | Check whether you have prescribed appropriate DBCs on structural interface DOFs.            |" << IO::endl;
//      IO::cout << "  +---------------------------------------------------------------------------------------------+" << IO::endl;
//    }
//  }
//
//#ifdef DEBUG
//  // check if removing Dirichlet conditions was successful
//  intersectionmaps.resize(0);
//  intersectionmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
//  intersectionmaps.push_back(StructureField()->Interface()->FSICondMap());
//  intersectionmap = LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);
//  if (intersectionmap->NumGlobalElements() != 0)
//    dserror("Could not remove structural interface Dirichlet conditions from structure DBC map.");
//#endif
//
//  sggtransform_  = Teuchos::rcp(new UTILS::MatrixRowColTransform);
//  sgitransform_  = Teuchos::rcp(new UTILS::MatrixRowTransform);
//  sigtransform_  = Teuchos::rcp(new UTILS::MatrixColTransform);
//  aigtransform_  = Teuchos::rcp(new UTILS::MatrixColTransform);
//  fmiitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
//  fmgitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
//  fsaigtransform_= Teuchos::rcp(new UTILS::MatrixColTransform);
//  fsmgitransform_= Teuchos::rcp(new UTILS::MatrixColTransform);
//
//  // const Teuchos::ParameterList& xdyn = DRT::Problem::Instance()->XFEMGeneralParams();
//  const Teuchos::ParameterList& xfluiddyn  = DRT::Problem::Instance()->XFluidDynamicParams();
//  Teuchos::RCP<Teuchos::ParameterList> fluidtimeparams = Teuchos::rcp(new Teuchos::ParameterList());
//
//  monolithic_approach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>
//                         (xfluiddyn.sublist("GENERAL"),"MONOLITHIC_XFFSI_APPROACH");
//
//  currentstep_ = 0;
//  relaxing_ale_ = (bool)DRT::INPUT::IntegralValue<int>(xfluiddyn.sublist("GENERAL"),"RELAXING_ALE");
//  relaxing_ale_every_ = xfluiddyn.sublist("GENERAL").get<int>("RELAXING_ALE_EVERY");
//
//  const Teuchos::ParameterList& adyn     = DRT::Problem::Instance()->AleDynamicParams();
//  int aletype = DRT::INPUT::IntegralValue<int>(adyn,"ALE_TYPE");
//
//  if ((aletype!=INPAR::ALE::incr_lin) and (monolithic_approach_!=INPAR::XFEM::XFFSI_Full_Newton))
//    dserror("Relaxing Ale Aprooach is just posiible with Ale-incr-lin!");
//
//  // Recovering of Lagrange multiplier happens on structure field
//  lambda_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap()));
//  ddiinc_ = Teuchos::null;
//  solipre_= Teuchos::null;
//  ddginc_ = Teuchos::null;
//  solgpre_= Teuchos::null;
//  ddgpred_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap(),true));
//
////fgpre_ = Teuchos::null;
////sgipre_ = Teuchos::null;
////sggpre_ = Teuchos::null;

  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  // enable debugging
  if (DRT::INPUT::IntegralValue<int>(fsidyn,"DEBUGOUTPUT")==1)
  {
    sdbg_ = Teuchos::rcp(new UTILS::DebugWriter(StructureField()->Discretization()));
    fdbg_ = Teuchos::rcp(new UTILS::DebugWriter(FluidField()->Discretization()));
  }

  std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
  s.append(".iteration");
  log_ = Teuchos::rcp(new std::ofstream(s.c_str()));
  itermax_ = fsidyn.get<int>("ITEMAX");
  normtypeinc_
    = DRT::INPUT::IntegralValue<INPAR::FSI::ConvNorm>(fsidyn,"NORM_INC");
  normtypefres_
    = DRT::INPUT::IntegralValue<INPAR::FSI::ConvNorm>(fsidyn,"NORM_RESF");
  combincfres_
    = DRT::INPUT::IntegralValue<INPAR::FSI::BinaryOp>(fsidyn,"NORMCOMBI_RESFINC");
  tolinc_ =  fsidyn.get<double>("CONVTOL");
  tolfres_ = fsidyn.get<double>("CONVTOL");

  TOL_DIS_RES_L2_  =  fsidyn.get<double>("TOL_DIS_RES_L2");
  TOL_DIS_RES_INF_ =  fsidyn.get<double>("TOL_DIS_RES_INF");
  TOL_DIS_INC_L2_  =  fsidyn.get<double>("TOL_DIS_RES_L2");
  TOL_DIS_INC_INF_ =  fsidyn.get<double>("TOL_DIS_RES_INF");
  TOL_PRE_RES_L2_  =  fsidyn.get<double>("TOL_PRE_RES_L2");
  TOL_PRE_RES_INF_ =  fsidyn.get<double>("TOL_PRE_RES_INF");
  TOL_PRE_INC_L2_  =  fsidyn.get<double>("TOL_PRE_RES_L2");
  TOL_PRE_INC_INF_ =  fsidyn.get<double>("TOL_PRE_RES_INF");
  TOL_VEL_RES_L2_  =  fsidyn.get<double>("TOL_VEL_RES_L2");
  TOL_VEL_RES_INF_ =  fsidyn.get<double>("TOL_VEL_RES_INF");
  TOL_VEL_INC_L2_  =  fsidyn.get<double>("TOL_VEL_RES_L2");
  TOL_VEL_INC_INF_ =  fsidyn.get<double>("TOL_VEL_RES_INF");
  // set tolerances for nonlinear solver

//  //Structural predictor step, initially filled with zeros
//  ddgpred_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap(),true));



  //----------------------------------------------------
  // XFSI solver: create a linear solver

  merge_fsi_blockmatrix_ = true;

  // get iterative solver
  if (merge_fsi_blockmatrix_ == false)
    CreateLinearSolver();
  // get direct solver, e.g. UMFPACK
  else  // (merge_fsi_blockmatrix_ == true)
  {

    if (Comm().MyPID() == 0)
      std::cout << "Merged FSI block matrix is used!\n" << std::endl;

    Teuchos::RCP<Teuchos::ParameterList> solverparams = Teuchos::rcp(new Teuchos::ParameterList);
    solverparams->set("solver","umfpack");

    solver_ = Teuchos::rcp(new LINALG::Solver(
                        solverparams,
                        Comm(),
                        DRT::Problem::Instance()->ErrorFile()->Handle()
                        )
                  );
  }  // end BlockMatrixMerge

  return;
}


/*----------------------------------------------------------------------*
 | full monolithic dof row map                             schott 07/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FSI::MonolithicXFEM::DofRowMap() const
{
  return Extractor()->FullMap();
}  // DofRowMap()


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupSystem()
{
  IO::cout << "FSI::MonolithicXFEM::SetupSystem() still has to be implemented" << IO::endl;


  //Extract parameter list FSI_DYNAMIC
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  //Extract information about the linear block solver (FSIAMG / PreconditionedKrylov)
  linearsolverstrategy_ = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

//  //Dimensionality of problem
//  const int ndim = DRT::Problem::Instance()->NDim();


  /*----------------------------------------------------------------------
  Create a combined map for Structure/Fluid-DOFs all in one!
  ----------------------------------------------------------------------*/
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;

  //Append the structural DOF map
  vecSpaces.push_back(StructureField()->DofRowMap());

  //Append the background fluid DOF map
  vecSpaces.push_back(FluidField()->DofRowMap());

  //solid maps empty??
  if (vecSpaces[0]->NumGlobalElements() == 0)
    dserror("No solid equations. Panic.");

  //fluid maps empty??
  if (vecSpaces[1]->NumGlobalElements() == 0)
    dserror("No fluid equations. Panic.");

  //The vector is complete, now fill the system's global block row map
  //with the maps previously set together!
  SetDofRowMaps(vecSpaces);

  //TODO: check this
//  // Use normal matrix for fluid equations but build (splitted) mesh movement
//  // linearization (shapederivatives) if requested in the input file
//  FluidField().UseBlockMatrix(true);


  /*----------------------------------------------------------------------*/
  // create block system matrix

  switch(linearsolverstrategy_)
  {
  case INPAR::FSI::PreconditionedKrylov:
  {
    /*----------------------------------------------------------------------*/
    // initialise XFSI-systemmatrix_
    systemmatrix_
      = Teuchos::rcp(
          new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                *Extractor(),
                *Extractor(),
                81,
                false,
                true
                )
          );
    break;
  }
//  case INPAR::FSI::FSIAMG:
//    systemmatrix_ = Teuchos::rcp(new OverlappingBlockMatrixFSIAMG(
//                                   Extractor(),
//                                   *StructureField(),
//                                   FluidField(),
//                                   AleField(),
//                                   false,
//                                   DRT::INPUT::IntegralValue<int>(fsidyn,"SYMMETRICPRECOND"),
//                                   blocksmoother,
//                                   schuromega,
//                                   pcomega,
//                                   pciter,
//                                   spcomega,
//                                   spciter,
//                                   fpcomega,
//                                   fpciter,
//                                   apcomega,
//                                   apciter,
//                                   DRT::INPUT::IntegralValue<int>(fsidyn,"FSIAMGANALYZE"),
//                                   linearsolverstrategy_,
//                                   DRT::Problem::Instance()->ErrorFile()->Handle()));
//    break;
  default:
    dserror("Unsupported type of monolithic solver");
    break;
  }


}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupRHS(Epetra_Vector& f, bool firstcall)
{


  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::SetupRHS");

  // here we add the structural part, the fluid part and the coupling rhs into the global rhs vector
  SetupVector(f,
              StructureField()->RHS(),
              FluidField()->RHS(),
              FluidField()->ResidualScaling());


  if (firstcall)
  {

    std::cout << "What to do in SetupRHS for the firstcall ?!" << std::endl;
//     // get structure matrix
//    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocks = StructureField()->BlockSystemMatrix();
//    if (blocks==Teuchos::null)
//      dserror("expect structure block matrix");
//
//    // get fluid shape derivatives matrix
//    Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();
//
//    //get ale matrix
//    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocka = AleField().BlockSystemMatrix();
//    if (blocka==Teuchos::null)
//      dserror("expect ale block matrix");
//
//    // extract structure and ale submatrices
//    LINALG::SparseMatrix& sig = blocks->Matrix(0,1);  // S_{I\Gamma}
//    LINALG::SparseMatrix& sgg = blocks->Matrix(1,1);  // S_{\Gamma\Gamma}
//    LINALG::SparseMatrix& aig = blocka->Matrix(0,1);  // A_{I\Gamma}
//
//    // get time integration parameters of structure an fluid time integrators
//    // to enable consistent time integration among the fields
//    const double stiparam = StructureField()->TimIntParam();
//    const double ftiparam = FluidField().TimIntParam();
//
//    // some scaling factors for fluid
//    const double scale     = FluidField().ResidualScaling();
//    const double dt        = FluidField().Dt();
//
//    // some often re-used vectors
//    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;
//
//    // ---------- inner structural DOFs
//    /* The following terms are added to the inner structural DOFs of right hand side:
//     *
//     * rhs_firstnewtonstep =
//     *
//     * (1)  S_{I\Gamma} * \Delta d_{\Gamma,p}
//     *
//     * (2)  - dt * S_{I\Gamma} * u^{n}_{\Gamma}
//     *
//     *  Remarks on all terms:
//     * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
//     *
//     */
//
//    // ----------adressing term 1
//    rhs = Teuchos::rcp(new Epetra_Vector(sig.RowMap(),true));
//    sig.Apply(*ddgpred_,*rhs);
//
//    Extractor().AddVector(*rhs,0,f);
//
//    // ----------adressing term 2
//    rhs = Teuchos::rcp(new Epetra_Vector(sig.RowMap(),true));
//    Teuchos::RCP<Epetra_Vector> fveln = FluidField().ExtractInterfaceVeln();
//    sig.Apply(*FluidToStruct(fveln),*rhs);
//    rhs->Scale(-dt);
//
//    Extractor().AddVector(*rhs,0,f);
//
//    // ---------- end of inner structural DOFs
//
//    // we need a temporary vector with the whole fluid dofs where we
//    // can insert the embedded dofrowmap into it
//    Teuchos::RCP<Epetra_Vector> fluidfluidtmp = LINALG::CreateVector(*FluidField().DofRowMap(),true);
//
//    // ---------- inner fluid DOFs
//    /* The following terms are added to the inner fluid DOFs of right hand side:
//     *
//     * rhs_firstnewtonstep =
//     *
//     * (1)  - dt * F^{G}_{I\Gamma} * u^{n}_{\Gamma}
//     *
//     */
//    // ----------addressing term 1
//    if (mmm!=Teuchos::null)
//    {
//      LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
//      rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap(),true));
//      fmig.Apply(*fveln,*rhs);
//      rhs->Scale(-dt);
//
//      rhs = FluidField().Interface()->InsertOtherVector(rhs);
//      xfluidfluidsplitter_->InsertFluidVector(rhs,fluidfluidtmp);
//
//      Extractor().AddVector(*fluidfluidtmp,1,f);
//    }
//
//    // ---------- end of inner fluid DOFs
//
//    // ---------- interface fluid DOFs
//    /* The following terms are added to the interface fluid DOFs of right hand side:
//     *
//     * rhs_firstnewtonstep =
//     *
//     * (1)  - dt * F^{G}_{\Gamma\Gamma} * u^{n}_{\Gamma}
//     *
//     * (2)  - dt * (1-ftiparam)/(1-stiparam) * S_{\Gamma\Gamma} * u^{n}_{\Gamma}
//     *
//     * (3)  + (1-ftiparam)/(1-stiparam) * S_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
//     *
//     * Remarks on all terms:
//     * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
//     *
//     */
//
//    // ----------addressing term 1
//    if (mmm!=Teuchos::null)
//    {
//      LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);
//      rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap(),true));
//      fmgg.Apply(*fveln,*rhs);
//      rhs->Scale(-dt);
//
//      rhs = FluidField().Interface()->InsertFSICondVector(rhs);
//
//      xfluidfluidsplitter_->InsertFluidVector(rhs,fluidfluidtmp);
//      fluidfluidtmp->PutScalar(0.0);
//
//      Extractor().AddVector(*fluidfluidtmp,1,f);
//    }
//
//    // ----------addressing term 2
//    rhs = Teuchos::rcp(new Epetra_Vector(sgg.RowMap(),true));
//    sgg.Apply(*FluidToStruct(fveln),*rhs);
//    rhs->Scale(-dt * (1.-ftiparam) / ((1-stiparam) * scale));
//
//    rhs = StructToFluid(rhs);
//    rhs = FluidField().Interface()->InsertFSICondVector(rhs);
//
//    fluidfluidtmp->PutScalar(0.0);
//    xfluidfluidsplitter_->InsertFluidVector(rhs,fluidfluidtmp);
//
//    Extractor().AddVector(*fluidfluidtmp,1,f);
//
//    // ----------addressing term 3
//    rhs = Teuchos::rcp(new Epetra_Vector(sgg.RowMap(),true));
//    sgg.Apply(*ddgpred_,*rhs);
//    rhs->Scale((1.-ftiparam) / ((1-stiparam) * scale));
//
//    rhs = StructToFluid(rhs);
//    rhs = FluidField().Interface()->InsertFSICondVector(rhs);
//
//    fluidfluidtmp->PutScalar(0.0);
//    xfluidfluidsplitter_->InsertFluidVector(rhs,fluidfluidtmp);
//
//    Extractor().AddVector(*fluidfluidtmp,1,f);
//    // ---------- end of interface fluid DOFs
//
//    // ---------- inner ale DOFs
//    /* The following terms are added to the inner ale DOFs of right hand side:
//     *
//     * rhs_firstnewtonstep =
//     *
//     * (1)  - dt * A_{I\Gamma} * u^{n}_{\Gamma}
//     *
//     */
//    // ----------addressing term 1
//    rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap(),true));
//     Teuchos::RCP<Epetra_Vector> sveln = FluidToStruct(fveln);
//    aig.Apply(*FluidToAleInterface(fveln),*rhs);
//    rhs->Scale(-dt);
//    Extractor().AddVector(*rhs,2,f);
//    // ---------- end of inner ale DOFs
//
//    // -----------------------------------------------------
//    // Now, all contributions/terms to rhs in the first Newton iteration are added.
//
//    // Apply Dirichlet boundary conditions
//    // structure
//    rhs = Extractor().ExtractVector(f,0);
//    Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(),true));
//    LINALG::ApplyDirichlettoSystem(rhs,zeros,*(StructureField()->GetDBCMapExtractor()->CondMap()));
//    Extractor().InsertVector(*rhs,0,f);
//
//    // fluid
//    rhs = Extractor().ExtractVector(f,1);
//    zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(),true));
//    LINALG::ApplyDirichlettoSystem(rhs,zeros,*(FluidField().FluidDirichMaps()));
//    Extractor().InsertVector(*rhs,1,f);
//
//    // ale
//    rhs = Extractor().ExtractVector(f,2);
//    zeros = Teuchos::rcp(new const Epetra_Vector(rhs->Map(),true));
//    LINALG::ApplyDirichlettoSystem(rhs,zeros,*(AleField().GetDBCMapExtractor()->CondMap()));
//    Extractor().InsertVector(*rhs,2,f);
  }

  // -----------------------------------------------------
  // Reset quantities for previous iteration step since they still store values from the last time step
  ddiinc_ = LINALG::CreateVector(*StructureField()->Interface()->OtherMap(),true);
  solipre_ = Teuchos::null;
  ddginc_ = LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(),true);
  solgpre_ = Teuchos::null;
  //fgcur_ = LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(),true);
  // sgicur_ = Teuchos::null;
  // sggcur_ = Teuchos::null;


  // store interface force onto the structure to know it in the next time step as previous force
  // in order to recover the Lagrange multiplier
  // fgpre_ = fgcur_;
  //fgcur_ = StructureField()->Interface()->ExtractFSICondVector(StructureField()->RHS());

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupSystemMatrix()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::SetupSystemMatrix");

  IO::cout << " implement the SetupSystemMatrix() " << IO::endl;

  // extract Jacobian matrices and put them into composite system
  // matrix W
//
  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField()->SystemMatrix();
  if (s==Teuchos::null)
    dserror("expect structure matrix");
  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField()->SystemMatrix();
  if (f==Teuchos::null)
    dserror("expect fluid matrix");

  // scaling factors for fluid
  const double scale     = FluidField()->ResidualScaling(); // inverse of the weighting of the quantities w.r.t the new timestep
  const double timescale = FluidField()->TimeScaling();           // inverse of FSI (1st order, 2nd order) scaling  1/(theta_FSI * dt)
//
//  // get time integration parameters of structure an fluid time integrators
//  // to enable consistent time integration among the fields

  // this is the interpolation weight for quantities from last time step ( alpha_f for genalpha and (1-theta) for OST )
  const double stiparam = StructureField()->TimIntParam();  // alpha_f for genalpha and (1-theta) for OST (weighting of the old timestep n for displacements)
  // this is the interpolation weight for quantities from last time step ( (1-alphaF_) for genalpha and 0.0 for OST and BDF2 )
//  const double ftiparam = FluidField()->TimIntParam();


//
//  // store parts of structural matrix to know them in the next iteration as previous iteration matricesu
//  sgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(s->Matrix(1,0)));
//  sggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(s->Matrix(1,1)));
//

  /*----------------------------------------------------------------------*/
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  /*----------------------------------------------------------------------*/
  // Structure diagonal block
  /*----------------------------------------------------------------------*/

  // TODO: check what this call does?!
  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  s->UnComplete();

  // scale the structure system matrix
  double scaling_S = 1.0/(1.0-stiparam);  // 1/(1-alpha_F) = 1/(1-weight^S_np)
  s->Scale(scaling_S);

  systemmatrix_->Assign(0,0,View,*s); // That's the structure


  /*----------------------------------------------------------------------*/
  // Coupling blocks C_sf, C_fs and C_ss
  /*----------------------------------------------------------------------*/

  Teuchos::RCP<ADAPTER::XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<ADAPTER::XFluidFSI>(FluidField(), true);

  Teuchos::RCP<LINALG::SparseMatrix> C_ss = xfluid->C_Struct_Struct_Matrix();
  Teuchos::RCP<LINALG::SparseMatrix> C_sf = xfluid->C_Struct_Fluid_Matrix();
  Teuchos::RCP<LINALG::SparseMatrix> C_fs = xfluid->C_Fluid_Struct_Matrix();

  C_sf->Scale(scale);            // 1/(theta_f*dt) = 1/weight(t^f_np)
  C_fs->Scale(scale*timescale);  // 1/(theta_f*dt) * 1/(theta_FSI*dt) = 1/weight(t^f_np) * 1/weight(t^FSI_np)
  //C_ss->Scale(scale*timescale);  // 1/(theta_f*dt) * 1/(theta_FSI*dt) = 1/weight(t^f_np) * 1/weight(t^FSI_np)

  // C_ss_block scaled with 1/(theta_f*dt) * 1/(theta_FSI*dt) = 1/weight(t^f_np) * 1/weight(t^FSI_np)
  LINALG::SparseMatrix& C_ss_block = (*systemmatrix_)(0,0);
  C_ss_block.Add(*C_ss, false, scale*timescale, 1.0);  // add the coupling part between fluid and structure on the diagonal block

  // Off-diagonal coupling blocks
  systemmatrix_->Assign(0,1,View,*C_sf);            // assign the off-diagonal coupling block
  systemmatrix_->Assign(1,0,View,*C_fs);            // assign the off-diagonal coupling block


  /*----------------------------------------------------------------------*/
  // Fluid diagonal block
  /*----------------------------------------------------------------------*/

  f->Complete();
  f->Scale(scale); //  1/(theta_f*dt) = 1/weight(t^f_np)

  // finally assign fluid block
  systemmatrix_->Assign(1,1,View,*f);


  // TODO: Apply Dirichlet to fluid and to structure
  //  f->ApplyDirichlet( *(FluidField().FluidDirichMaps()),true);


  /*----------------------------------------------------------------------*/
  // Complete the system matrix
  /*----------------------------------------------------------------------*/

  // done. make sure all blocks are filled.
  systemmatrix_->Complete();

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps)
{

  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  Extractor()->Setup(*fullmap,maps);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
//  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::InitialGuess");
//
//  SetupVector(*ig,
//              StructureField()->InitialGuess(),
//              FluidField().InitialGuess(),
//              AleField().InitialGuess(),
//              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
//  //should we scale the system?
//  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
//  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(fsidyn,"INFNORMSCALING");
//
//  if (scaling_infnorm)
//  {
//    // The matrices are modified here. Do we have to change them back later on?
//
//    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
//    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
//    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
//    A->InvRowSums(*srowsum_);
//    A->InvColSums(*scolsum_);
//    if (A->LeftScale(*srowsum_) or
//        A->RightScale(*scolsum_) or
//        mat.Matrix(0,1).EpetraMatrix()->LeftScale(*srowsum_) or
//        mat.Matrix(0,2).EpetraMatrix()->LeftScale(*srowsum_) or
//        mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_) or
//        mat.Matrix(2,0).EpetraMatrix()->RightScale(*scolsum_))
//      dserror("structure scaling failed");
//
//    A = mat.Matrix(2,2).EpetraMatrix();
//    arowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
//    acolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
//    A->InvRowSums(*arowsum_);
//    A->InvColSums(*acolsum_);
//    if (A->LeftScale(*arowsum_) or
//        A->RightScale(*acolsum_) or
//        mat.Matrix(2,0).EpetraMatrix()->LeftScale(*arowsum_) or
//        mat.Matrix(2,1).EpetraMatrix()->LeftScale(*arowsum_) or
//        mat.Matrix(0,2).EpetraMatrix()->RightScale(*acolsum_) or
//        mat.Matrix(1,2).EpetraMatrix()->RightScale(*acolsum_))
//      dserror("ale scaling failed");
//
//    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b,0);
//    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b,2);
//
//    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0))
//      dserror("structure scaling failed");
//    if (ax->Multiply(1.0, *arowsum_, *ax, 0.0))
//      dserror("ale scaling failed");
//
//    Extractor().InsertVector(*sx,0,b);
//    Extractor().InsertVector(*ax,2,b);
//  }
}

/*----------------------------------------------------------------------*
 |  map containing the dofs with Dirichlet BC
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> FSI::MonolithicXFEM::CombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map > scondmap = StructureField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map > fcondmap = FluidField()->GetDBCMapExtractor()->CondMap();

  Teuchos::RCP<Epetra_Map> condmap = LINALG::MergeMap(scondmap, fcondmap, false);

  return condmap;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
//  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
//  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(fsidyn,"INFNORMSCALING");
//
//  if (scaling_infnorm)
//  {
//    Teuchos::RCP<Epetra_Vector> sy = Extractor().ExtractVector(x,0);
//    Teuchos::RCP<Epetra_Vector> ay = Extractor().ExtractVector(x,2);
//
//    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0))
//      dserror("structure scaling failed");
//    if (ay->Multiply(1.0, *acolsum_, *ay, 0.0))
//      dserror("ale scaling failed");
//
//    Extractor().InsertVector(*sy,0,x);
//    Extractor().InsertVector(*ay,2,x);
//
//    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b,0);
//    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b,2);
//
//    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0))
//      dserror("structure scaling failed");
//    if (ax->ReciprocalMultiply(1.0, *arowsum_, *ax, 0.0))
//      dserror("ale scaling failed");
//
//    Extractor().InsertVector(*sx,0,b);
//    Extractor().InsertVector(*ax,2,b);
//
//    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
//    srowsum_->Reciprocal(*srowsum_);
//    scolsum_->Reciprocal(*scolsum_);
//    if (A->LeftScale(*srowsum_) or
//        A->RightScale(*scolsum_) or
//        mat.Matrix(0,1).EpetraMatrix()->LeftScale(*srowsum_) or
//        mat.Matrix(0,2).EpetraMatrix()->LeftScale(*srowsum_) or
//        mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_) or
//        mat.Matrix(2,0).EpetraMatrix()->RightScale(*scolsum_))
//      dserror("structure scaling failed");
//
//    A = mat.Matrix(2,2).EpetraMatrix();
//    arowsum_->Reciprocal(*arowsum_);
//    acolsum_->Reciprocal(*acolsum_);
//    if (A->LeftScale(*arowsum_) or
//        A->RightScale(*acolsum_) or
//        mat.Matrix(2,0).EpetraMatrix()->LeftScale(*arowsum_) or
//        mat.Matrix(2,1).EpetraMatrix()->LeftScale(*arowsum_) or
//        mat.Matrix(0,2).EpetraMatrix()->RightScale(*acolsum_) or
//        mat.Matrix(1,2).EpetraMatrix()->RightScale(*acolsum_))
//      dserror("ale scaling failed");
//
//  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupVector(Epetra_Vector &f,
                                      Teuchos::RCP<const Epetra_Vector> sv,
                                      Teuchos::RCP<const Epetra_Vector> fv,
                                      double fluidscale)
{
  // get time integration parameters of structure an fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
//  const double ftiparam = FluidField()->TimIntParam();


  /*----------------------------------------------------------------------*/
  // structure part of rhs vector
  /*----------------------------------------------------------------------*/

  // add fluid interface coupling values to structure vector
  // sv: structure dofs

  Teuchos::RCP<ADAPTER::XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<ADAPTER::XFluidFSI>(FluidField(), true);

  // this vector is based on the boundary dis which is part of the structure dis
  Teuchos::RCP<Epetra_Vector> rhs_C_s = xfluid->RHS_Struct_Vec();

  Teuchos::RCP<Epetra_Vector> sv_cast = Teuchos::rcp_const_cast< Epetra_Vector >(sv);

  // scaling for rhs_C_s
  double scaling_S = 1.0/(1.0-stiparam);

  // adding rh_C_s_ to structure residual
  for (int iter=0; iter<rhs_C_s->MyLength();++iter)
  {
    // gid of the rhs entry
    int rhsdgid = rhs_C_s->Map().GID(iter);
    // is this entry a row entry on this proc?
    if (rhs_C_s->Map().MyGID(rhsdgid) == false) dserror("rhs_C_s should be distributed to all processors");
    // is this entry also available as row on the underlying structural dis, this has to be!
    if (sv_cast->Map().MyGID(rhsdgid))
      (*sv_cast)[sv_cast->Map().LID(rhsdgid)]=(*sv_cast)[sv_cast->Map().LID(rhsdgid)] + (*rhs_C_s)[rhs_C_s->Map().LID(rhsdgid)] * scaling_S;
    else dserror("boundary (cut) dof not on structure discret available! Something wrong!");
  }

  //TODO: how to do this here ?
  // CONTRIBUTION FROM LAGRANGE_MUTLIPLIER FROM PREVIOUS TIME STEP
  // add contribution of Lagrange multiplier from previous time step
//    if (lambda_ != Teuchos::null)
//    {
//      Teuchos::RCP<Epetra_Vector> lambdaglobal = FluidField().Interface()->InsertFSICondVector(StructToFluid(lambda_));
//      modfv->Update((-ftiparam+stiparam*(1.0-ftiparam)/(1.0-stiparam))/fluidscale, *lambdaglobal, 1.0);
//    }

  // when the structure part is filled we can insert it into the global vector
  Extractor()->InsertVector(*sv_cast,0,f);

  //TODO: where to add the structure DBC?


  //------------------------------------------------------------------
  // fluid part of rhs vector
  //------------------------------------------------------------------
  if (fluidscale!=0)
  {
//   // modfv: whole embedded fluid map but entries at fsi dofs
//    Teuchos::RCP<Epetra_Vector> modfv = FluidField().Interface()->InsertFSICondVector(StructToFluid(scv));
//
//    // modfv = modfv * 1/fluidscale * d/b
//    modfv->Scale( (1./fluidscale)*(1.0-ftiparam)/(1.0-stiparam));
//
//
//    // we need a temporary vector with the whole fluid dofs where we
//    // can insert veln which has the embedded dofrowmap into it
//    Teuchos::RCP<Epetra_Vector> fluidtmp = LINALG::CreateVector(*FluidField()->DofRowMap(),true);
//    xfluidfluidsplitter_->InsertFluidVector(modfv,fluidtmp);
//
//    Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(fluidtmp->Map(),true));
//    LINALG::ApplyDirichlettoSystem(fluidtmp,zeros,*(FluidField().FluidDirichMaps()));

    // all fluid dofs
    Teuchos::RCP<Epetra_Vector> fvfluid = Teuchos::rcp_const_cast< Epetra_Vector >(fv);

    //TODO: do we need this scaling?

    fvfluid->Scale(fluidscale); // scale with FluidField()->ResidualScaling()

//    // adding modfv to fvfluid
//    fvfluid->Update(1.0,*fluidtmp,1.0);

    // no fluid scaling, then just insert the fluid vector into the global vector
    Extractor()->InsertVector(*fvfluid,1,f);
  }
  else
  {
    // no fluid scaling, then just insert the fluid vector into the global vector
    Extractor()->InsertVector(*fv,1,f);
  }

}



/*----------------------------------------------------------------------*
 | time loop of the monolithic system                      schott 08/13 |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Timeloop()
{
  // time loop
  while (NotFinished())
  {
    // counter and print header
    // predict solution of both field (call the adapter)
    PrepareTimeStep();

    // outer iteration loop when active fluid dofsets change
    // calls inner Newton-Raphson iterations within each outer iteration
    Solve();

//    // calculate stresses, strains, energies
//    PrepareOutput();

    // update all single field solvers
    Update();

    // write output to screen and files
    Output();


  }  // NotFinished
}  // TimeLoop()


/*----------------------------------------------------------------------*
 | prepare the time step for fluid and structure           schott 08/13 |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::PrepareTimeStep()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::PrepareTimeStep");

  IncrementTimeAndStep();

  PrintHeader();

  //--------------------------------------------
  // Structure PrepareTimeStep
  //--------------------------------------------
  // * apply structural predictor
  // * apply Dirichlet conditions and
  // * print the residual based on the structural predictor
  StructureField()->PrepareTimeStep();

  //--------------------------------------------
  // Fluid PrepareTimeStep
  //--------------------------------------------
  // * set time integrator
  // * set time parameters
  FluidField()    ->PrepareTimeStep();

}


/*----------------------------------------------------------------------*
 | outer iteration loop when active fluid dofsets change   schott 08/13 |
 | calls inner Newton-Raphson iterations within each outer iteration    |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Solve()
{

  iter_outer_ = 1; // allow restarts of the Newton scheme in case of changing fluid maps

  // We want to make sure, that the outer loop is entered at least once!
  // We exit the outer loop if either the inner Newton loop is converged
  // (checked in Converged() method) OR the maximum number of outer iterations is exceeded!
  while ( (iter_outer_ == 1) or ( iter_outer_ <= itermax_outer_) )
  {
    if(NewtonFull()) // stop since the main inner Newton loop converged
    {
      if( Comm().MyPID() == 0)
      {
        IO::cout << "-------------------------------------- Outer loop finished with converged NewtonLoop ---------------------------------------" << IO::endl;
      }
      break;
    }
    else
    {
      if( Comm().MyPID() == 0)
      {
        IO::cout << "------------------------------------------ Restart NewtonLoop - DOF-sets changed -------------------------------------------" << IO::endl;
      }
    }

    iter_outer_ += 1;
  }

  if(iter_outer_ > itermax_outer_)
  {
    if( Comm().MyPID() == 0)
    {
      IO::cout << "------------------------------ Outer loop did not converge ! DOF-sets have changed too often -------------------------------" << IO::endl;
    }
  }
}


/*----------------------------------------------------------------------*
 | inner iteration loop (Newton-Raphson iteration)         schott 08/13 |
 | return if converged or not (in case of changing fluid maps)          |
 *----------------------------------------------------------------------*/
bool FSI::MonolithicXFEM::NewtonFull()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::NewtonFull::NewtonFull()");

  if (Comm().MyPID() == 0) std::cout << "FSI::MonolithicXFEM::NewtonFull()" << std::endl;

  //--------------------------------------------------------
  // Extract the current interface displacement
  // For a new timestep we predicted the structural displacement in PerpareTimeStep()
  // REMARK: be aware of using const Dis predictor!
  //         This results in zero disp_incr and u^n+1 = -u^n for second order interface
  //                                        and u^n+1 = 0 for second order interface
  Teuchos::RCP<Epetra_Vector> idispnp = StructureField()->ExtractInterfaceDispnp();

  //--------------------------------------------------------
  // Apply the current interface displacement idispnp to the boundary discretization of the fluid field
  FluidField()->ApplyMeshDisplacement(idispnp);


  /*----------------------------------------------------------------------
  Prepare the fluid field for a new start of the Newton scheme (for a new timestep or after changing fluid maps):
  * cut with the current interface position and create new state class ( + vectors + fluid-sysmat) and
  * perform the xfluid time integration

  * set SetOldPartOfRighthandside
  * set Dirichlet values in fluid velnp vector and compute Neumann loads

  REMARK:
  * applying Dirichlet values to rhs and to system has to be done
  * for the global system after Evaluate routines have been called for each field
  * in xfluid we set Dirichlet values into velnp, therefore the stepinc used in evaluate adds zeros for Dirichlet values
    not to modify the set Dirichlet values
  ----------------------------------------------------------------------*/

  FluidField()->PrepareSolve();

  // setup a new system since the size of the system has changed when NewtonFull is called
  SetupNewSystem();

  /*----------------------------------------------------------------------
  Initialization
  ----------------------------------------------------------------------*/
  //Iteration counter
  iter_ = 1;

  //Create a new Epetra_Vector with a given row map and initialize it!
  //DofRowMap() contains all DOFS for the monolithic system w.r.t the current interface position and returns the
  //DofRowMap of the blockrowdofmap_, a MapExtractor object for a DOF map, split into blocks.

  //Increment sum vector (that's the total increment w.r.t the old time step t^n)
  x_sum_ = LINALG::CreateVector(*DofRowMap(), true);

  // Solution vector for linear solve = iteration increment Delta x = x^n+1_i+1 - x^n+1_i
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);

  // Residual vector
  rhs_ = LINALG::CreateVector(*DofRowMap(), true);

  zeros_ = LINALG::CreateVector(*DofRowMap(), true);

  //TODO: do we really need this?
  // Flag for special treatment of RHS-setup for i==0
  firstcall_ = true;

//  /*----------------------------------------------------------------------
//  Extract predictor increments
//  ----------------------------------------------------------------------*/
//  //Increment of structural interface displacement --> structural predictor!!
//  ddgpred_=Teuchos::rcp(new Epetra_Vector(*(StructureField()->ExtractInterfaceDispnp())));
//  ddgpred_->Update(-1.0,*StructureField()->ExtractInterfaceDispn(),1.0);
//
//  /*----------------------------------------------------------------------*/
//  //Initialize the increment vectors, they are updated in Evaluate(...)->ExtractFieldVectors(...)
//  //at every Newton iteration!
//
//  //Initialization for 1st Newton call
//  //structural interface predictor
//  ddginc_    = Teuchos::rcp(new Epetra_Vector(*ddgpred_));
//  ddialeinc_ = Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->OtherMap()),true);
//  duiinc_    = Teuchos::rcp(new Epetra_Vector(*Extractor().ExtractVector(iterinc_,1)),true);


  //We want to make sure, that the loop is entered at least once!
  //We exit the loop if either the convergence criteria are met (checked in
  //Converged() method) OR the maximum number of iterations is exceeded!
  while ( (iter_ == 1) or ( (not Converged()) and (iter_ <= itermax_) ) )
  {
    std::cout << "Evaluate-Call " << "iter_ " << iter_ << "/" << itermax_ << std::endl;


    //Evaluate()- call
    //This function calls evaluate methods for all fields,
    //assembles and transfers data
    Evaluate(iterinc_);

    //Did the fluid map change after increment?
    //Problem here: The fluid FSI DOFs are removed from the system
    //they are not in FSI algorithm's global DOF-map anymore -
    //the increment vector map is free of FSI DOFs.
    //But as they are still available in the internal fluid field DOF map,
    //a comparison of the 2 maps would always return 'false'.
    //Solution:
    //We need a splitter object to remove fluid FSI DOFs from the fluid DOF map
    //before the maps are compared!
    //We want to control, if the fluid map itself has really changed!
//
//    Teuchos::RCP<Epetra_Map> fsimap=Teuchos::rcp(new Epetra_Map(*FluidField().Interface()->FSICondMap()));
//    Teuchos::RCP<LINALG::MapExtractor> ffsextractor=Teuchos::rcp(new LINALG::MapExtractor(*(FluidField().DofRowMap()),fsimap,true));
//    Teuchos::RCP<const Epetra_Map> innerfluidmap=ffsextractor->OtherMap();
//
//    bool isnewfluidmap=(innerfluidmap->SameAs(Extractor().ExtractVector(iterinc_,1)->Map()));
//


    /*----------------------------------------------------------------------*/
    //Build the linear system
    //J^{n+1}(x_i) \Delta^{n+1}_{x_i+1}=-r^{n+1}(x_i)
    //i: Newton iteration counter
    //J: Jacobian
    //r: RHS-vector
    /*----------------------------------------------------------------------*/
    SetupSystemMatrix();

    if ( not systemmatrix_->Filled() )
      dserror("Unfilled system matrix! Fatal error!");

    // Create the RHS consisting of the field residuals and
    // predictor-related terms added at the first Newton iteration only.
    SetupRHS(*rhs_,firstcall_);

    //Solver call
    LinearSolve();

    //Adapt solver tolerance (important if Aztec solver is picked)
    solver_->ResetTolerance();

    //Build residual and incremental norms,
    //count the DOFs!
    BuildCovergenceNorms();

    //Give some output
    PrintNewtonIter();

    //Increment loop index
    iter_+=1;

    //firstcall_ = false;
  }//End of Newton loop!

  //After the loop exit, the iteration counter is 1 higher than the true no. of
  //iterations! Correct that:
  iter_-=1;

//  /*----------------------------------------------------------------------*/
//  //Compute the increments needed for recovery of Lagrange Multiplier!
//  //After the last Newton iteration, the increments are not updated.
//  //We need the last increment for the recovery of lambda.
//  /*----------------------------------------------------------------------*/
//  //Fluid
//  duiinc_->Update(1.0,*Extractor().ExtractVector(iterinc_,1),0.0);
//  //Structure
//  Teuchos::RCP<Epetra_Vector> ddinc=Extractor().ExtractVector(iterinc_,0);
//  ddginc_->Update(1.0,*StructureField()->Interface()->ExtractFSICondVector(ddinc),0.0);
//  //ALE
//  ddialeinc_->Update(1.0,*Extractor().ExtractVector(iterinc_,2),0.0);

  if ( Converged() )
  {
    if(Comm().MyPID() == 0)
    {
      IO::cout << "-------------------------------------------------------Newton Converged ! --------------------------------------------------" << IO::endl;
    }
    return true;
  }
  else if ( iter_ >= itermax_ )
  {
    if(Comm().MyPID() == 0)
    {
      IO::cout << "----------------------------------------- Newton unconverged in ITEMAX iterations ! ----------------------------------------" << IO::endl;
    }
    return false;
  }

  return false;
}  // NewtonFull()



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Evaluate(Teuchos::RCP<const Epetra_Vector> iterinc)
{
   TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::Evaluate");
   Teuchos::RCP<const Epetra_Vector> sx;
   Teuchos::RCP<const Epetra_Vector> fx;

   // TODO: do we have to set x_sum for the 1st iteration due to the predictor?!

   double inc_norm = 0.0;
   iterinc->Norm2(&inc_norm);

   if( iter_ == 1 and fabs(inc_norm) > 1e-12) dserror("non-zero increment in the first Newton-Call, what to do with x_sum?");


   // structure and fluid fields expects the step increment (x^n+1_i+1 - x^n).
   // So we add all of the increments together to build the step increment.
   //
   // The update of the latest step increment with iteration increments:
   // x^n+1_i+1 = x^n+1_i + iterinc with x the current step increment

   x_sum_->Update(1.0,*iterinc,1.0);

   // sx contains the current step increment w.r.t. t^n for the structure block
   // fx contains the current step increment w.r.t. t^n for the fluid block
   ExtractFieldVectors(x_sum_,sx,fx);

   if (sdbg_!=Teuchos::null)
   {
     sdbg_->NewIteration();
     sdbg_->WriteVector("x",*StructureField()->Interface()->ExtractFSICondVector(sx));
   }


   //-------------------
   // structure field
   //-------------------
   {
     // structural field
     Epetra_Time ts(Comm());

     // call the structure evaluate with the current time step increment
     StructureField()->Evaluate(sx);

     IO::cout  << "structure time: " << ts.ElapsedTime() << IO::endl;
   }


   // ------------------------------------------------------------------
   // set the current interface displacement for the fluid field
   // cut, perform time integration and
   // ------------------------------------------------------------------

   if(iter_ > 1)
   {

     // get the current structural displacement
     Teuchos::RCP<Epetra_Vector> idispnp = StructureField()->ExtractInterfaceDispnp();

     //--------------------------------------------------------
     // apply the current interface displacement idispnp to the boundary discretization of the fluid field
     FluidField()->ApplyMeshDisplacement(idispnp);

     //--------------------------------------------------------
     // cut with new interface position
     FluidField()->CutAndSetStateVectors(true); // tell the xfluid time integration that it is a new (not first) step of the Newton scheme

   }


   //--------------------------------------------------------
   // set current interface velocities
   //--------------------------------------------------------
   {

     //--------------------------------------------------------
     // new vector for interface velocities
     Teuchos::RCP<Epetra_Vector> ivelnp = LINALG::CreateVector( *FluidField()->Interface()->FSICondMap(),true);

     // now ivelnp contains the interface displacements
     ivelnp->Update(1.0, *StructureField()->ExtractInterfaceDispnp(), -1.0, *StructureField()->ExtractInterfaceDispn(), 0.0 );

     // now idispnp is converted to the interface velocity
     FluidField()->DisplacementToVelocity( ivelnp );

     ivelnp->Update( 1.0, *(FluidField()->ExtractInterfaceVeln()), 1.0);

     // set ivelnp in Xfluid
     FluidField()->ApplyInterfaceVelocities(ivelnp);
   }

   // ------------------------------------------------------------------
   // Call all fields evaluate method and assemble rhs and matrices
   // ------------------------------------------------------------------

   //-------------------
   // fluid field
   //-------------------
   {
     // fluid field
     Epetra_Time tf(Comm());

     // call the fluid evaluate with the current time step increment
     FluidField()->Evaluate(fx);

     IO::cout << "fluid time : " << tf.ElapsedTime() << IO::endl;
   }

}


/*----------------------------------------------------------------------
 |                                                        schott 08/13 |
 |   extract the two field vectors from a given composed vector        |
 ----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector>  x,
                                              Teuchos::RCP<const Epetra_Vector>& sx,
                                              Teuchos::RCP<const Epetra_Vector>& fx )
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::ExtractFieldVectors");

  /*----------------------------------------------------------------------*/
  //Process structure unknowns
  /*----------------------------------------------------------------------*/
  // Extract whole structure field vector
  sx=Extractor()->ExtractVector(x,0);

//  //Structural part of FSI interface
//  Teuchos::RCP<Epetra_Vector> scx=StructureField()->Interface()->ExtractFSICondVector(sx);

  /*----------------------------------------------------------------------*/
  //Process fluid unknowns
  /*----------------------------------------------------------------------*/
  // Extract vector of fluid unknowns from x
  fx=Extractor()->ExtractVector(x,1);

}


/*----------------------------------------------------------------------*
 | solve linear FSI system                                 schott 07/13 |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::LinearSolve()
{
  if (Comm().MyPID() == 0)
    std::cout << " FSI::MonolithicXFEM::LinearSolve()" <<  std::endl;

  // Solve for inc_ = [disi_,tempi_]
  // Solve K_Teffdyn . IncX = -R  ===>  IncX_{n+1} with X=[d,(u,p)]
  // \f$x_{i+1} = x_i + \Delta x_i\f$
  if (solveradapttol_ and (iter_ > 1))
  {
    double worst = normrhs_;
    double wanted = tolrhs_;
    solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
  }

  // Dirichlet boundary conditions are already applied to FSI system, i.e.
  // FSI system is prepared for solve, i.e. FSI systemmatrix, FSI rhs, FSI inc
  // --> in PrepareSystemForNewtonSolve(): done for rhs and diagonal blocks
  // --> in SetupSystemMatrix() done for off-diagonal blocks k_st, k_ts

  // apply Dirichlet BCs to system of equations
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

  // default: use block matrix
  if (merge_fsi_blockmatrix_ == false)
  {
    if (Comm().MyPID() == 0)
    {
      std::cout << " DBC applied to TSI system on proc" << Comm().MyPID() <<  std::endl;
    }

    // Infnormscaling: scale system before solving
    ScaleSystem(*systemmatrix_,*rhs_);

    // solve the problem, work is done here!
    solver_->Solve(
        systemmatrix_->EpetraOperator(),
        iterinc_,
        rhs_,
        true,
        iter_==1
    );

    // Infnormscaling: unscale system after solving
    UnscaleSolution(*systemmatrix_,*iterinc_,*rhs_);
  }  // use block matrix
  else // (merge_fsi_blockmatrix_ == true)
  {
    // merge blockmatrix to SparseMatrix and solve
    Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

    // apply combined Dirichlet to System
    LINALG::ApplyDirichlettoSystem(
      sparse,
      iterinc_,
      rhs_,
      Teuchos::null,
      zeros_,
      *CombinedDBCMap()
      );

    // standard solver call
    solver_->Solve(
               sparse->EpetraOperator(),
               iterinc_,
               rhs_,
               true,
               iter_==1
               );
  }  // MergeBlockMatrix

  if (Comm().MyPID() == 0) { std::cout << " Solved" <<  std::endl; }

}  // LinearSolve()



/*----------------------------------------------------------------------*
 | create linear solver                                   wiesner 07/11 |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::CreateLinearSolver()
{
  // get the solver number used for linear TSI solver
  const int linsolvernumber = fsidyn_.get<int>("LINEAR_SOLVER");
  // check if the TSI solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for monolithic TSI. Please set LINEAR_SOLVER in TSI DYNAMIC to a valid number!");

  // get parameter list of structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // use solver blocks for structure
  // get the solver number used for structural solver
  const int slinsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (slinsolvernumber == (-1))
    dserror("no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

  // get parameter list of thermal dynamics
  const Teuchos::ParameterList& tdyn = DRT::Problem::Instance()->ThermalDynamicParams();
  // use solver blocks for temperature (thermal field)
  // get the solver number used for thermal solver
  const int tlinsolvernumber = tdyn.get<int>("LINEAR_SOLVER");
  // check if the TSI solver has a valid solver number
  if (tlinsolvernumber == (-1))
    dserror("no linear solver defined for thermal field. Please set LINEAR_SOLVER in THERMAL DYNAMIC to a valid number!");

  // get solver parameter list of linear TSI solver
  const Teuchos::ParameterList& tsisolverparams
    = DRT::Problem::Instance()->SolverParams(linsolvernumber);

  const int solvertype
    = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(
        tsisolverparams,
        "SOLVER"
        );

  if ( (solvertype != INPAR::SOLVER::aztec_msr) and (solvertype != INPAR::SOLVER::belos) )
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << " Note: the BGS2x2 preconditioner now "                  << std::endl;
    std::cout << " uses the structural solver and thermal solver blocks"  << std::endl;
    std::cout << " for building the internal inverses"                    << std::endl;
    std::cout << " Remove the old BGS PRECONDITIONER BLOCK entries "      << std::endl;
    std::cout << " in the dat files!"                                     << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    dserror("aztec solver expected");
  }
  const int azprectype
    = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(
        tsisolverparams,
        "AZPREC"
        );

  // plausibility check
  switch (azprectype)
  {
  case INPAR::SOLVER::azprec_BGS2x2:
    break;
  case INPAR::SOLVER::azprec_BGSnxn:
  case INPAR::SOLVER::azprec_TekoSIMPLE:
  {
#ifdef HAVE_TEKO
    // check if structural solver and thermal solver are Stratimikos based (Teko expects stratimikos)
    int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(slinsolvernumber), "SOLVER");
    if ( (solvertype != INPAR::SOLVER::stratimikos_amesos) and
         (solvertype != INPAR::SOLVER::stratimikos_aztec) and
         (solvertype != INPAR::SOLVER::stratimikos_belos)
       )
    dserror("Teko expects a STRATIMIKOS solver object in STRUCTURE SOLVER");

    solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(tlinsolvernumber), "SOLVER");
    if ( (solvertype != INPAR::SOLVER::stratimikos_amesos) and
         (solvertype != INPAR::SOLVER::stratimikos_aztec) and
         (solvertype != INPAR::SOLVER::stratimikos_belos)
       )
      dserror("Teko expects a STRATIMIKOS solver object in thermal solver %3d",tlinsolvernumber);
#else
    dserror("Teko preconditioners only available with HAVE_TEKO flag for TRILINOS_DEV (>Q1/2011)");
#endif
    break;
  }
  default:
    dserror("Block Gauss-Seidel BGS2x2 preconditioner expected");
    break;
  }

  solver_ = Teuchos::rcp(new LINALG::Solver(
                               tsisolverparams,
                               // ggfs. explizit Comm von STR wie lungscatra
                               Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()
                               )
              );

  // use solver blocks for structure and temperature (thermal field)
  const Teuchos::ParameterList& ssolverparams = DRT::Problem::Instance()->SolverParams(slinsolvernumber);
  const Teuchos::ParameterList& tsolverparams = DRT::Problem::Instance()->SolverParams(tlinsolvernumber);

  solver_->PutSolverParamsToSubParams("Inverse1", ssolverparams);
  solver_->PutSolverParamsToSubParams("Inverse2", tsolverparams);

  // prescribe rigid body modes
  StructureField()->Discretization()->ComputeNullSpaceIfNecessary(
                                        solver_->Params().sublist("Inverse1")
                                        );
  FluidField()->Discretization()->ComputeNullSpaceIfNecessary(
                                     solver_->Params().sublist("Inverse2")
                                     );
}  // CreateLinearSolver()



/*----------------------------------------------------------------------*
 | check convergence of Newton iteration (public)          schott 08/13 |
 *----------------------------------------------------------------------*/
bool FSI::MonolithicXFEM::Converged()
{
  // check for single norms
  bool convinc  = false;
  bool convfres = false;

  // residual increments
  switch (normtypeinc_)
  {
    case INPAR::FSI::convnorm_abs:
      convinc = norminc_ < tolinc_;
      break;
    case INPAR::FSI::convnorm_rel:
      convinc = (((normstrincL2_/ns_)     < TOL_DIS_INC_L2_)    and
                 ((normstrincInf_)        < TOL_DIS_INC_INF_)   and
                 ((normflvelincL2_/nfv_)  < TOL_VEL_INC_L2_)    and
                 ((normflvelincInf_)      < TOL_VEL_INC_INF_)   and
                 ((normflpresincL2_/nfp_) < TOL_PRE_INC_L2_)    and
                 ((normflpresincInf_)     < TOL_PRE_INC_INF_));
      break;
    case INPAR::FSI::convnorm_mix:
      dserror("not implemented!");
      break;
  default:
      dserror("Cannot check for convergence of residual values!"); break;
  }

  // structural and fluid residual forces
  switch (normtypefres_)
  {
  case INPAR::FSI::convnorm_abs:
    convfres = normrhs_ < tolfres_;
    break;
  case INPAR::FSI::convnorm_rel:
    convfres =  (((normstrrhsL2_/ns_)     < TOL_DIS_RES_L2_)    and
                 ((normstrrhsInf_)        < TOL_DIS_RES_INF_)   and
                 ((normflvelrhsL2_/nfv_)  < TOL_VEL_RES_L2_)    and
                 ((normflvelrhsInf_)      < TOL_VEL_RES_INF_)   and
                 ((normflpresrhsL2_/nfp_) < TOL_PRE_RES_L2_)    and
                 ((normflpresrhsInf_)     < TOL_PRE_RES_INF_));
    break;
  case INPAR::FSI::convnorm_mix:
    dserror("not implemented!");
    break;
  default:
    dserror("Cannot check for convergence of residual forces!"); break;
  }

  // combined
  bool conv = false;
  if (combincfres_==INPAR::FSI::bop_and)
     conv = convinc and convfres;
   else
     dserror("Something went wrong!");

  return conv;
}  // Converged()


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Update()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::Update");

  //TODO:
  // recover Lagrange multiplier \lambda_\Gamma at the interface at the end of each time step
  // (i.e. condensed forces onto the structure) needed for rhs in next time step
  //RecoverLagrangeMultiplier();

  StructureField()->Update();
  FluidField()->Update();

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Output()
{
  StructureField()->Output();
//  // output Lagrange multiplier
//  {
//    Teuchos::RCP<Epetra_Vector> lambdafull = StructureField()->Interface()->InsertFSICondVector(lambda_);
//
//    const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
//    const int uprestart = fsidyn.get<int>("RESTARTEVRY");
//    const int upres = fsidyn.get<int>("UPRES");
//    if((uprestart != 0 && FluidField().Step() % uprestart == 0) || FluidField().Step() % upres == 0)
//      StructureField()->DiscWriter()->WriteVector("fsilambda", lambdafull);
//  }

  FluidField()->Output();
  FluidField()->LiftDrag();
  if (StructureField()->GetConstraintManager()->HaveMonitor())
  {
    StructureField()->GetConstraintManager()->ComputeMonitorValues(StructureField()->Dispnp());
    if(Comm().MyPID() == 0)
      StructureField()->GetConstraintManager()->PrintMonitorValues();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::ReadRestart(int step)
{
  IO::cout << "FSI::MonolithicXFEM::ReadRestart still has to be checked!!!" << IO::endl;
  dserror("check the implementation now!");
//  // read Lagrange multiplier
//  {
//    Teuchos::RCP<Epetra_Vector> lambdafull = Teuchos::rcp(new Epetra_Vector(*StructureField()->DofRowMap(),true));
//    IO::DiscretizationReader reader = IO::DiscretizationReader(StructureField()->Discretization(),step);
//    reader.ReadVector(lambdafull, "fsilambda");
//    lambda_ = StructureField()->Interface()->ExtractFSICondVector(lambdafull);
//  }
//
  StructureField()->ReadRestart(step);
  FluidField()->ReadRestart(step);

  SetTimeStep(FluidField()->Time(),FluidField()->Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupNewSystem()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::SetupNewSystem()");

  /*----------------------------------------------------------------------
  Create a combined map for Structure/Fluid-DOFs all in one!
  ----------------------------------------------------------------------*/
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;

  //Append the structural DOF map
  vecSpaces.push_back(StructureField()->DofRowMap());

  //Append the background fluid DOF map
  vecSpaces.push_back(FluidField()->DofRowMap());

  //solid maps empty??
  if (vecSpaces[0]->NumGlobalElements() == 0)
    dserror("No solid equations. Panic.");

  //fluid maps empty??
  if (vecSpaces[1]->NumGlobalElements() == 0)
    dserror("No fluid equations. Panic.");

  //The vector is complete, now fill the system's global block row map
  //with the maps previously set together!
  SetDofRowMaps(vecSpaces);

  /*----------------------------------------------------------------------*/
  // create block system matrix

  switch(linearsolverstrategy_)
  {
  case INPAR::FSI::PreconditionedKrylov:
  {
    /*----------------------------------------------------------------------*/
    // initialise XFSI-systemmatrix_
    systemmatrix_
      = Teuchos::rcp(
          new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                *Extractor(),
                *Extractor(),
                81+54, // maybe adapt this number for XFluid!
                false,
                true
                )
          );
    break;
  }
//  case INPAR::FSI::FSIAMG:
//    systemmatrix_ = Teuchos::rcp(new OverlappingBlockMatrixFSIAMG(
//                                   Extractor(),
//                                   *StructureField(),
//                                   FluidField(),
//                                   AleField(),
//                                   false,
//                                   DRT::INPUT::IntegralValue<int>(fsidyn,"SYMMETRICPRECOND"),
//                                   blocksmoother,
//                                   schuromega,
//                                   pcomega,
//                                   pciter,
//                                   spcomega,
//                                   spciter,
//                                   fpcomega,
//                                   fpciter,
//                                   apcomega,
//                                   apciter,
//                                   DRT::INPUT::IntegralValue<int>(fsidyn,"FSIAMGANALYZE"),
//                                   linearsolverstrategy_,
//                                   DRT::Problem::Instance()->ErrorFile()->Handle()));
//    break;
  default:
    dserror("Unsupported type of monolithic solver");
    break;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Newton()
{

//  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::Newton()");
//
//  // initialise equilibrium loop
//  iter_ = 1;
//
//  x_sum_ = LINALG::CreateVector(*DofRowMap(),true);
//  x_sum_->PutScalar(0.0);
//
//  // incremental solution vector with length of all FSI dofs
//  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
//  iterinc_->PutScalar(0.0);
//
//  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
//  zeros_->PutScalar(0.0);
//
//  // residual vector with length of all FSI dofs
//  rhs_ = LINALG::CreateVector(*DofRowMap(), true);
//  rhs_->PutScalar(0.0);
//
//  firstcall_ = true;
//
//  // equilibrium iteration loop (loop over k)
//  while ( (iter_ ==  1) or((not Converged()) and (iter_ <= itermax_)) )
//  {
//    // compute residual forces #rhs_ and tangent #tang_
//    // build linear system stiffness matrix and rhs/force
//    // residual for each field
//
//    //TODO: implement Evaluate
//    Evaluate(iterinc_);
//
//    if (not FluidField().DofRowMap()->SameAs(Extractor().ExtractVector(iterinc_,1)->Map()))
//    {
//      IO::cout << " New Map!! " << IO::endl;
//      // save the old x_sum
//      Teuchos::RCP<Epetra_Vector> x_sum_n =  LINALG::CreateVector(*DofRowMap(), true);
//      *x_sum_n = *x_sum_;
//      Teuchos::RCP<const Epetra_Vector> sx_n;
//      Teuchos::RCP<const Epetra_Vector> ax_n;
//      sx_n = Extractor().ExtractVector(x_sum_n,0);
//      ax_n = Extractor().ExtractVector(x_sum_n,2);
//
//      SetupNewSystem();
//      xfluidfluidsplitter_ = FluidField().XFluidFluidMapExtractor();
//      rhs_ = LINALG::CreateVector(*DofRowMap(), true);
//      iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
//      zeros_ = LINALG::CreateVector(*DofRowMap(), true);
//      x_sum_ = LINALG::CreateVector(*DofRowMap(),true);
//
//      // build the new iter_sum
//      Extractor().InsertVector(sx_n,0,x_sum_);
//      Extractor().InsertVector(FluidField().Stepinc(),1,x_sum_);
//      Extractor().InsertVector(ax_n,2,x_sum_);
//      nf_ = (*(FluidField().RHS())).GlobalLength();
//    }
//
//    // create the linear system
//    // J(x_i) \Delta x_i = - R(x_i)
//    // create the systemmatrix
//    SetupSystemMatrix();
//
//    // check whether we have a sanely filled tangent matrix
//    if (not systemmatrix_->Filled())
//    {
//      dserror("Effective tangent matrix must be filled here");
//    }
//
//    SetupRHS(*rhs_,firstcall_);
//
//    LinearSolve();
//
//    // reset solver tolerance
//    solver_->ResetTolerance();
//
//    // build residual and incremental norms
//    // for now use for simplicity only L2/Euclidian norm
//    BuildCovergenceNorms();
//
//    // print stuff
//    PrintNewtonIter();
//
//    // increment equilibrium loop index
//    iter_ += 1;
//
//    firstcall_ = false;
//
//  }// end while loop
//
//  // correct iteration counter
//  iter_ -= 1;
//
//  // test whether max iterations was hit
//  if ( (Converged()) and (Comm().MyPID()==0) )
//  {
//  IO::cout << IO::endl;
//    IO::cout << "  Newton Converged! " <<  IO::endl;
//  }
//  else if (iter_ >= itermax_)
//  {
//  IO::cout << IO::endl;
//    IO::cout << " Newton unconverged in "<< iter_ << " iterations " <<  IO::endl;
//  }
}


///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//void FSI::MonolithicXFEM::SetDefaultParameters(const Teuchos::ParameterList& fsidyn,
//                                               Teuchos::ParameterList& list)
//{
//  // Get the top level parameter list
//  Teuchos::ParameterList& nlParams = list;
//
//  nlParams.set<std::string>("Nonlinear Solver", "Line Search Based");
//  //nlParams.set("Preconditioner", "None");
//  //nlParams.set("Norm abs F", fsidyn.get<double>("CONVTOL"));
//
//  nlParams.set("Max Iterations", fsidyn.get<int>("ITEMAX"));
//
//  nlParams.set("Norm abs pres", fsidyn.get<double>("CONVTOL"));
//  nlParams.set("Norm abs vel",  fsidyn.get<double>("CONVTOL"));
//  nlParams.set("Norm abs disp", fsidyn.get<double>("CONVTOL"));
//
//  // sublists
//
//  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
//  Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
//  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
//  //Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");
//
//
//  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
//  dirParams.set<std::string>("Method","User Defined");
////   Teuchos::RCP<NOX::Direction::UserDefinedFactory> newtonfactory = Teuchos::rcp(this,false);
////   dirParams.set("User Defined Direction Factory",newtonfactory);
//
//
//  // status tests are expensive, but instructive
//  //solverOptions.set<std::string>("Status Test Check Type","Minimal");
//  solverOptions.set<std::string>("Status Test Check Type","Complete");
//
//  // be explicit about linear solver parameters
//  lsParams.set<std::string>("Aztec Solver","GMRES");
//  //lsParams.set<std::string>("BiCGStab","GMRES");
//  lsParams.set<std::string>("Orthogonalization","Modified");
//
//  // "r0", "rhs", "norm", "no scaling", "sol"
//  lsParams.set<std::string>("Convergence Test","r0");
//
//  lsParams.set<int>("Size of Krylov Subspace",50);
//  lsParams.set<int>("Max Iterations",1000);
//  lsParams.set<std::string>("Preconditioner","User Defined");
//  lsParams.set<int>("Output Frequency",10);
//  lsParams.set<bool>("Output Solver Details",true);
//
//  // adaptive tolerance settings for linear solver
//  lsParams.set<double>("base tolerance",fsidyn.get<double>("BASETOL")); // relative tolerance
//  lsParams.set<double>("adaptive distance",fsidyn.get<double>("ADAPTIVEDIST")); // adaptive distance
//}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::BuildCovergenceNorms()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::BuildCovergenceNorms()");

  // build map extractors for velocity and pressure dofs
  std::vector<Teuchos::RCP<const Epetra_Map> > fluidvelpres;
  fluidvelpres.push_back(FluidField()->VelocityRowMap());
  fluidvelpres.push_back(FluidField()->PressureRowMap());
  LINALG::MultiMapExtractor fluidvelpresextract(*(FluidField()->DofRowMap()),fluidvelpres);



  //-------------------------------
  // build residual norms
  //-------------------------------

  // build residual norms (the whole vector)
  rhs_->Norm2(&normrhs_);

  // structural Dofs
  Extractor()->ExtractVector(rhs_,0)->Norm2(&normstrrhsL2_);
  Extractor()->ExtractVector(rhs_,0)->NormInf(&normstrrhsInf_);

  // fluid velocity Dofs
  fluidvelpresextract.ExtractVector(Extractor()->ExtractVector(rhs_,1),0)->Norm2(&normflvelrhsL2_);
  fluidvelpresextract.ExtractVector(Extractor()->ExtractVector(rhs_,1),0)->NormInf(&normflvelrhsInf_);

  // fluid pressure Dofs
  fluidvelpresextract.ExtractVector(Extractor()->ExtractVector(rhs_,1),1)->Norm2(&normflpresrhsL2_);
  fluidvelpresextract.ExtractVector(Extractor()->ExtractVector(rhs_,1),1)->NormInf(&normflpresrhsInf_);

  //-------------------------------
  // build solution increment norms
  //-------------------------------

  // build increment norm
  iterinc_->Norm2(&norminc_);

  // structural Dofs
  Extractor()->ExtractVector(iterinc_,0)->Norm2(&normstrincL2_);
  Extractor()->ExtractVector(iterinc_,0)->NormInf(&normstrincInf_);

  // fluid velocity Dofs
  fluidvelpresextract.ExtractVector(Extractor()->ExtractVector(iterinc_,1),0)->Norm2(&normflvelincL2_);
  fluidvelpresextract.ExtractVector(Extractor()->ExtractVector(iterinc_,1),0)->NormInf(&normflvelincInf_);

  // fluid pressure Dofs
  fluidvelpresextract.ExtractVector(Extractor()->ExtractVector(iterinc_,1),1)->Norm2(&normflpresincL2_);
  fluidvelpresextract.ExtractVector(Extractor()->ExtractVector(iterinc_,1),1)->NormInf(&normflpresincInf_);

  //get length of the structural, fluid and ale vector
  ns_  = (*(Extractor()->ExtractVector(rhs_,0))).GlobalLength();                                      //structure
  nf_  = (*(Extractor()->ExtractVector(rhs_,1))).GlobalLength();                                      //fluid
  nfv_ = (*(fluidvelpresextract.ExtractVector(Extractor()->ExtractVector(rhs_,1),0))).GlobalLength(); //fluid velocity
  nfp_ = (*(fluidvelpresextract.ExtractVector(Extractor()->ExtractVector(rhs_,1),1))).GlobalLength(); //fluid pressure
  nall_ = (*rhs_).GlobalLength();                                                                     //all


}








/*----------------------------------------------------------------------*/
/* Recover the Lagrange multiplier at the interface                     */
/*----------------------------------------------------------------------*/
 void FSI::MonolithicXFEM::RecoverLagrangeMultiplier()
 {
//   TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::RecoverLagrangeMultiplier");
//
//   // get time integration parameters of structural time integrator
//   // to enable consistent time integration among the fields
//   const double stiparam = StructureField()->TimIntParam();
//
//  // some often re-used vectors
//  Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::null;     // stores intermediate result of terms (3)-(8)
//  Teuchos::RCP<Epetra_Vector> auxvec = Teuchos::null;     // just for convenience
//  Teuchos::RCP<Epetra_Vector> auxauxvec = Teuchos::null;  // just for convenience
//
//  /* Recovery of Lagrange multiplier lambda^{n+1} is done by the following
//   * condensation expression:
//   *
//   * lambda^{n+1} =
//   *
//   * (1)  - stiparam / (1.-stiparam) * lambda^{n}
//   *
//   * (2)  + 1. / (1.-stiparam) * tmpvec
//   *
//   * with tmpvec =
//   *
//   * (3)    r_{\Gamma}^{S,n+1}
//   *
//   * (4)  + S_{\Gamma I} * \Delta d_{I}^{S,n+1}
//   *
//   * (5)  + tau * S_{\Gamma\Gamma} * \Delta u_{\Gamma}^{F,n+1}
//   *
//   * (6)  + dt * S_{\Gamma\Gamma} * u_{\Gamma}^n]
//   *
//   * Remark on term (6):
//   * Term (6) has to be considered only in the first Newton iteration.
//   * Hence, it will usually not be computed since in general we need more
//   * than one nonlinear iteration until convergence.
//   *
//   * Remarks on all terms:
//   * +  Division by (1.0 - stiparam) will be done in the end
//   *    since this is common to all terms
//   * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
//   * +  neglecting terms (4)-(5) should not alter the results significantly
//   *    since at the end of the time step the solution increments tend to zero.
//   *
//   *                                                 Matthias Mayr (10/2012)
//   */
//
//  // ---------Addressing term (1)
//  lambda_->Update(-stiparam,*lambda_,0.0);
//  // ---------End of term (1)
//
//  // ---------Addressing term (3)
//  Teuchos::RCP<Epetra_Vector> structureresidual = StructureField()->Interface()->ExtractFSICondVector(StructureField()->RHS());
//  structureresidual->Scale(-1.0); // invert sign to obtain residual, not rhs
//  tmpvec = Teuchos::rcp(new Epetra_Vector(*structureresidual));
//  // ---------End of term (3)
//
//  /* You might want to comment out terms (4) to (6) since they tend to
//   * introduce oscillations in the Lagrange multiplier field for certain
//   * material properties of the structure.
//   *                                                    Matthias Mayr 11/2012
//  // ---------Addressing term (4)
//  auxvec = Teuchos::rcp(new Epetra_Vector(sgicur_->RangeMap(),true));
//  sgicur_->Apply(*ddiinc_,*auxvec);
//  tmpvec->Update(1.0,*auxvec,1.0);
//  // ---------End of term (4)
//
//  // ---------Addressing term (5)
//  auxvec = Teuchos::rcp(new Epetra_Vector(sggcur_->RangeMap(),true));
//  sggcur_->Apply(*ddginc_,*auxvec);
//  tmpvec->Update(1.0/timescale,*auxvec,1.0);
//  // ---------End of term (5)
//
//  //---------Addressing term (6)
//  if (firstcall_)
//  {
//    auxvec = Teuchos::rcp(new Epetra_Vector(sggprev_->RangeMap(),true));
//    sggprev_->Apply(*FluidToStruct(FluidField().ExtractInterfaceVeln()),*auxvec);
//    tmpvec->Update(Dt(),*auxvec,1.0);
//  }
//  // ---------End of term (6)
//  *
//  */
//
//  // ---------Addressing term (2)
//  lambda_->Update(1.0,*tmpvec,1.0);
//  // ---------End of term (2)
//
//  // finally, divide by -(1.-stiparam) which is common to all terms
//  lambda_->Scale(1./(1.0-stiparam));
//
//  // Finally, the Lagrange multiplier 'lambda_' is recovered here.
//  // It represents nodal forces acting onto the structure.
//
//
////------ old version changed on 5/12/12  --------------
////   // compute the product S_{\Gamma I} \Delta d_I
////   Teuchos::RCP<Epetra_Vector> sgiddi = LINALG::CreateVector(*StructureField()->Interface()->OtherMap(),true); // store the prodcut 'S_{\GammaI} \Delta d_I^{n+1}' in here
////   (sgicur_->EpetraMatrix())->Multiply(false, *ddiinc_, *sgiddi);
//
////    // compute the product S_{\Gamma\Gamma} \Delta d_\Gamma
////    Teuchos::RCP<Epetra_Vector> sggddg = LINALG::CreateVector(*StructureField()->Interface()->OtherMap(),true); // store the prodcut '\Delta t / 2 * S_{\Gamma\Gamma} \Delta u_\Gamma^{n+1}' in here
////    (sggcur_->EpetraMatrix())->Multiply(false, *ddginc_, *sggddg);
//
////    // Update the Lagrange multiplier:
////    /* \lambda^{n+1} =  - a/b*\lambda^n - f_\Gamma^S
////     *                  - S_{\Gamma I} \Delta d_I - S_{\Gamma\Gamma} \Delta d_\Gamma
////     */
////    lambda_->Update(1.0, *fgcur_, -stiparam);
////    //lambda_->Update(-1.0, *sgiddi, -1.0, *sggddg, 1.0);
////    lambda_->Scale(1/(1.0-stiparam)); //entire Lagrange multiplier it divided by (1.-stiparam)
////

   return;
}

 /*----------------------------------------------------------------------*/
 /*  print Newton-Raphson iteration to screen and error file             */
 /*----------------------------------------------------------------------*/
 void FSI::MonolithicXFEM::PrintNewtonIter()
 {
   // print to standard out
   if ( Comm().MyPID()==0 )
   {
     if (iter_== 1)
       PrintNewtonIterHeader();
     PrintNewtonIterText();
   }
 }
 /*----------------------------------------------------------------------*/
 /* print Newton-Raphson iteration to screen and error file              */
 /*----------------------------------------------------------------------*/
 void FSI::MonolithicXFEM::PrintNewtonIterHeader()
 {
   IO::cout << "CONVTOL: " << tolfres_ << IO::endl;

   // open outstringstream
   //std::ostringstream oss;

   IO::cout << "============================================================================================================================="<< IO::endl;

   // enter converged state etc
   IO::cout << "|nit|";

   // different style due relative or absolute error checking
   // displacement
   switch ( normtypefres_ )
   {
   case INPAR::FSI::convnorm_abs :
     IO::cout <<"            "<< "abs-res-norm  |";
     break;
   case INPAR::FSI::convnorm_rel :
    IO::cout << "str-rs-l2|"  << "flv-rs-l2|" << "flp-rs-l2|" ;
    IO::cout << "str-rs-li|"  << "flv-rs-li|" << "flp-rs-li|" ;
     break;
   case INPAR::FSI::convnorm_mix :
     dserror("not implemented");
     break;
   default:
     dserror("You should not turn up here."); break;
   }

   switch ( normtypeinc_ )
   {
   case INPAR::FSI::convnorm_abs :
     IO::cout <<"                  "<< "abs-inc-norm";
     break;
   case INPAR::FSI::convnorm_rel :
      IO::cout << "str-in-l2|"  << "flv-in-l2|" << "flp-in-l2|" ;
      IO::cout << "str-in-li|"  << "flv-in-li|" << "flp-in-li|" ;
     break;
   case INPAR::FSI::convnorm_mix :
     dserror("not implemented");
     break;
   default:
     dserror("You should not turn up here."); break;
   }

   // add solution time
   IO::cout << IO::endl;
   IO::cout << "============================================================================================================================="<< IO::endl;

 }

 /*---------------------------------------------------------------------*/
 /*  print Newton-Raphson iteration to screen                           */
 /*---------------------------------------------------------------------*/
 void FSI::MonolithicXFEM::PrintNewtonIterText()
 {
   // enter converged state etc
   IO::cout << " " << iter_ << "/" << itermax_;

   // different style due relative or absolute error checking
   // displacement
   switch ( normtypefres_ )
   {
   case INPAR::FSI::convnorm_abs :
     IO::cout << "             " << (normrhs_) << IO::endl;
     break;
   case INPAR::FSI::convnorm_rel :
     IO::cout << "|" << (normstrrhsL2_/ns_)
              << "|" << (normflvelrhsL2_/nfv_)
              << "|" << (normflpresrhsL2_/nfp_)
              << "|" << (normstrrhsInf_)
              << "|" << (normflvelrhsInf_)
              << "|" << (normflpresrhsInf_);
     break;
   case INPAR::FSI::convnorm_mix :
     dserror("not implemented!");
     break;
   default:
     dserror("You should not turn up here."); break;
  }

   switch ( normtypeinc_ )
   {
   case INPAR::FSI::convnorm_abs :
     IO::cout << "             " << (norminc_) << IO::endl;
     break;
   case INPAR::FSI::convnorm_rel :
     IO::cout << "|" << (normstrincL2_/ns_)
              << "|" << (normflvelincL2_/nfv_)
              << "|" << (normflpresincL2_/nfp_)
              << "|" << (normstrincInf_)
              << "|" << (normflvelincInf_)
              << "|" << (normflpresincInf_)
              << "|" << IO::endl;
     break;
   case INPAR::FSI::convnorm_mix :
     dserror("not implemented!");
     break;
   default:
     dserror("You should not turn up here."); break;
   }
 }
