/*!----------------------------------------------------------------------
\file fsi_xfem_monolithic.cpp
\brief Control routine for monolithic FSI (XFSI) solved via a classical Newton scheme
       taking into account changing fluid dofsets

\level 2

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*----------------------------------------------------------------------*/

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_poro_wrapper.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_poroelast/poroelast_monolithic.H"
#include "../drt_adapter/ad_ale_fpsi.H"
#include "../drt_fluid_xfluid/xfluid.H"

#include "../drt_fsi/fsi_debugwriter.H"

#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_solver.H"
#include "../drt_inpar/inpar_fsi.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_colors.H"

#include "../drt_structure/stru_aux.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_io/io_pstream.H"

#include "XFScoupling_manager.H"
#include "XFAcoupling_manager.H"
#include "XFPcoupling_manager.H"

#include  "../drt_xfem/xfem_condition_manager.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>

#include "fsi_xfem_monolithic.H"
/*----------------------------------------------------------------------*/
// constructor
/*----------------------------------------------------------------------*/
FSI::MonolithicXFEM::MonolithicXFEM(const Epetra_Comm& comm,
                                    const Teuchos::ParameterList& timeparams,
                                    const ADAPTER::FieldWrapper::Fieldtype type)
    : AlgorithmXFEM(comm, timeparams, type),
    fsidyn_(DRT::Problem::Instance()->FSIDynamicParams()),
    fsimono_(fsidyn_.sublist("MONOLITHIC SOLVER")),
    solveradapttol_(true),
    solveradaptolbetter_(fsimono_.get<double>("ADAPTIVEDIST")), // adaptive distance
    merge_fsi_blockmatrix_(false),
    scaling_infnorm_((bool)DRT::INPUT::IntegralValue<int>(fsimono_,"INFNORMSCALING")),
//  erroraction_(erroraction_stop),
    log_(Teuchos::null),
//  logada_(Teuchos::null)
    /// tolerance and for linear solver
    tolrhs_(fsimono_.get<double>("BASETOL")),     // absolute tolerance for full residual for adapting the linear solver
    /// iteration counter
    iter_         (0),
    iter_outer_   (0),
    itermin_      (0),
    itermax_      (fsimono_.get<int>("ITEMAX")),
    itermax_outer_(5),
    /// Convergence criterion and convergence tolerances for Newton scheme
    normtypeinc_  (DRT::INPUT::IntegralValue<INPAR::FSI::ConvNorm>(fsimono_,"NORM_INC")),
    normtypefres_ (DRT::INPUT::IntegralValue<INPAR::FSI::ConvNorm>(fsimono_,"NORM_RESF")),
    combincfres_  (DRT::INPUT::IntegralValue<INPAR::FSI::BinaryOp>(fsimono_,"NORMCOMBI_RESFINC")),
    tolinc_       (fsimono_.get<double>("CONVTOL")),
    tolfres_      (fsimono_.get<double>("CONVTOL")),
    /// set tolerances for nonlinear solver
    /// tolerances for structural displacements
    TOL_DIS_RES_L2_  (fsimono_.get<double>("TOL_DIS_RES_L2")),
    TOL_DIS_RES_INF_ (fsimono_.get<double>("TOL_DIS_RES_INF")),
    TOL_DIS_INC_L2_  (fsimono_.get<double>("TOL_DIS_INC_L2")),
    TOL_DIS_INC_INF_ (fsimono_.get<double>("TOL_DIS_INC_INF")),
    /// tolerances for fluid pressure
    TOL_PRE_RES_L2_  (fsimono_.get<double>("TOL_PRE_RES_L2")),
    TOL_PRE_RES_INF_ (fsimono_.get<double>("TOL_PRE_RES_INF")),
    TOL_PRE_INC_L2_  (fsimono_.get<double>("TOL_PRE_INC_L2")),
    TOL_PRE_INC_INF_ (fsimono_.get<double>("TOL_PRE_INC_INF")),
    /// tolerances for fluid velocity
    TOL_VEL_RES_L2_  (fsimono_.get<double>("TOL_VEL_RES_L2")),
    TOL_VEL_RES_INF_ (fsimono_.get<double>("TOL_VEL_RES_INF")),
    TOL_VEL_INC_L2_  (fsimono_.get<double>("TOL_VEL_INC_L2")),
    TOL_VEL_INC_INF_ (fsimono_.get<double>("TOL_VEL_INC_INF"))
{
  //TODO set some of these flags via the input file


//  const Teuchos::ParameterList& xdyn       = DRT::Problem::Instance()->XFEMGeneralParams();
//  const Teuchos::ParameterList& xfluiddyn  = DRT::Problem::Instance()->XFluidDynamicParams();

  //-------------------------------------------------------------------------
  // enable debugging
  //-------------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<int>(fsidyn_,"DEBUGOUTPUT")==1)
  {
    // debug writer for structure field
    sdbg_ = Teuchos::rcp(new UTILS::DebugWriter(StructurePoro()->Discretization()));
    // debug writer for fluid field
    fdbg_ = Teuchos::rcp(new UTILS::DebugWriter(FluidField()->Discretization()));
  }
  //-------------------------------------------------------------------------
  // write files
  //-------------------------------------------------------------------------

  // write iterations-file
  std::string fileiter = DRT::Problem::Instance()->OutputControlFile()->FileName();
  fileiter.append(".iteration");
  log_ = Teuchos::rcp(new std::ofstream(fileiter.c_str()));

  // write energy-file
  if (DRT::INPUT::IntegralValue<int>(fsidyn_.sublist("MONOLITHIC SOLVER"), "ENERGYFILE") == 1)
  {
    dserror("writing energy not supported yet");
//  TODO
//    std::string fileiter2 = DRT::Problem::Instance()->OutputControlFile()->FileName();
//    fileiter2.append(".fsienergy");
//    logenergy_ = Teuchos::rcp(new std::ofstream(fileiter2.c_str()));
  }


  //-------------------------------------------------------------------------
  // time step size adaptivity
  //-------------------------------------------------------------------------
  const bool timeadapton = DRT::INPUT::IntegralValue<bool>(fsidyn_.sublist("TIMEADAPTIVITY"),"TIMEADAPTON");

  if (timeadapton)
  {
    dserror("FSI - TimeIntAdaptivity not supported for XFEM yet");
    //InitTimIntAda(fsidyn);
  }


  //-------------------------------------------------------------------------
  // Create direct or iterative solver for XFSI system
  //-------------------------------------------------------------------------
  CreateLinearSolver();


  //-------------------------------------------------------------------------
  // validate parameters for monolithic approach
  //-------------------------------------------------------------------------
  ValidateParameters();

  //-------------------------------------------------------------------------
  // Setup Coupling Objects
  //-------------------------------------------------------------------------
  SetupCouplingObjects();

  //build ale system matrix in splitted system
  if (HaveAle())
    AleField()->CreateSystemMatrix(AleField()->Interface());


  //-------------------------------------------------------------------------
  // Finish standard FluidField()->Init()!
  // REMARK: We don't want to do this at the beginning, to be able to use std ADAPTER::Coupling for FA-Coupling
  //-------------------------------------------------------------------------
  FluidField()->CreateInitialState();

  //Todo: move that somewhere else
  {
    // set initial field by given function
    // we do this here, since we have direct access to all necessary parameters
    const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
    INPAR::FLUID::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn,"INITIALFIELD");
    if(initfield != INPAR::FLUID::initfield_zero_field)
    {
      int startfuncno = fdyn.get<int>("STARTFUNCNO");
      if (initfield != INPAR::FLUID::initfield_field_by_function and
          initfield != INPAR::FLUID::initfield_disturbed_field_from_function)
      {
        startfuncno=-1;
      }
      FluidField()->SetInitialFlowField(initfield,startfuncno);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | SetupCouplingObjects                                      ager 06/16 |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupCouplingObjects()
{
  int coup_idx = 0;
  std::vector<int> idx;

  //Just Add coupling object in case there is an FSI Interface
  if (FluidField()->GetConditionManager()->GetMeshCoupling("XFEMSurfFSIMono") != Teuchos::null)
  {
    idx.push_back(structp_block_);
    idx.push_back(fluid_block_);
    coup_man_[coup_idx] = Teuchos::rcp( new XFEM::XFSCoupling_Manager(FluidField()->GetConditionManager(),
        StructurePoro()->StructureField(), FluidField(),idx));
  }

  if (HaveAle())
  {
    ++coup_idx;
    idx.clear();
    idx.push_back(fluid_block_);
    idx.push_back(ale_i_block_);
    coup_man_[coup_idx] = Teuchos::rcp (new XFEM::XFACoupling_Manager(FluidField(),AleField(), idx));
  }

  if (StructurePoro()->isPoro())
  {
    //Just Add coupling object in case there is an FPI Interface
    if (FluidField()->GetConditionManager()->GetMeshCoupling("XFEMSurfFPIMono_ps_ps") != Teuchos::null)
    {
      ++coup_idx;
      idx.clear();
      idx.push_back(structp_block_);
      idx.push_back(fluid_block_);
      coup_man_[coup_idx] = Teuchos::rcp( new XFEM::XFPCoupling_Manager(FluidField()->GetConditionManager(), StructurePoro()->PoroField(), FluidField(),idx));
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | validate the input parameter combinations               schott 07/14 |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::ValidateParameters()
{
  // check for reasonable input parameter combinations!

  //Check for the timestepsize
  if (fabs(FluidField()->Dt() - StructurePoro()->StructureField()->Dt()) > 1e-16)
    dserror("ValidateParameters(): Timestep of fluid and structure not equal (%d != %d)!", FluidField()->Dt(), StructurePoro()->StructureField()->Dt());
  if (HaveAle())
    if (fabs(FluidField()->Dt() - AleField()->Dt()) > 1e-16)
      dserror("ValidateParameters(): Timestep of fluid and ale not equal (%d != %d)!", FluidField()->Dt(), AleField()->Dt());
  if (StructurePoro()->isPoro())
    if (fabs(FluidField()->Dt() - StructurePoro()->PoroField()->Dt()) > 1e-16)
      dserror("ValidateParameters(): Timestep of fluid and poro not equal (%d != %d)!", FluidField()->Dt(), StructurePoro()->PoroField()->Dt());

  //TODO
  // REMARK: be aware of using const Dis predictor!
  //         This results in zero disp_incr and u^n+1 = -u^n for second order disp_to_vel interface conversion
  //                                        and u^n+1 = 0    for first order disp_to_vel interface conversion
return;
}

/*----------------------------------------------------------------------*
 | setup of the monolithic XFSI system,                    schott 08/14 |
 | setup a new combined block row map and a new block matrix            |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupSystem()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::SetupSystem()");

  /*----------------------------------------------------------------------
   Create a combined map for Structure/Fluid-DOFs all in one!
   ----------------------------------------------------------------------*/

  CreateCombinedDofRowMap();


  /*----------------------------------------------------------------------
    Initialise XFSI-systemmatrix_
   ----------------------------------------------------------------------*/
  CreateSystemMatrix();

}

/*----------------------------------------------------------------------*
 | setup of the monolithic XFSI system,                    schott 08/14 |
 | setup a new combined block row map and a new block matrix            |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::CreateSystemMatrix()
{
  if(Comm().MyPID() == 0)
    std::cout << "Create a new global systemmatrix (BlockSparseMatrix)" << std::endl;

  //TODO: check the savegraph option and explicit Dirichlet flag for the matrix!
  //TODO: check if it is okay to use a BlockSparseMatrix without the FE-flag for the fluid-submatrix and the coupling submatrices?
  //TODO: do we add a already communicated (completed) fluid matrix?!
  //TODO: check the number of non-zeros predicted
  /*----------------------------------------------------------------------*/

  if(systemmatrix_.strong_count() > 1)
    dserror("deleting systemmatrix does not work properly, the number of RCPs pointing to it is %i", systemmatrix_.strong_count());

  // do not want to have two sysmats in memory at the same time
  if(systemmatrix_ != Teuchos::null)
  {
    if(Comm().MyPID() == 0)
      std::cout << "Delete the global systemmatrix (BlockSparseMatrix)" << std::endl;
    systemmatrix_ = Teuchos::null;
  }

  systemmatrix_
  = Teuchos::rcp(
      new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          Extractor(),
          Extractor(),
          0,
          false, // explicit dirichlet, do not change the graph and do not create a new matrix when applying Dirichlet values
          false  // savegraph (used when submatrices will be reset), we create new fluid sysmats anyway
      )
  );
}


/*----------------------------------------------------------------------*
 * setup composed system matrix from field solvers, complete the global system matrix
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupSystemMatrix()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::SetupSystemMatrix");

  // reset the block system matrix
  // note: Zero() is not sufficient for the coupling blocks,
  // as the couplings between fluid and structure can change (structure moves between iterations)
  // while the fluid dofsets remain unchanged
//  systemmatrix_->Reset();

  /*----------------------------------------------------------------------*/
  // extract Jacobian matrices and put them into composite system
  Teuchos::RCP<LINALG::SparseMatrix> s = StructurePoro()->SystemMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField()->SystemMatrix();

#ifdef DEBUG
  if (s==Teuchos::null)
    dserror("expect structure matrix");
  if (f==Teuchos::null)
    dserror("expect fluid matrix");
#endif

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // get time integration parameters of structure and fluid time integrators
  // as well as of the FSI time-integration
  // to enable consistent time integration among the fields
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/

  /*----------------------------------------------------------------------*/
  // scaling factors for fluid terms/blocks
  // inverse of the weighting of the quantities w.r.t the new time step
  const double scaling_F   = FluidField()->ResidualScaling(); // 1/(theta * dt) = 1/weight^F_np

  /*----------------------------------------------------------------------*/
  // this is the interpolation weight for quantities from last time step
  // alpha_f for genalpha and (1-theta) for OST (weighting of the old time step n for displacements)
  const double stiparam = StructurePoro()->TimIntParam();    // (1-theta) for OST and alpha_f for Genalpha
  // scale factor for the structure system matrix w.r.t the new time step
  const double scaling_S = 1.0/(1.0-stiparam);  // 1/(1-alpha_F) or 1/theta = 1/weight^S_np


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/


  /*----------------------------------------------------------------------*/
  // Structure diagonal block (structural system matrix)
  /*----------------------------------------------------------------------*/

  // Uncomplete structure matrix to be able to deal with slightly defective interface meshes.
  //
  // The additional coupling block C_ss can contain additional non-zero entries,
  // e.g. from DBCs which are already applied to s in the structural evaluate, however, not
  // to the coupling block C_ss yet
  s->UnComplete();

  // NOTE: UnComplete creates a new Matrix and a new matrix graph as well which is not allocated with staticprofile
  // Therefore, the savegraph = true option set in structural timint has no effect, as a new graph is created
  // whenever UnComplete is called
  // then, due to memory fragmentation the evaluate time for the structure can vary a lot!
  // UPDATE: actually after the next complete with store the graph

  // scale the structure system matrix
  s->Scale(scaling_S);

  // assign the structure sysmat diagonal block
  systemmatrix_->Assign(0,0,LINALG::View,*s);

  /*----------------------------------------------------------------------*/
  // Fluid diagonal block
  /*----------------------------------------------------------------------*/

  // scale the fluid diagonal block
  f->Scale(scaling_F); //<  1/(theta_f*dt) = 1/weight(t^f_np)

  // assign the fluid diagonal block
  systemmatrix_->Assign(1,1,LINALG::View,*f);

  //Add Coupling Sysmat
  for (std::map<int, Teuchos::RCP<XFEM::Coupling_Manager> >::iterator coupit = coup_man_.begin(); coupit != coup_man_.end(); ++coupit)
    coupit->second->AddCouplingMatrix(*systemmatrix_,scaling_F);

  /*----------------------------------------------------------------------*/
  // Complete the global system matrix
  /*----------------------------------------------------------------------*/

  // done. make sure all blocks are filled.
  systemmatrix_->Complete();
}


/*----------------------------------------------------------------------*/
// setup composed right hand side from field solvers
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupRHS()
{

  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::SetupRHS");

  // We want to add into a zero vector
  rhs_->PutScalar(0.0);

  // contributions of single field residuals
  SetupRHSResidual(*rhs_);

  //Add Coupling RHS
  const double scaling_F   = FluidField()->ResidualScaling();
  for (std::map<int, Teuchos::RCP<XFEM::Coupling_Manager> >::iterator coupit = coup_man_.begin(); coupit != coup_man_.end(); ++coupit)
    coupit->second->AddCouplingRHS(rhs_,Extractor(),scaling_F);
}

/*----------------------------------------------------------------------*/
// setup RHS contributions based on single field residuals
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupRHSResidual(Epetra_Vector &f)
{
  /*----------------------------------------------------------------------*/
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  /*----------------------------------------------------------------------*/

  /*----------------------------------------------------------------------*/
  // scaling factors for fluid terms/blocks
  // inverse of the weighting of the quantities w.r.t the new time step
  const double scaling_F   = FluidField()->ResidualScaling(); // 1/(theta * dt) = 1/weight^F_np

  /*----------------------------------------------------------------------*/
  // this is the interpolation weight for quantities from last time step
  // alpha_f for genalpha and (1-theta) for OST (weighting of the old time step n for displacements)
  const double stiparam = StructurePoro()->StructureField()->TimIntParam();    // (1-theta) for OST and alpha_f for Genalpha

  // scale factor for the structure system matrix w.r.t the new time step
  const double scaling_S = 1.0/(1.0-stiparam);  // 1/(1-alpha_F) = 1/weight^S_np


  /*----------------------------------------------------------------------*/
  // get single field residuals
  Teuchos::RCP<Epetra_Vector> sv = Teuchos::rcp(new Epetra_Vector(*StructurePoro()->RHS()));
  Teuchos::RCP<Epetra_Vector> fv = Teuchos::rcp(new Epetra_Vector(*FluidField()->RHS()));

  // scale the structural rhs
  sv->Scale(scaling_S);

  // scale the fluid rhs
  fv->Scale(scaling_F); // scale with FluidField()->ResidualScaling()

  // put the single field residuals together
  CombineFieldVectors(f,sv,fv);

return;
}

/*----------------------------------------------------------------------*/
// apply Dirichlet boundary conditions to XFSI system
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::ApplyDBC()
{

  // note, the structural evaluate applies DBC already to the structure block,
  // however, the additional C_ss coupling block is added and we have to apply DBC again to the sum
  // no DBC are applied for in the fluid-evaluate for the ff-diagonal block or the fluid coupling blocks
  // therefore we apply BCS to the whole system.
  // Just the internal DBCs via the structural evaluate (PrepareSystemForNewtonSolve()) are applied twice

#if(1)
  // apply combined Dirichlet to whole XFSI system
  LINALG::ApplyDirichlettoSystem(
      systemmatrix_,
      iterinc_,
      rhs_,
      Teuchos::null, //possible trafo?!
      zeros_,
      *CombinedDBCMap()
  );

#else
  // TODO: transformations for locsys
  // be aware of the rhs_Cs-coupling block when using localsys

  // alternatively, when trafos for the structure have to be used, we apply DBC for each field separately
  ((*systemmatrix_)(0,0)).ApplyDirichlet(*StructureField()->GetDBCMapExtractor()->CondMap(),true);  // set diagonal 1.0 for SS-Block
  ((*systemmatrix_)(0,1)).ApplyDirichlet(*StructureField()->GetDBCMapExtractor()->CondMap(),false); // do not set diagonal 1.0 for SF-Block
  ((*systemmatrix_)(1,0)).ApplyDirichlet(*FluidField()->GetDBCMapExtractor()->CondMap(),false);     // do not set diagonal 1.0 for FS-Block
  ((*systemmatrix_)(1,1)).ApplyDirichlet(*FluidField()->GetDBCMapExtractor()->CondMap(),true);      // set diagonal 1.0 for FF-Block

  LINALG::ApplyDirichlettoSystem(rhs_, zeros_, *CombinedDBCMap());

#endif

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  //TODO: what to do with this function???
//  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::InitialGuess");
//
//  SetupVector(*ig,
//              StructureField()->InitialGuess(),
//              FluidField().InitialGuess(),
//              AleField().InitialGuess(),
//              0.0);
}


/*----------------------------------------------------------------------*/
// Create the combined DOF row map for the FSI problem;
// row maps of structure and xfluid to an global FSI DOF row map
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::CreateCombinedDofRowMap()
{
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;

  // Append the structural DOF map
  vecSpaces.push_back(StructurePoro()->DofRowMap());

  // Append the background fluid DOF map
  vecSpaces.push_back(FluidField()->DofRowMap());

  // solid maps empty??
  if (vecSpaces[0]->NumGlobalElements() == 0)
    dserror("No solid equations. Panic.");

  // fluid maps empty??
  if (vecSpaces[1]->NumGlobalElements() == 0)
    dserror("No fluid equations. Panic.");

  // Append the background fluid DOF map
  if (HaveAle())
  {
   vecSpaces.push_back(AleField()->Interface()->OtherMap());

   // ale maps empty??
   if (vecSpaces[2]->NumGlobalElements() == 0)
     dserror("No ale equations. Panic.");
  }

  // The vector is complete, now fill the system's global block row map
  // with the maps previously set together!
  SetDofRowMaps(vecSpaces);

  return;
}




/*----------------------------------------------------------------------*/
// set full monolithic dof row map
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps)
{

  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  blockrowdofmap_.Setup(*fullmap,maps);
}



/*----------------------------------------------------------------------*
 | Put two field vectors together to a monolithic vector
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::CombineFieldVectors(
      Epetra_Vector &v,                       ///< composed vector containing all field vectors
      Teuchos::RCP<const Epetra_Vector> sv,   ///< structural DOFs
      Teuchos::RCP<const Epetra_Vector> fv    ///< fluid DOFs
)
{
  Extractor().AddVector(*sv, 0, v);
  Extractor().AddVector(*fv, 1, v);
}



/*----------------------------------------------------------------------
 |                                                        schott 08/14 |
 |   extract the two field vectors from a given composed vector        |
 ----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector>  x,
                                              Teuchos::RCP<const Epetra_Vector>& sx,
                                              Teuchos::RCP<const Epetra_Vector>& fx,
                                              Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::ExtractFieldVectors");

  /*----------------------------------------------------------------------*/
  //Process structure unknowns
  /*----------------------------------------------------------------------*/
  // Extract whole structure field vector
  sx = Extractor().ExtractVector(x,structp_block_);

  /*----------------------------------------------------------------------*/
  //Process fluid unknowns
  /*----------------------------------------------------------------------*/
  // Extract vector of fluid unknowns from x
  fx = Extractor().ExtractVector(x,fluid_block_);

  // Extract vector of ale unknowns from x
  if (HaveAle())
    ax = Extractor().ExtractVector(x,ale_i_block_);
}



/*----------------------------------------------------------------------*
 | time loop of the monolithic system                      schott 08/14 |
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

    //TODO: check this function
    // TODO: erst update und dann prepare output? oder anders rum?
    // calculate stresses, strains, energies
    PrepareOutput();

    // update all single field solvers
    Update();

    // write output to screen and files
    Output();

  }  // NotFinished
}  // TimeLoop()


/*----------------------------------------------------------------------*
 | prepare the time step for fluid and structure           schott 08/14 |
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
  StructurePoro()->PrepareTimeStep();

  //--------------------------------------------
  // Fluid PrepareTimeStep
  //--------------------------------------------
  // * set time integrator
  // * set time parameters
  // * apply fluid predictor (before the cut is performed the first time in the new time step, see first evaluate call)
  // * DBCs will be applied to predicted solution in a PrepareNonlinearSolve-call within the evaluate
  //   after velnp has been mapped to the new interface position. DBCs will be applied again for each iteration
  //   as the CUT is performed for each increment and therefore the DBCs have to be set again
  FluidField()->PrepareTimeStep();

  if (HaveAle())
    AleField()->PrepareTimeStep();

}


/*----------------------------------------------------------------------*
 | outer iteration loop, restarts inner Newton-Raphson iterations       |
 | when fluid dofsets changes                              schott 08/14 |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Solve()
{
  // initialize outer loop iteration index which allows for restarts of the Newton scheme in case of changing fluid maps
  iter_outer_ = 1;

  // reset the single step-increments for structural and fluid field
  sx_sum_ = Teuchos::null;
  fx_sum_ = Teuchos::null;
  ax_sum_ = Teuchos::null;


  // We want to make sure, that the outer loop is entered at least once!
  // We exit the outer loop if either the inner Newton loop is converged (checked in Converged() method)
  // OR the maximum number of outer iterations is exceeded!
  while ( (iter_outer_ == 1) or ( iter_outer_ <= itermax_outer_) )
  {
    //--------------------------------------------------------
    // call the inner Newton loop and check for convergence
    //--------------------------------------------------------
    if(Newton()) // stop since the main inner Newton loop converged
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
        IO::cout << "---------------------------------------- Restart Newton-Raphson - DOF-sets changed -----------------------------------------" << IO::endl;
      }
    }

    iter_outer_++;
  }

  if(iter_outer_ > itermax_outer_)
  {
    if( Comm().MyPID() == 0)
    {
      IO::cout << "-------------------------- Maximum number of restarts reached - Fluid DOF-sets have changed too often ----------------------" << IO::endl;
    }
  }
}


/*----------------------------------------------------------------------*
 | recover Lagrange multiplier (structural forces) needed for rhs in    |
 | next time step and update single fields                 schott 08/14 |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Update()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::Update");

  const double scaling_F   = FluidField()->ResidualScaling(); // 1/(theta * dt) = 1/weight^F_np
  for (std::map<int, Teuchos::RCP<XFEM::Coupling_Manager> >::iterator coupit = coup_man_.begin(); coupit != coup_man_.end(); ++coupit)
    coupit->second->Update(scaling_F);

  // update the single fields
  StructurePoro()->Update();
  FluidField()->Update();
  if (HaveAle()) AleField()->Update();
}



/*----------------------------------------------------------------------*
 | write output                                            schott 08/14 |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Output()
{
  //--------------------------------
  // output for structural field
  //--------------------------------
  StructurePoro()->Output();

  //--------------------------------
  // output for Lagrange multiplier field (ie forces onto the structure, Robin-type forces
  // consisting of fluid forces and the Nitsche penalty term contribution)
  //--------------------------------
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  const int uprestart = fsidyn.get<int>("RESTARTEVRY");
  const int upres = fsidyn.get<int>("RESULTSEVRY");
  if((uprestart != 0 && FluidField()->Step() % uprestart == 0) || FluidField()->Step() % upres == 0) //Fluid desides about restart, write output
  {
      for (std::map<int, Teuchos::RCP<XFEM::Coupling_Manager> >::iterator coupit = coup_man_.begin(); coupit != coup_man_.end(); ++coupit)
      coupit->second->Output(*StructurePoro()->StructureField()->DiscWriter());
  }

  //--------------------------------
  // output for fluid field - writes the whole GMSH output if switched on
  //--------------------------------
  FluidField()->Output();
  FluidField()->LiftDrag();


  //--------------------------------
  if (StructurePoro()->GetConstraintManager()->HaveMonitor())
  {
    StructurePoro()->GetConstraintManager()->ComputeMonitorValues( StructurePoro()->Dispnp());
    if(Comm().MyPID() == 0)
       StructurePoro()->GetConstraintManager()->PrintMonitorValues();
  }

    if (HaveAle())
      AleField()->Output();

  return;
}

/*----------------------------------------------------------------------*
 | inner iteration loop (Newton-Raphson scheme)            schott 08/14 |
 | return "true" if converged or                                        |
 | "false" if unconverged or in case of changing fluid dof maps         |
 *----------------------------------------------------------------------*/
bool FSI::MonolithicXFEM::Newton()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::Newton");


  //--------------------------------------------------------
  // Perform Newton-Raphson iterations as long as the size of the global system does not change between iterations.
  // * Due to the different structural interface positions the fluid dofsets can change:
  //   -> in case that the number of dofsets per node changes or a simple copying of dofs between std and multiple ghost dofsets
  //      of a node is not possible anymore we have to restart the Newton-scheme
  // * During the Newton iterations the sorting of fluid-dofsets can change. However, for the update of vectors
  //   in a Newton scheme we have to use vectors that do not change ordering of dofs. Therefore permutation of the
  //   increments and other vectors before and after solving the linear systems can be necessary
  // * Then, the convergence of the Newton scheme is observed based on vectors whose dofsets do not change during this restart
  //   call.
  // * When restarting the Newton, we have to ensure that the step-increments for structure and fluid are initialized
  //   properly, see below.
  // * The fluid and structure evaluate routines expect the full step-increment w.r.t old solution t^n, therefore we have
  //   to sum up the Newton-increments to full step-increment
  // * Applying Dirichlet values to rhs and to system has to be done
  //   for the global system after Evaluate routines for the single fields have been called and have summed up to
  //   the global residual
  // * in xfluid we set Dirichlet values into velnp, therefore the stepinc used in evaluate adds zeros for Dirichlet entries
  //   not to modify the already set Dirichlet values
  //--------------------------------------------------------


  /*----------------------------------------------------------------------*/
  // Initialization
  /*----------------------------------------------------------------------*/
  //Iteration counter
  iter_ = 1;

  permutation_map_.clear();
  permutation_.clear();


  /*----------------------------------------------------------------------*/
  // Create a new solver, important in particular when fluid null space changes
  /*----------------------------------------------------------------------*/
  CreateLinearSolver();

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // Newton-Raphson iteration with unchanging, however, permuting dofsets
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // We want to make sure, that the loop is entered at least once!
  // We exit the loop if either the convergence criteria are met (checked in
  // Converged() method) OR the maximum number of inner iterations is exceeded!
  while ( (iter_ == 1) or ( (not Converged()) and (iter_ <= itermax_) ) )
  {
//    std::cout << "Evaluate-Call " << "iter_ " << iter_ << "/" << itermax_ << std::endl;

    /*----------------------------------------------------------------------*/
    //Evaluate()- call
    // * calls evaluate methods for single fields
    // * assembles single field rhs and system-matrices and fluid-structure coupling blocks
    // * check if dofsets between last two Newton iterations drastically changed or maybe simply permuted
    /*----------------------------------------------------------------------*/

    if(systemmatrix_.strong_count() > 1)
      dserror("deleting block sparse matrix does not work properly, the number of RCPs pointing to it is %i", systemmatrix_.strong_count());

    // reduce counter in the rcp-pointers pointing to underlying epetra matrix objects which actually hold large chunks of memory
    // this ensures that the single field matrices can be really deleted (memory can be freed)
    // before we can create a new state class in fluid's Evaluate
    // NOTE: the blocksparsematrix' sparse matrices hold strong RCP's to the single-fields EpetraMatrix objects
    // NOTE: fluid's Evaluate will create a new LINALG::SparseMatrix and coupling matrices anyway
    systemmatrix_ = Teuchos::null;
    //TODO: can we delete the solver here? this is done in solver_->Reset after solving the last system
//    solver_ = Teuchos::null;

    const bool changed_fluid_dofsets = Evaluate();


    //-------------------
    // store structure step-increment
    //-------------------
    //TODO: Zwischenspeicherung von Structure-step increment
//    w.r.t. first interface position of the respective Newton-restart call k: (U,P)^n_(0,k)
        // set at the beginning of the first outer iteration using the structural PrepareTimeStep-call
//    if(iter_ == 1) sx_sum_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->DofRowMap()));
//
//    sx_sum_->Update(1.0, *Extractor().ExtractVector(x_sum_,0), 0.0);
//
//    //create a copy from the structural part of x_sum_
//    Extractor().ExtractVector(x_sum_,0,sx_sum_);
//    sx_sum_ = sx; // Teuchos::rcp(new Epetra_Vector(*sx));

    //-------------------
    // store fluid step-increment
    //-------------------
    // this step-increment is based on the old solution from t^n (veln) mapped to the interface position of the current
    // Newton iteration (via XFEM-timeintegration) and based on the permutation of dofsets of the current Newton iteration
    // in contrast to the fluid-block of the global x_sum_-step increment (see below)
    // note: current iteration velnp (velnp_ip) and veln have been already mapped to the current interface position
    //       via the Evaluate-call above
    // safe velnp_i+1 - veln in fx_sum_ as start for the Newton as Fluid-Evaluate expects the full step-increment

    fx_sum_ = Teuchos::rcp(new Epetra_Vector(*FluidField()->DofRowMap()));
    fx_sum_->Update(1.0,*FluidField()->Velnp(), -1.0, *FluidField()->Veln(), 0.0);


    /*----------------------------------------------------------------------*/
    // Perform at least one solve after the first Evaluate-call.
    // For further iterations, decide if a Newton-restart is required based on the first computed Newton increment
    /*----------------------------------------------------------------------*/
    if(iter_ > 1)
    {
      if(changed_fluid_dofsets) return false;

      SetupSystem(); // set new blockdofrowmap and create a new global block sparse matrix
    }
    else if(iter_ == 1) // the first run
    {
      //-------------------
      // initialize a new system after the first evaluate-run
      //-------------------
      // note: after the first fluid-evaluate run, the size of the system is now known and new global vectors can be created
      //       until a new restart has to be performed, the dofsets are assumed not to change, however, the std-ghost dofsets
      //       of a node can permute as a a priori sorting of sets is not possible when the interface position changes

      //-------------------
      // setup a new system since the size of the system has changed when NewtonFull is called
      //-------------------
      SetupSystem(); // set new blockdofrowmap and create a new global block sparse matrix


      //-------------------
      // Create a new Epetra_Vector with a given row map and initialize it!
      // DofRowMap() contains all DOFS for the monolithic system w.r.t the respective interface position and returns the
      // DofRowMap of the blockrowdofmap_, a MapExtractor object for a DOF map, split into blocks.
      //-------------------

      //-------------------
      // new global vectors based on UNCHANGING/NON-PERMUTING dofsets
      //-------------------

      // Global increment sum vector (= Global step-increment vector)
      // That's the total increment (structure + fluid) w.r.t the old time step t^n.
      // When using a structural predictor or when restarting the Newton, we have to initialize this vector
      // such that it contains the increment w.r.t the old time step t^n and not w.r.t the restarted solution
      // NOTE:
      //   - the structural dofs do not change, x_sum(0) (structure step-increment) is given w.r.t predictor-solution
      //     from the beginning of the new time-step
      //   - the fluid dofs can permute between Newton iterations and dofsets can completely change such that a Newton restart
      //     is required. x_sum(1) (fluid step-increment) is given w.r.t the old solution at t^n (veln)
      //     mapped/permuted to the interface position of first evaluate-call after (re-)starting the Newton scheme
      //     (in contrast to fx_sum_!) and so the global dofset of this vector remains unchanged/non-permuted
      //     until the next restart has to be performed. In order to update this vector and to use it for further evaluate-calls
      //     permutations backward/forward of the fluid block have to be applied
      x_sum_ = LINALG::CreateVector(*DofRowMap(), true);


      //-------------------
      // new global vectors based on UNCHANGING/NON-PERMUTING dofsets before and after the solve,
      // however, PERMUTING dofsets during the solve itself
      //-------------------

      // Global solution vector for linear solve = iteration increment Delta x = x^n+1_i+1 - x^n+1_i
      // based on non-changing, however, permuting dofsets during the Newton.
      // During the solve we use permuting fluid vectors, depending on the std/ghost-dofset order of each iteration.
      // After the solve we directly have to permute the fluid block Newton increment backwards to the reference ordering
      // (see x_sum) to update the step-increment (xsum) which does not change between the Newton iterations
      // as long as no Newton restart is necessary (see LinearSolve())
      // note: during the solve the increment vector can have permuted dofsets,
      //       directly after the LinearSolve the vector is permuted backwards to the initial ordering of creation
      iterinc_ = LINALG::CreateVector(*DofRowMap(), true);

      // Global residual vector, unchanged/permuting dofsets (the same as for the iterinc_ vector)
      // note: for the assembly and during the solve the residual vector can have permuted dofsets,
      //       directly after the LinearSolve the vector is permuted backwards to the initial ordering of creation
      rhs_ = LINALG::CreateVector(*DofRowMap(), true);

      // Global zero vector for DBCs, unchanged/permuting dofsets (the same as for the iterinc_ vector)
      // note: this vector is just used during the solve and is NOT permuted backwards as iterinc or rhs
      zeros_ = LINALG::CreateVector(*DofRowMap(), true);


    }
    else dserror("the Newton iteration index is assumed to be >= 0");


    /*----------------------------------------------------------------------*/
    // initialize the structural and fluid part of the global step-increment
    /*----------------------------------------------------------------------*/
    if(iter_ == 1)
    {
      //-------------------
      // initialize the structural part of the global step-increment
      //-------------------
      // it is based on the structural predictor-solution D^(n+1)_(pred,k=0)
      // and has been set at the beginning of the first outer iteration using the structural PrepareTimeStep-call

      // if not a new time-step, take the step-increment which has been summed up so far
      if(sx_sum_ != Teuchos::null) Extractor().AddVector(sx_sum_, structp_block_, x_sum_);

      //-------------------
      // initialize the fluid part of the global step-increment
      //-------------------
      // note: during the Newton fx_sum_ and x_sum_(1) can have different dofset orderings.
      //       However, in the first iteration, the two dofsets are equal and permutation is not necessary at the beginning.
      //       After a restart, the interface position w.r.t which we stored fx_sum_ before the restart and
      //       the interface position for the first fluid-evaluate call are the same (then we can ensure the same dofset sorting).
      //       At the beginning of a new time step fx_sum_ and x_sum are created based on the same DofRowMap.
      Extractor().InsertVector(fx_sum_,fluid_block_, x_sum_);

      // if not a new time-step, take the step-increment which has been summed up so far
      if (HaveAle() && ax_sum_ != Teuchos::null) Extractor().AddVector(ax_sum_, ale_i_block_, x_sum_);
    }

    /*----------------------------------------------------------------------*/
    // Setup the new linear system and solve it and check convergence
    /*----------------------------------------------------------------------*/

    //-------------------
    //Build the linear system
    //J^{n+1}(x_i) \Delta^{n+1}_{x_i+1}=-r^{n+1}(x_i)
    //i: Newton iteration counter
    //J: Jacobian
    //r: RHS-vector
    //-------------------
    SetupSystemMatrix();

    if ( not systemmatrix_->Filled() ) dserror("Unfilled system matrix! Fatal error!");

    // Create the RHS consisting of the field residuals and coupling term residuals
    SetupRHS();


    //-------------------
    // Apply Dirichlet BCs to the whole system
    //-------------------
    ApplyDBC();


    //-------------------
    // Solver call
    //-------------------
    LinearSolve();


    //-------------------
    //Build residual and incremental norms, count the DOFs!
    //-------------------
    BuildCovergenceNorms();


    //-------------------
    //Give some output
    //-------------------
    PrintNewtonIter();



    //Increment loop index
    iter_+=1;

  }//End of Newton loop!

  //After the loop exit, the iteration counter is 1 higher than the true no. of
  //iterations! Correct that:
  iter_-=1;

  // Note:
  // the iteration increment computed at the latest will not be added finally,
  // in doing so, also the lambda-forces correspond to the iteration before.
  // When the residual was small enough and the new iteration increment is also sufficient small this is fine.
  // The lambda-forces and the current iterations of fluid and solid correspond to each other and yield to a global residual
  // which is smaller than the required tolerance


  /*----------------------------------------------------------------------*/
  // print converged/non-converged info
  /*----------------------------------------------------------------------*/
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
      IO::cout << "----------------------------------------- Newton not converged in ITEMAX iterations ! --------------------------------------" << IO::endl;

      if(iter_outer_ < itermax_outer_) // just in case that another restart will be performed!
        IO::cout << "- WARNING: increase the number nonlinear Newton-iterations, the additional restart does not help but solves the same system twice!!! -" << IO::endl;
    }
    return false;
  }

  return false;
}  // NewtonFull()



/*----------------------------------------------------------------------*/
// compute all norms used for convergence check
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

  // build full residual norms
  rhs_->Norm2(&normrhs_);

  // structural Dofs
  Extractor().ExtractVector(rhs_,0)->Norm2(&normstrrhsL2_);
  Extractor().ExtractVector(rhs_,0)->NormInf(&normstrrhsInf_);

  // fluid velocity Dofs
  fluidvelpresextract.ExtractVector(Extractor().ExtractVector(rhs_,1),0)->Norm2(&normflvelrhsL2_);
  fluidvelpresextract.ExtractVector(Extractor().ExtractVector(rhs_,1),0)->NormInf(&normflvelrhsInf_);

  // fluid pressure Dofs
  fluidvelpresextract.ExtractVector(Extractor().ExtractVector(rhs_,1),1)->Norm2(&normflpresrhsL2_);
  fluidvelpresextract.ExtractVector(Extractor().ExtractVector(rhs_,1),1)->NormInf(&normflpresrhsInf_);


  //-------------------------------
  // build solution increment norms
  //-------------------------------

  // build full increment norm
  iterinc_->Norm2(&norminc_);

  // structural Dofs
  Extractor().ExtractVector(iterinc_,0)->Norm2(&normstrincL2_);
  Extractor().ExtractVector(iterinc_,0)->NormInf(&normstrincInf_);

  // fluid velocity Dofs
  fluidvelpresextract.ExtractVector(Extractor().ExtractVector(iterinc_,1),0)->Norm2(&normflvelincL2_);
  fluidvelpresextract.ExtractVector(Extractor().ExtractVector(iterinc_,1),0)->NormInf(&normflvelincInf_);

  // fluid pressure Dofs
  fluidvelpresextract.ExtractVector(Extractor().ExtractVector(iterinc_,1),1)->Norm2(&normflpresincL2_);
  fluidvelpresextract.ExtractVector(Extractor().ExtractVector(iterinc_,1),1)->NormInf(&normflpresincInf_);


  //-------------------------------
  //get length of the structural and fluid vector
  //-------------------------------
  ns_   = (*(Extractor().ExtractVector(rhs_,0))).GlobalLength();                                      //structure
  nf_   = (*(Extractor().ExtractVector(rhs_,1))).GlobalLength();                                      //fluid
  nfv_  = (*(fluidvelpresextract.ExtractVector(Extractor().ExtractVector(rhs_,1),0))).GlobalLength(); //fluid velocity
  nfp_  = (*(fluidvelpresextract.ExtractVector(Extractor().ExtractVector(rhs_,1),1))).GlobalLength(); //fluid pressure
  nall_ = (*rhs_).GlobalLength();                                                                     //all


}



/*----------------------------------------------------------------------*
 * update the global step-increment, evaluate the single fields with
 * x^n+1 with x^n+1 = x_n + stepinc and return if the fluid dofsets
 * between the two last iterations changed and
 * a Newton restart is necessary                            schott 08/14
 *----------------------------------------------------------------------*/
bool FSI::MonolithicXFEM::Evaluate()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::Evaluate");


  // ------------------------------------------------------------------
  // ------------------------------------------------------------------
  // Update the global step-increment from the last Newton-solve and extract the single field step-increments
  // ------------------------------------------------------------------
  // ------------------------------------------------------------------

  // Structure and fluid fields (evaluate-calls) expect the respective step-increment (x^n+1_i+1 - x^n).
  // So we add all of the increments together to build the step increment.
  //
  // The update of the latest step increment with iteration increments:
  // x^n+1_i+1 = x^n+1_i + iterinc with x the current step increment

  // step-increments for single fields:
  // sx contains the current step increment w.r.t. Pred(disp(t^n)) for the structure block
  // fx contains the current step increment w.r.t. t^n from the last Newton restart for the fluid block
  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> ax;


  // update the whole step-increment vector
  // note: for iter_ = 1 the global x_sum_ is not available yet, then we take the single step-increments
  if(iter_ > 1)
  {
    // update the step-increment
    x_sum_->Update(1.0,*iterinc_,1.0);

    // extract the single field step-increments from the global step-increment
    ExtractFieldVectors(x_sum_,sx,fx,ax);
  }
  else
  {
    sx = sx_sum_; // take the sx from before the restart
    fx = fx_sum_; // take the fx from before the restart
    if (HaveAle())
      ax = ax_sum_; // take the ax from before the restart
  }

  // TODO: check this update, can we place it anywhere such that it is more consistent with the fluid-specific fx_sum?
  // save the current structural step-increment
  sx_sum_ = sx;
  if (HaveAle())
    ax_sum_ = ax;

  //TODO:
  if (sdbg_!=Teuchos::null)
  {
    sdbg_->NewIteration();
    sdbg_->WriteVector("x",*StructurePoro()->Interface()->ExtractFSICondVector(sx));
  }



  // ------------------------------------------------------------------
  // ------------------------------------------------------------------
  // Call all fields evaluate method and assemble fields rhs vectors and matrices
  // ------------------------------------------------------------------
  // ------------------------------------------------------------------

  //StructureField()->writeGmshStrucOutputStep();
  //-------------------
  // structure field
  //-------------------
  {
    Comm().Barrier();

    // structural field
    Epetra_Time ts(Comm());

    // ------------------------------------------------------------------
    // ------------------------------------------------------------------
    // Set Field State Section, here we should set the state with the step increments in all fields (atm this is just done for ALE)
    // ------------------------------------------------------------------
    // ------------------------------------------------------------------
    if (HaveAle() && ax != Teuchos::null) //we should move this into the ALE Field!
    {
      Teuchos::RCP<Epetra_Vector> DispnpAle = Teuchos::rcp(new Epetra_Vector(*AleField()->DofRowMap()),true);
      DispnpAle->Update(1.0,*AleField()->Interface()->InsertOtherVector(ax),1.0,*AleField()->Dispn(),0.0); //update ale disp here...
      AleField()->GetDBCMapExtractor()->InsertOtherVector(AleField()->GetDBCMapExtractor()->ExtractOtherVector(DispnpAle),AleField()->WriteAccessDispnp()); //just update displacements which are not on dbc condition
    }

    // call the structure evaluate with the current step increment  Delta d = d^(n+1,i+1) - d^n
    if (sx == Teuchos::null)
      sx = Teuchos::rcp(new Epetra_Vector(*StructurePoro()->DofRowMap(),true));
    StructurePoro()->Evaluate(sx, iter_ == 1);

    if(Comm().MyPID() == 0) IO::cout  << "structure time: " << ts.ElapsedTime() << IO::endl;


  }

  //--------------------------------------------------------
  // permute the fluid step-inc (ordered w.r.t. restart state) to current dofset-state
  // nothing has to be permuted for the first run as the following first evaluate call will fix the reference dofset for this Newton loop
  // nothing has to be permuted before we call the fluid-evaluate the second time, since the second call we determine if dofsets have permuted
  // the first potentially valid permutation is set and available after the second fluid-evaluate call
  //--------------------------------------------------------

  Teuchos::RCP<Epetra_Vector> fx_permuted = Teuchos::null;

  if(fx != Teuchos::null)
  {
    fx_permuted = Teuchos::rcp(new Epetra_Vector(*fx));

    PermuteFluidDOFSForward(fx_permuted);
  }



  // update fluid field increments
  FluidField()->UpdateByIncrements(fx_permuted);

  //StructurePoro()->StructureField()->writeGmshStrucOutputStep();

  // ------------------------------------------------------------------
  // set the current interface displacement to the fluid field to be used in the cut
  // ------------------------------------------------------------------


  // update coupling objects and conditionmanager

  for (std::map<int, Teuchos::RCP<XFEM::Coupling_Manager> >::iterator coupit = coup_man_.begin(); coupit != coup_man_.end(); ++coupit)
    coupit->second->SetCouplingStates();

  // update ALE
  if (HaveAle())
    AleField()->Evaluate();


  //-------------------
  // fluid field
  //-------------------
  // * update ivelnp-vector with the step-increment
  //
  // PrepareSolve();
  // * cut at new interface position
  // * create new state-vectors and systemmatrix and
  // * perform time-integration to obtain a new reference solution of veln at t^n w.r.t. new interface position
  // * perform a pseudo-time-integration to map the current iteration velnp (u^(n+1,i+1)) to new interface position
  // * TODO: possibly perform fluid predictor
  // * update old-rhs and
  // * evaluate Neumann and DBCs
  //
  // Evaluate:
  // * call evaluate routine to assemble fluid rhs and systemmatrix
  //
  {
    Comm().Barrier();

    //Teuchos::TimeMonitor::zeroOutTimers();

    // fluid field
    Epetra_Time tf(Comm());

    // call the fluid evaluate with the current time-step-increment w.r.t. u^n from the restart:
    // Delta(u,p) = (u,p)^(n+1,i+1) - u^n
    // For the first call of a time-step, call Evaluate with a null-pointer
    // note: call the fluid with the permuted step-increment vector as the fluid-dofsets can permute during the Newton
    // whereas the x_sum_ has to preserve the order of dofs during the Newton
    FluidField()->Evaluate();


    if(Comm().MyPID() == 0) IO::cout << "fluid time : " << tf.ElapsedTime() << IO::endl;

    //Teuchos::TimeMonitor::summarize();

  }

  //--------------------------------------------------------
  // update permutation cycles how to permute fluid dofs during the Newton
  //--------------------------------------------------------

  // update the permutation map used for permuting fluid-dofs between Dofset after restarting the Newton.
  // Build a vector of cycles how to permute dofs between the reference dofset from restarting the Newton and current dofset
  // note: No permutation for the first call since restarting the Newton - this is the reference dofset
  //       the first potentially valid permutation is set and available after the second fluid-evaluate call
  if(iter_ > 1)
  {
    UpdatePermutationMap(*FluidField()->GetPermutationMap());
  }

  //-------------------
  // check for changing dofsets compared to the last Newton iteration to decide if the Newton has to get restarted or continued
  //-------------------


  if(FluidField()->NewtonRestartMonolithic()) return true;



  return false; // continue with the setup of the new system and solving the system
}


/*----------------------------------------------------------------------*
 | check convergence of Newton iteration (public)          schott 08/14 |
 *----------------------------------------------------------------------*/
bool FSI::MonolithicXFEM::Converged()
{
  // check for single norms (increment, residual)
  bool convinc  = false; // increment converged?
  bool convfres = false; // residual converged?

  //---------------------------------------------
  // structural and fluid increments
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

  //---------------------------------------------
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

  //---------------------------------------------
  // combined increment + residual check?
  bool converged = false;

  if (combincfres_==INPAR::FSI::bop_and)
    converged = (convinc and convfres);
  else
    dserror("Just binary operator and for convergence check of Newton increment and residual supported!");

  return converged;
}  // Converged()


/*----------------------------------------------------------------------
 | update the permutation map between using the recent permutations    |
 | between the last two Newton iterations                 schott 08/14 |
 ----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::UpdatePermutationMap(
    std::map<int, int> permutation_map /// permutation map between last two Newton iterations, by copy, do not call by reference
)
{
  //TODO: remove the counter and the screen output
  int count_updates = 0;
  int removed_permutations = 0;

  //--------------------------------
  // update the permutation map
  //--------------------------------

  // first, look if one of the already existing permutations have to be updated
  for(std::map<int,int>::iterator p =  permutation_map_.begin();
      p != permutation_map_.end();
      p++)
  {
    // check if there is an update in the recent permutation map available
    std::map<int,int>::iterator p_recent_it = permutation_map.find(p->second);
    // update available
    if(p_recent_it != permutation_map.end())
    {
      // update the target of the permutation
      p->second = p_recent_it->second;
      // delete the update permutation, this has not to be considered again
      permutation_map.erase(p_recent_it);

      if(p->first == p->second)
      {
        // permutation became obsolete, remove it
        std::map<int,int>::iterator tmp_it = p; // save the current iterator before removing and increase it
        tmp_it++;

        // erase the permutation as it is done
        permutation_map_.erase(p);
        removed_permutations++;

        p=tmp_it;
      }

      count_updates++;
    }
  }

#if(1)
  std::cout << " adapted permutations        " << count_updates          << std::endl;
  std::cout << " removed permutations        " << removed_permutations   << std::endl;
  std::cout << " new additional permutations " << permutation_map.size() << std::endl;
#endif

  // second, we have to add all new additional permutations
  permutation_map_.insert(permutation_map.begin(), permutation_map.end());

  //--------------------------------
  // build new permutation cycles based on the updated permutation map
  //--------------------------------

  // finally build the new permutation cycles used for permuting the dofs forward and backward
  BuildFluidPermutation();
}


/*----------------------------------------------------------------------
 | build the new permutation cycles used for permuting the             |
 | fluid dofs forward and backward                        schott 08/14 |
 ----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::BuildFluidPermutation()
{

  /// vector of permutation cycles, one cycle consists of a series of global dof ids (gids)
  permutation_.clear();

  // make a copy of the internal permutation_map which we will modify here
  std::map<int, int> permutation_map = permutation_map_;

  // the order in a permutation cycle describes how dofs w.r.t old interface position can be mapped to new interface position
  // by reading the cycle from beginning to the end


  // permutation cycle of ghost dofs (global ghost DOF ids (gids) w.r.t one node: a b c d)
  //      a                              b
  //      b           ---  P --->>       d       -----> (1, 3, 4, 2) (,1) = cycle
  //      c         <<---P^(-1)---       a
  //      d                              c
  //
  // old gids ---> permuation_map ---> new gids ( forward  permutation )
  // new gids <--- permuation_map <--- new gids ( backward permutation )

  // build permutation cycles from the permutation map
  // in the permutation map one-to-one relations between gids are stored (key = old gid, value = new gid)
  for(std::map<int,int>::iterator map_it=permutation_map.begin(); map_it!=permutation_map.end(); map_it++)
  {

    if(map_it->second == -1) continue; // mapping already done by another cycle

    // create a new permutation cycle
    std::vector<int> new_cycle;

    int start_gid = map_it->first;

    // initialize the next-iterator for iterating the cycle to be created by the current start iteration of the permutation map
    std::map<int,int>::iterator next_it = map_it;

    // create the cycle
    while(next_it != permutation_map.end())
    {
      // add the current stored gid to the cycle
      new_cycle.push_back(next_it->first);
      // the next stored gid in the cycle
      int next_gid = next_it->second;

      // mark the entry as done, do not consider this single permutation again
      next_it->second = -1;

      // we reached the end of the cycle, stop the new cycle here
      if(next_gid == start_gid) break;

      // jump to the next entry
      next_it = permutation_map.find(next_gid);
    }

#if(0)
    // print cycle

    std::cout << " new cycle: " << std::endl;

    for(int i=0; i< (int)new_cycle.size(); i++)
      std::cout << " " << new_cycle[i];

    std::cout << std::endl;
#endif

    if((int)new_cycle.size() > 2) dserror("this is the first time that we permute more than two ghost dofsets! Check if the implementation works properly!");

    // new cycle
    permutation_.push_back(new_cycle);
  }
}


/*----------------------------------------------------------------------
 | forward permutation of fluid dofs -                                 |
 | transform vectors (based on dofsets) w.r.t old interface position   |
 | forward to a vector (based on dofsets) w.r.t. new interface position|
 |                                                        schott 08/14 |
 ----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::PermuteFluidDOFSForward(
    Teuchos::RCP< Epetra_Vector>&  fx
)
{

  //---------------------------------
  // forward permutation of dofsets
  //---------------------------------
  // transform vectors (based on dofsets) w.r.t old interface position forward to a vector (based on dofsets) w.r.t. new interface position

  // loop all permutation cycles
  for(std::vector<std::vector<int> >::iterator i=permutation_.begin(); i!=permutation_.end(); i++)
  {
    // permutation cycle of ghost dofs (global ghost DOF ids (gids) w.r.t one node: a b c d)
    //      a                              b
    //      b           ---  P --->>       d       -----> (1, 3, 4, 2) (,1) = cycle
    //      c                              a
    //      d                              c
    //
    // old gids ---> permuation_map ---> new gids ( forward  permutation )

    std::vector<int>& p_cycle = *i;

    double tmp_value = 0.0;

    // first  value -- to -- second position
    // second value -- to -- third  position
    // ...
    // last   value -- to -- first position
    for(std::vector<int>::iterator key=p_cycle.begin(); key!=p_cycle.end(); key++)
    {
      if(key+1 != p_cycle.end()) // standard during the cycle
      {
        tmp_value = (*fx)[fx->Map().LID(*(key+1))]; // save the value before it will be overwritten
        (*fx)[fx->Map().LID(*(key+1))] = (*fx)[fx->Map().LID(*(key))]; // set current value to next position
        //std::cout << "copy value from gid " << *(key) << " to " << *(key+1) << std::endl;
      }
      else // last value in cycle reached
      {
        (*fx)[fx->Map().LID(*p_cycle.begin())] = tmp_value;
        //std::cout << "copy value from tmp to " << *p_cycle.begin() << std::endl;
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------
 | backward permutation of fluid dofs -                                |
 | transform vectors (based on dofsets) w.r.t new interface position   |
 | backward to a vector (based on dofsets) w.r.t. new interface        |
 | position                                               schott 08/14 |
 ----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::PermuteFluidDOFSBackward(
    Teuchos::RCP< Epetra_Vector>&  fx
)
{

  //---------------------------------
  // backward permutation of dofsets
  //---------------------------------
  // transform vectors (based on dofsets) w.r.t new interface position backward to a vector (based on dofsets) w.r.t. old interface position


  // loop all permutation cycles
  for(std::vector<std::vector<int> >::iterator i=permutation_.begin(); i!=permutation_.end(); i++)
  {
    // permutation cycle of ghost dofs (global ghost DOF ids (gids) w.r.t one node: a b c d)
    //      a                              b
    //      b                              d       -----> (1, 3, 4, 2) (,1) = cycle
    //      c         <<---P^(-1)---       a
    //      d                              c
    //
    // new gids <--- permuation_map <--- new gids ( backward permutation )

    std::vector<int>& p_cycle = *i;

    double tmp_value = 0.0;

    //  last    value -- to -- (last-1) position
    // (last-1) value -- to -- (last-2) position
    // ...
    //  first   value -- to -- last position
    for(std::vector<int>::iterator key_reverse=p_cycle.end()-1; (key_reverse+1)!=p_cycle.begin(); key_reverse--)
    {

      if(key_reverse != p_cycle.begin()) // standard during the cycle
      {
        tmp_value = (*fx)[fx->Map().LID(*(key_reverse-1))]; // save the value before it will be overwritten
        (*fx)[fx->Map().LID(*(key_reverse-1))] = (*fx)[fx->Map().LID(*(key_reverse))]; // set current value to position before
        //std::cout << "copy value from gid " << *(key_reverse) << " to " << *(key_reverse-1) << std::endl;
      }
      else
      {
        (*fx)[fx->Map().LID(*(p_cycle.end()-1))] = tmp_value;
        //std::cout << "copy value from tmp to " << *(p_cycle.end()-1) << std::endl;
      }
    }
  }

  return;
}



/*----------------------------------------------------------------------*
 | create linear solver                            schott/wiesner 10/14 |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::CreateLinearSolver()
{
  // get the solver number used for linear XFSI solver
//  const int linsolvernumber = fsidyn_.get<int>("LINEAR_SOLVER");
  //TODO: get via input file, no LINEAR_SOLVER in FSI-Dynamic so far...
  const int linsolvernumber = 1;
  // check if the XFSI solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for monolithic XFSI. Please set LINEAR_SOLVER in XFSI DYNAMIC to a valid number!");

  // get solver parameter list of linear XFSI solver
  const Teuchos::ParameterList& xfsisolverparams = DRT::Problem::Instance()->SolverParams(linsolvernumber);

  // safety check if the hard-coded solver number is the XFSI-solver
  if(xfsisolverparams.get<std::string>("NAME") != "XFSI_SOLVER")
    dserror("check whether solver with number 1 is the XFSI_SOLVER and has this name!");


  const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>( xfsisolverparams, "SOLVER" );

  //----------------------------------------------
  // create direct solver for merged block matrix
  //----------------------------------------------
  if( solvertype == INPAR::SOLVER::umfpack )
  {
    if (Comm().MyPID() == 0) std::cout << "Merged XFSI block matrix is used!\n" << std::endl;

    merge_fsi_blockmatrix_ = true;

    Teuchos::RCP<Teuchos::ParameterList> solverparams = Teuchos::rcp(new Teuchos::ParameterList);
    solverparams->set("solver","umfpack");

    solver_ = Teuchos::rcp(new LINALG::Solver( solverparams, Comm(), DRT::Problem::Instance()->ErrorFile()->Handle() ));

    return;
  }

  //----------------------------------------------
  // create iterative solver for XFSI block matrix
  //----------------------------------------------

  if ( (solvertype != INPAR::SOLVER::aztec_msr) and (solvertype != INPAR::SOLVER::belos) )
  {
    dserror("aztec solver expected");
  }


  // get parameter list of structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // use solver blocks for structure
  // get the solver number used for structural solver
  const int slinsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (slinsolvernumber == (-1))
    dserror("no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

  // get parameter list of thermal dynamics
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  // use solver blocks for temperature (thermal field)
  // get the solver number used for thermal solver
  const int flinsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
  // check if the XFSI solver has a valid solver number
  if (flinsolvernumber == (-1))
    dserror("no linear solver defined for fluid field. Please set FLUID_SOLVER in FLUID DYNAMIC to a valid number!");


  const int azprectype = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>( xfsisolverparams, "AZPREC" );

  // plausibility check
  switch (azprectype)
  {
  case INPAR::SOLVER::azprec_BGS2x2 :
    break;
  case INPAR::SOLVER::azprec_BGSnxn :
  case INPAR::SOLVER::azprec_TekoSIMPLE :
  {
#ifdef HAVE_TEKO
    // check if structural solver and thermal solver are Stratimikos based (Teko expects stratimikos)
    int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(slinsolvernumber), "SOLVER");
    if ( (solvertype != INPAR::SOLVER::stratimikos_amesos) and
         (solvertype != INPAR::SOLVER::stratimikos_aztec) and
         (solvertype != INPAR::SOLVER::stratimikos_belos)
       )
    dserror("Teko expects a STRATIMIKOS solver object in STRUCTURE SOLVER");

    solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(flinsolvernumber), "SOLVER");
    if ( (solvertype != INPAR::SOLVER::stratimikos_amesos) and
         (solvertype != INPAR::SOLVER::stratimikos_aztec) and
         (solvertype != INPAR::SOLVER::stratimikos_belos)
       )
      dserror("Teko expects a STRATIMIKOS solver object in fluid solver %3d",flinsolvernumber);
#else
    dserror("Teko preconditioners only available with HAVE_TEKO flag for TRILINOS_DEV (>Q1/2011)");
#endif
    break;
  }
  case INPAR::SOLVER::azprec_MueLuAMG_sym :
  case INPAR::SOLVER::azprec_AMGnxn :
  case INPAR::SOLVER::azprec_CheapSIMPLE :
  {
    // no plausibility checks here
    // if you forget to declare an xml file you will get an error message anyway
  }
  break;
  default:
    dserror("Block Gauss-Seidel BGS2x2 preconditioner expected. Alternatively you can define your own AMG block preconditioner (using an xml file). This is experimental.");
    break;
  }


  // prepare linear solvers and preconditioners
  switch (azprectype)
  {
  case INPAR::SOLVER::azprec_BGS2x2 :
  case INPAR::SOLVER::azprec_BGSnxn :
  case INPAR::SOLVER::azprec_TekoSIMPLE :
  case INPAR::SOLVER::azprec_AMGnxn :
  case INPAR::SOLVER::azprec_CheapSIMPLE :
  {
    // This should be the default case (well-tested and used)
    solver_ = Teuchos::rcp(new LINALG::Solver(
                                 xfsisolverparams,
                                 // ggfs. explizit Comm von STR wie lungscatra
                                 Comm(),
                                 DRT::Problem::Instance()->ErrorFile()->Handle()
                                 )
                );

    // use solver blocks for structure and fluid
    const Teuchos::ParameterList& ssolverparams = DRT::Problem::Instance()->SolverParams(slinsolvernumber);
    const Teuchos::ParameterList& fsolverparams = DRT::Problem::Instance()->SolverParams(flinsolvernumber);

    solver_->PutSolverParamsToSubParams("Inverse1", ssolverparams);
    solver_->PutSolverParamsToSubParams("Inverse2", fsolverparams);

    // prescribe rigid body modes
    StructurePoro()->Discretization()->ComputeNullSpaceIfNecessary(
        solver_->Params().sublist("Inverse1")
        );
    FluidField()->Discretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse2")
      );

    if(azprectype==INPAR::SOLVER::azprec_CheapSIMPLE)
    {
    // Tell to the LINALG::SOLVER::SimplePreconditioner that we use the general implementation
      solver_->Params().set<bool>("GENERAL",true);
    }

    break;
  }
  case INPAR::SOLVER::azprec_MueLuAMG_sym:
  {
    solver_ = Teuchos::rcp(new LINALG::Solver(
                                 xfsisolverparams,
                                 // ggfs. explizit Comm von STR wie lungscatra
                                 Comm(),
                                 DRT::Problem::Instance()->ErrorFile()->Handle()
                                 )
                );

    // use solver blocks for structure and fluid
    const Teuchos::ParameterList& ssolverparams = DRT::Problem::Instance()->SolverParams(slinsolvernumber);
    const Teuchos::ParameterList& fsolverparams = DRT::Problem::Instance()->SolverParams(flinsolvernumber);

    // This is not very elegant:
    // first read in solver parameters. These have to contain ML parameters such that...
    solver_->PutSolverParamsToSubParams("Inverse1", ssolverparams);
    solver_->PutSolverParamsToSubParams("Inverse2", fsolverparams);

    // ... BACI calculates the null space vectors. These are then stored in the sublists
    //     Inverse1 and Inverse2 from where they...
    StructurePoro()->Discretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse1")
      );
    FluidField()->Discretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse2")
      );

    // ... are copied from here to ...
    const Teuchos::ParameterList& inv1source = solver_->Params().sublist("Inverse1").sublist("ML Parameters");
    const Teuchos::ParameterList& inv2source = solver_->Params().sublist("Inverse2").sublist("ML Parameters");

    // ... here. The "MueLu Parameters" sublists "Inverse1" and "Inverse2" only contain the basic
    //     information about the corresponding null space vectors, which are actually copied ...
    Teuchos::ParameterList& inv1 = solver_->Params().sublist("MueLu Parameters").sublist("Inverse1");
    Teuchos::ParameterList& inv2 = solver_->Params().sublist("MueLu Parameters").sublist("Inverse2");

    // ... here.
    inv1.set<int>("PDE equations", inv1source.get<int>("PDE equations"));
    inv2.set<int>("PDE equations", inv2source.get<int>("PDE equations"));
    inv1.set<int>("null space: dimension", inv1source.get<int>("null space: dimension"));
    inv2.set<int>("null space: dimension", inv2source.get<int>("null space: dimension"));
    inv1.set<double*>("null space: vectors", inv1source.get<double*>("null space: vectors"));
    inv2.set<double*>("null space: vectors", inv2source.get<double*>("null space: vectors"));
    inv1.set<Teuchos::RCP<std::vector<double> > >("nullspace", inv1source.get<Teuchos::RCP<std::vector<double> > >("nullspace"));
    inv2.set<Teuchos::RCP<std::vector<double> > >("nullspace", inv2source.get<Teuchos::RCP<std::vector<double> > >("nullspace"));

    //TODO: muelu for XFSI similar to TSI?
    dserror("MueLu for XFSI?");
    solver_->Params().sublist("MueLu Parameters").set("TSI",true);
    break;
  }
  default:
    dserror("Block Gauss-Seidel BGS2x2 preconditioner expected");
    break;
  }
}  // CreateLinearSolver()


/*----------------------------------------------------------------------*
 | solve linear FSI system                                 schott 07/13 |
 *----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::LinearSolve()
{
  if (Comm().MyPID() == 0) std::cout << " FSI::MonolithicXFEM::LinearSolve()" <<  std::endl;

  // Solve for inc_ = [disi_,tempi_]
  // Solve K_Teffdyn . IncX = -R  ===>  IncX_{n+1} with X=[d,(u,p)]
  // \f$x_{i+1} = x_i + \Delta x_i\f$

  // apply Dirichlet BCs to system of equations
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

  // default: use block matrix
  if (merge_fsi_blockmatrix_ == false)
  {
    // adapt solver tolerance
    if (solveradapttol_ and (iter_ > 1))
    {
      double worst = normrhs_;
      double wanted = tolrhs_;
      solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
    }

    // Infnormscaling: scale system before solving
    ScaleSystem(*systemmatrix_,*rhs_);

    // solve the problem, work is done here!
    solver_->Solve(
        systemmatrix_->EpetraOperator(),
        iterinc_,
        rhs_,
        true,     // refactorize the preconditioner?
        iter_==1  // build completely new solver including preconditioner
    );

    // Infnormscaling: unscale system after solving
    UnscaleSolution(*systemmatrix_,*iterinc_,*rhs_);


    //Adapt solver tolerance (important if Aztec solver is picked)
    //TODO: does or how does this work for changing Newton systems
    solver_->ResetTolerance();

  }  // use block matrix
  else // (merge_fsi_blockmatrix_ == true)
  {
    if(scaling_infnorm_) dserror("infnorm-scaling of FSI-system not supported for direct solver");

    //------------------------------------------
    // merge blockmatrix to SparseMatrix and solve
    Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

    //------------------------------------------
    // standard solver call
    solver_->Solve(
        sparse->EpetraOperator(),
        iterinc_,
        rhs_,
        true,
        iter_==1
    );
  }  // MergeBlockMatrix

  //TODO: can we do this?!
  // reset the solver (frees the pointer to the LINALG:: matrix' EpetraOperator and vectors also!)
  //std::cout << "reset the solver" << std::endl;
  solver_->Reset();

//  {
//
//  std::cout << "------------------ before permuting the solution--------------------" << std::endl;
//
//  Teuchos::RCP<Epetra_Vector> f_rhs_tmp = Extractor().ExtractVector(rhs_, 1);
//
//  std::cout << "rhs at dof 33174: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33174)] << std::endl;
//  std::cout << "rhs at dof 33175: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33175)] << std::endl;
//  std::cout << "rhs at dof 33176: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33176)] << std::endl;
//  std::cout << "rhs at dof 33177: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33177)] << std::endl;
//
//  std::cout << "rhs at dof 33178: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33178)] << std::endl;
//  std::cout << "rhs at dof 33179: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33179)] << std::endl;
//  std::cout << "rhs at dof 33180: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33180)] << std::endl;
//  std::cout << "rhs at dof 33181: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33181)] << std::endl;
//
//  Teuchos::RCP<Epetra_Vector> f_iterinc_tmp = Extractor().ExtractVector(iterinc_, 1);
//
//  std::cout << "incr at dof 33174: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33174)] << std::endl;
//  std::cout << "incr at dof 33175: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33175)] << std::endl;
//  std::cout << "incr at dof 33176: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33176)] << std::endl;
//  std::cout << "incr at dof 33177: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33177)] << std::endl;
//
//  std::cout << "incr at dof 33178: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33178)] << std::endl;
//  std::cout << "incr at dof 33179: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33179)] << std::endl;
//  std::cout << "incr at dof 33180: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33180)] << std::endl;
//  std::cout << "incr at dof 33181: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33181)] << std::endl;
//
//  }

  //TODO: check the GMSH-Debug-output (control via input file)
#if 0
    {
      // DEBUG output
      Teuchos::RCP<ADAPTER::XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<ADAPTER::XFluidFSI>(FluidField(), true);

      // Remark that GMSH-output of rhs and iterinc have to be written w.r.t current fluid
      // dofsets based on the respective cut/vc-data before the vectors are permuted to the reference dofset ordering
      xfluid->GmshOutput("DEBUG_MONOLITIC_RESIDUAL",  Step(), iter_ , Extractor().ExtractVector(rhs_,1), Teuchos::null);
      xfluid->GmshOutput("DEBUG_MONOLITIC_INC",  Step(), iter_ , Extractor().ExtractVector(iterinc_,1), Teuchos::null);
    }
#endif



  //---------------------------------------------
  // permute the increment and rhs vector back to the reference configuration w.r.t which iterinc_ and rhs are defined

  Teuchos::RCP<Epetra_Vector> f_iterinc_permuted = Extractor().ExtractVector(iterinc_, 1);
  PermuteFluidDOFSBackward(f_iterinc_permuted);
  Extractor().InsertVector(*f_iterinc_permuted, 1, *iterinc_);


  Teuchos::RCP<Epetra_Vector> f_rhs_permuted = Extractor().ExtractVector(rhs_, 1);
  PermuteFluidDOFSBackward(f_rhs_permuted);
  Extractor().InsertVector(*f_rhs_permuted, 1, *rhs_);

  //---------------------------------------------

//  {
//  std::cout << "------------------ after permuting the solution--------------------" << std::endl;
//
//  Teuchos::RCP<Epetra_Vector> f_rhs_tmp = Extractor().ExtractVector(rhs_, 1);
//
//  std::cout << "rhs at dof 33174: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33174)] << std::endl;
//  std::cout << "rhs at dof 33175: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33175)] << std::endl;
//  std::cout << "rhs at dof 33176: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33176)] << std::endl;
//  std::cout << "rhs at dof 33177: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33177)] << std::endl;
//
//  std::cout << "rhs at dof 33178: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33178)] << std::endl;
//  std::cout << "rhs at dof 33179: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33179)] << std::endl;
//  std::cout << "rhs at dof 33180: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33180)] << std::endl;
//  std::cout << "rhs at dof 33181: " <<  (*f_rhs_tmp)[f_rhs_tmp->Map().LID(33181)] << std::endl;
//
//  Teuchos::RCP<Epetra_Vector> f_iterinc_tmp = Extractor().ExtractVector(iterinc_, 1);
//
//  std::cout << "incr at dof 33174: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33174)] << std::endl;
//  std::cout << "incr at dof 33175: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33175)] << std::endl;
//  std::cout << "incr at dof 33176: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33176)] << std::endl;
//  std::cout << "incr at dof 33177: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33177)] << std::endl;
//
//  std::cout << "incr at dof 33178: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33178)] << std::endl;
//  std::cout << "incr at dof 33179: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33179)] << std::endl;
//  std::cout << "incr at dof 33180: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33180)] << std::endl;
//  std::cout << "incr at dof 33181: " <<  (*f_iterinc_tmp)[f_iterinc_tmp->Map().LID(33181)] << std::endl;
//
//  }


  if (Comm().MyPID() == 0) { std::cout << " Solved" <<  std::endl; }

}  // LinearSolve()


/*----------------------------------------------------------------------*/
// apply infnorm scaling to linear block system            schott 10/14 |
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{

  if (scaling_infnorm_)
  {
    // The matrices are modified here. Do we have to change them back later on?

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);

    if (A->LeftScale(*srowsum_) or
        A->RightScale(*scolsum_) or
        mat.Matrix(0,1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");


    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b,0);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0))
      dserror("structure scaling failed");

    Extractor().InsertVector(*sx,0,b);
  }
}



/*----------------------------------------------------------------------*/
// undo infnorm scaling from scaled solution               schott 10/14 |
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  if (scaling_infnorm_)
  {

    Teuchos::RCP<Epetra_Vector> sy = Extractor().ExtractVector(x,0);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0))
      dserror("structure scaling failed");

    Extractor().InsertVector(*sy,0,x);

    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b,0);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0))
      dserror("structure scaling failed");

    Extractor().InsertVector(*sx,0,b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if (A->LeftScale(*srowsum_) or
        A->RightScale(*scolsum_) or
        mat.Matrix(0,1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");

  }
}



/*----------------------------------------------------------------------*
 | create combined Dirichlet boundary condition map,                    |
 | map containing the dofs with Dirichlet BC                            |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> FSI::MonolithicXFEM::CombinedDBCMap()
{
  Teuchos::RCP<const Epetra_Map > scondmap = StructurePoro()->CombinedDBCMap();
  const Teuchos::RCP<const Epetra_Map > fcondmap = FluidField()->GetDBCMapExtractor()->CondMap();

  Teuchos::RCP<Epetra_Map> condmap = LINALG::MergeMap(scondmap, fcondmap, false);

  return condmap;
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

   IO::cout << "============================================================================================================================="<< IO::endl;

   // enter converged state etc
   IO::cout << "|outerit";
   IO::cout << "|  nit  |";

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
//   if (myrank_ == 0)
//   {
//     if (itnum>0)
//     {
//       printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
//         itnum,itmax,ittol,vresnorm_,presnorm_,incvelnorm_L2_/velnorm_L2_,
//         incprenorm_L2_/prenorm_L2_);

   //TODO: komplette Überarbeitung von Rel vs Abs notwendig!!! siehe abs vs rel z.B. in Fluid-Code


   IO::cout << " " << iter_outer_ << "/" << itermax_outer_;
   IO::cout << " " << iter_       << "/" << itermax_;

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



 /*----------------------------------------------------------------------*/
 // read restart data for monolithic XFSI system
 /*----------------------------------------------------------------------*/
 void FSI::MonolithicXFEM::ReadRestart(int step)
 {
   //--------------------------------
   // read structural field
   StructurePoro()->ReadRestart(step);

   //--------------------------------
   // read fluid field
   FluidField()->ReadRestart(step);

   //--------------------------------
   // read ale field
   if (HaveAle())
     AleField()->ReadRestart(step);

   //--------------------------------
   // setup a new system as dofrowmaps could have been changed!
   SetupSystem();

   //--------------------------------
   // NOTE: do the following after StructureField()->ReadRestart and after FluidField()->ReadRestart
   // as ReadMesh can change the discretization and the dofrowmaps!!!

   // read Lagrange multiplier (ie forces onto the structure, Robin-type forces
   // consisting of fluid forces and the Nitsche penalty term contribution)
   IO::DiscretizationReader reader = IO::DiscretizationReader(StructurePoro()->Discretization(),step);
   for (std::map<int, Teuchos::RCP<XFEM::Coupling_Manager> >::iterator coupit = coup_man_.begin(); coupit != coup_man_.end(); ++coupit)
     coupit->second->ReadRestart(reader);
   //

   SetTimeStep(FluidField()->Time(),FluidField()->Step());
   return;
 }
