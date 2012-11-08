/*----------------------------------------------------------------------*/
/*!
 \file poroelast_monolithic.cpp

 \brief  Basis of all monolithic poroelasticity algorithms

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                                          |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  vuong 01/12 |
 *----------------------------------------------------------------------*/
#include "poroelast_monolithic.H"
#include "poroelast_defines.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_inpar/inpar_solver.H"

#include <Teuchos_TimeMonitor.hpp>
// needed for PrintNewton
#include <sstream>

// include this header for coupling stiffness terms
#include "../drt_lib/drt_assemblestrategy.H"

#include "../drt_io/io_control.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_adapter/ad_fld_poro.H"

#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_structure/stru_aux.H"

/*----------------------------------------------------------------------*
 | monolithic                                              vuong 01/12  |
 *----------------------------------------------------------------------*/
POROELAST::Monolithic::Monolithic(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams
    ) :
    PoroBase(comm,timeparams),
    printscreen_(true), // ADD INPUT PARAMETER
    printiter_(true), // ADD INPUT PARAMETER
    printerrfile_(true), // ADD INPUT PARAMETER FOR 'true'
    errfile_(DRT::Problem::Instance()->ErrorFile()->Handle()),
    zeros_(Teuchos::null),
    timer_(comm)
{
  const Teuchos::ParameterList& sdynparams
  = DRT::Problem::Instance()->StructuralDynamicParams();

  //some solver paramaters are red form the structure dynamic list (this is not the best way to do it ...)
  solveradapttol_= (DRT::INPUT::IntegralValue<int>(sdynparams, "ADAPTCONV") == 1);
  solveradaptolbetter_ = (sdynparams.get<double> ("ADAPTCONV_BETTER"));

  strmethodname_ = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams,"DYNAMICTYP");

  //  solver
#ifdef POROELASTBLOCKMATRIXMERGE
  // create a linear solver
  // get UMFPACK...
  // get dynamic section of poroelasticity
  const Teuchos::ParameterList& poroelastdyn = DRT::Problem::Instance()->PoroelastDynamicParams();
  // get the solver number used for linear poroelasticity solver
  const int linsolvernumber = poroelastdyn.get<int>("LINEAR_SOLVER");
  // check if the poroelasticity solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for poroelasticity. Please set LINEAR_SOLVER in POROELASTICITY DYNAMIC to a valid number!");
  const Teuchos::ParameterList& solverparams =
    DRT::Problem::Instance()->SolverParams(linsolvernumber);
  const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(
    solverparams, "SOLVER");
  if (solvertype != INPAR::SOLVER::umfpack)
    dserror("umfpack solver expected");

  solver_ = rcp(new LINALG::Solver( solverparams,
                                    Comm(),
                                    DRT::Problem::Instance()->ErrorFile()->Handle())
               );

#else
  //dserror("implicit solver not implemented");
  // create a linear solver
  CreateLinearSolver();
#endif

  // Get the parameters for the Newton iteration
  itermax_ = poroelastdyn.get<int> ("ITEMAX");
  itermin_ = poroelastdyn.get<int> ("ITEMIN");
  normtypeinc_ = DRT::INPUT::IntegralValue<INPAR::POROELAST::ConvNorm>(
      poroelastdyn, "NORM_INC");
  normtypefres_ = DRT::INPUT::IntegralValue<INPAR::POROELAST::ConvNorm>(
      poroelastdyn, "NORM_RESF");
  combincfres_ = DRT::INPUT::IntegralValue<INPAR::POROELAST::BinaryOp>(
      poroelastdyn, "NORMCOMBI_RESFINC");

  tolinc_ = poroelastdyn.get<double> ("INCTOL");
  tolfres_ = poroelastdyn.get<double> ("RESTOL");

}


/*----------------------------------------------------------------------*
 | solve time step                    vuong 01/12     |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::Solve()
{
  // counter and print header
  // predict solution of both field (call the adapter)
  PrepareTimeStep();

  // Newton-Raphson iteration
  NewtonFull();

  // calculate stresses, strains, energies
  PrepareOutput();

  // update all single field solvers
  Update();

  // write output to screen and files
  Output();

} // TimeLoop

/*----------------------------------------------------------------------*
 | solution with full Newton-Raphson iteration            vuong 01/12   |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::NewtonFull()
{

  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #rhs_ is the positive force residuum
  // --> On #systemmatrix_ is the effective dynamic tangent matrix

  //initialize norms and increment
  SetupNewton();

  //---------------------------------------------- iteration loop

  // equilibrium iteration loop (loop over k)
  while (((not Converged()) and (iter_ <= itermax_)) or (iter_ <= itermin_))
  {
    timer_.ResetStartTime();

    // compute residual forces #rhs_ and tangent #tang_
    // whose components are globally oriented
    // build linear system stiffness matrix and rhs/force residual for each
    // field, here e.g. for structure field: field want the iteration increment
    // 1.) Update(iterinc_),
    // 2.) EvaluateForceStiffResidual(),
    // 3.) PrepareSystemForNewtonSolve()
    Evaluate(iterinc_);

    // create the linear system
    // \f$J(x_i) \Delta x_i = - R(x_i)\f$
    // create the systemmatrix
    SetupSystemMatrix();

    // check whether we have a sanely filled tangent matrix
    if (not systemmatrix_->Filled())
    {
      dserror("Effective tangent matrix must be filled here");
    }

    // create full monolithic rhs vector
    SetupRHS(iter_==1);

    // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
    // is done in PrepareSystemForNewtonSolve() within Evaluate(iterinc_)
    LinearSolve();

    // reset solver tolerance
    solver_->ResetTolerance();

    //build norms
    BuildCovergenceNorms();

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;
  } // end equilibrium loop

  //---------------------------------------------- iteration loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ( (Converged()) and (Comm().MyPID()==0) )
  {
    PrintNewtonConv();
  }
  else if (iter_ >= itermax_)
  {
    dserror("Newton unconverged in %d iterations", iter_);
  }

  // recover Lagrange multiplier \lambda_\Gamma at the interface at the end of each time step
  // (i.e. condensed forces onto the structure) needed for rhs in next time step
  RecoverLagrangeMultiplier();

}    // NewtonFull()

/*----------------------------------------------------------------------*
 | evaluate the single fields                              vuong 01/12   |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::Evaluate");

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;

  // if an increment vector exists
  if (x != Teuchos::null)
  {
    // extract displacement sx and fluid fx incremental vector of global
    // unknown incremental vector x
    ExtractFieldVectors(x, sx, fx,iter_==1);
  }

  // Newton update of the fluid field
  // update velocities and pressures before passed to the structural field
  //  UpdateIterIncrementally(fx),
  FluidField().UpdateNewton(fx);

  // call all elements and assemble rhs and matrices
  /// structural field

  // structure Evaluate (builds tangent, residual and applies DBC)
  //Epetra_Time timerstructure(Comm());

  // apply current velocity and pressures to structure
  //fluidveln_ ->Update(1.0, *FluidField().Velnp(),0.0);
  SetFluidSolution();

  // Monolithic Poroelasticity accesses the linearised structure problem:
  //   UpdaterIterIncrementally(sx),
  //   EvaluateForceStiffResidual()
  //   PrepareSystemForNewtonSolve()
  StructureField()->Evaluate(sx);
  //cout << "  structure time for calling Evaluate: " << timerstructure.ElapsedTime() << "\n";

  /// fluid field

  // fluid Evaluate
  // (builds tangent, residual and applies DBC and recent coupling values)
  //Epetra_Time timerfluid(Comm());

  SetStructSolution();

  // monolithic Poroelasticity accesses the linearised fluid problem
  //   EvaluateRhsTangResidual() and
  //   PrepareSystemForNewtonSolve()
  FluidField().Evaluate(Teuchos::null);
  //cout << "  fluid time for calling Evaluate: " << timerfluid.ElapsedTime() << "\n";

} // Evaluate()

/*----------------------------------------------------------------------*
 | extract field vectors for calling Evaluate() of the       vuong 01/12|
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::ExtractFieldVectors(Teuchos::RCP<
    const Epetra_Vector> x, Teuchos::RCP<const Epetra_Vector>& sx,
    Teuchos::RCP<const Epetra_Vector>& fx,
    bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::ExtractFieldVectors");

  // process structure unknowns of the first field
  sx = Extractor().ExtractVector(x, 0);

  // process fluid unknowns of the second field
  fx = Extractor().ExtractVector(x, 1);
}


/*----------------------------------------------------------------------*
 | calculate velocities                                     vuong 01/12 |
 | like InterfaceVelocity(disp) in FSI::DirichletNeumann                |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> POROELAST::Monolithic::CalcVelocity(Teuchos::RCP<
    const Epetra_Vector> sx)
{
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::null;
  // copy D_n onto V_n+1
  vel = rcp(new Epetra_Vector(*(StructureField()->ExtractDispn())));
  // calculate velocity with timestep Dt()
  //  V_n+1^k = (D_n+1^k - D_n) / Dt
  vel->Update(1. / Dt(), *sx, -1. / Dt());

  return vel;
} // CalcVelocity()


/*----------------------------------------------------------------------*
 | setup system (called in porolast.cpp)                 vuong 01/12    |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::SetupSystem()
{
  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;

  // use its own DofRowMap, that is the 0th map of the discretization
  //
  // when using constraints applied via Lagrange-Multipliers there is a
  // difference between StructureField()->DofRowMap() and StructureField()->DofRowMap(0).
  // StructureField()->DofRowMap(0) returns the DofRowMap
  // known to the discretization (without lagrange multipliers)
  // while StructureField()->DofRowMap() returns the DofRowMap known to
  // the constraint manager (with lagrange multipliers)
  // In the constrained case we want the "whole" RowDofMap,
  // otherwise both calls are equivalent

  vecSpaces.push_back(StructureField()->DofRowMap());
  vecSpaces.push_back(FluidField().DofRowMap(0));

  if (vecSpaces[0]->NumGlobalElements() == 0)
    dserror("No structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements()==0)
    dserror("No fluid equation. Panic.");

  SetDofRowMaps(vecSpaces);

  // initialize Poroelasticity-systemmatrix_
  systemmatrix_ = rcp(new LINALG::BlockSparseMatrix<
      LINALG::DefaultBlockMatrixStrategy>(Extractor(), Extractor(), 81, false,
      true));

  k_sf_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(StructureField()->DofRowMap()), 81, true, true));
  k_fs_ = Teuchos::rcp(new LINALG::SparseMatrix(
                        *(FluidField().Discretization()->DofRowMap(0)),
                        //*(FluidField().DofRowMap()),
                        81, true, true));

  } // SetupSystem()

/*----------------------------------------------------------------------*
 | put the single maps to one full                                      |
 | Poroelasticity map together                              vuong 01/12 |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::SetDofRowMaps(const std::vector<Teuchos::RCP<
    const Epetra_Map> >& maps)
{
  fullmap_ = LINALG::MultiMapExtractor::MergeMaps(maps);

  // full Poroelasticity-blockmap
  blockrowdofmap_.Setup(*fullmap_, maps);
}

/*----------------------------------------------------------------------*
 | setup system matrix of poroelasticity                   vuong 01/12  |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::SetupSystemMatrix");

  // pure structural part k_ss (3nx3n)

  // build pure structural block k_ss
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  Teuchos::RCP<LINALG::SparseMatrix> k_ss = StructureField()->SystemMatrix();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  k_ss->UnComplete();

  // assign structure part to the TSI matrix
  mat.Assign(0, 0, View, *k_ss);

  /*----------------------------------------------------------------------*/
  // structural part k_sf (3nxn)
  // build mechanical-fluid block

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> k_sf = StructFluidCouplingMatrix();

  // call the element and calculate the matrix block
  ApplyStrCouplMatrix(k_sf);

  // Uncomplete mechanical-fluid matrix to be able to deal with slightly
  // defective interface meshes.
  k_sf->UnComplete();

  // assign coupling part to the Poroelasticity matrix
  mat.Assign(0, 1, View, *(k_sf));

  /*----------------------------------------------------------------------*/
  // pure fluid part k_ff ( (3n+1)x(3n+1) )

  // build pure fluid block k_ff
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  Teuchos::RCP<LINALG::SparseMatrix> k_ff = FluidField().SystemMatrix();

  if(nopencond_.size())
  {
    //Evaluate poroelasticity specific conditions
    EvaluateCondition(k_ff, Teuchos::null);
  }

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  k_ff->UnComplete();

  // assign fluid part to the poroelasticity matrix
  mat.Assign(1, 1, View, *(k_ff));

  /*----------------------------------------------------------------------*/
  // fluid part k_fs ( (3n+1)x3n )
  // build fluid-mechanical block

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> k_fs = FluidStructCouplingMatrix();

  // call the element and calculate the matrix block
  ApplyFluidCouplMatrix(k_fs);

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  k_fs->UnComplete();

  // assign coupling part to the Poroelasticity matrix
  mat.Assign(1, 0, View, *(k_fs));

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();

} // SetupSystemMatrix

/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                                 vuong 01/12  |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::SetupRHS(bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::SetupRHS");

  // create full monolithic rhs vector
  rhs_ = rcp(new Epetra_Vector(*DofRowMap(), true));

  // fill the Poroelasticity rhs vector rhs_ with the single field rhss
  SetupVector(*rhs_, StructureField()->RHS(), FluidField().RHS());
} // SetupRHS()


/*----------------------------------------------------------------------*
 | Solve linear Poroelasticity system                      vuong 01/12   |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::LinearSolve()
{
  // Solve for inc_ = [disi_,tempi_]
  // Solve K_Teffdyn . IncX = -R  ===>  IncX_{n+1} with X=[d,T]
  // \f$x_{i+1} = x_i + \Delta x_i\f$
  if (solveradapttol_ and (iter_ > 1))
  {
    double worst = normrhs_;
    double wanted = tolfres_;
    solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
  }

#ifdef POROELASTBLOCKMATRIXMERGE
  // merge blockmatrix to SparseMatrix and solve
  Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

  // apply Dirichlet BCs to system of equations
  iterinc_->PutScalar(0.0); // Useful? depends on solver and more
  LINALG::ApplyDirichlettoSystem(
      sparse,
      iterinc_,
      rhs_,
      Teuchos::null,
      zeros_,
      *CombinedDBCMap()
      );
  //  if ( Comm().MyPID()==0 ) { cout << " DBC applied to system" << endl; }

  if(nopencond_.size())
  {
    const Teuchos::RCP<const Epetra_Map >& nopenetrationmap = nopenetration_->Map(1);
    LINALG::ApplyDirichlettoSystem(iterinc_,rhs_,cond_rhs_,*nopenetrationmap);
  }

  // standard solver call
  solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, true, iter_ == 1);
  //  if ( Comm().MyPID()==0 ) { cout << " Solved" << endl; }

#else // use bgs2x2_operator
  //dserror("implicit solver not implemented");

  // apply Dirichlet BCs to system of equations
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

  // in case of inclined boundary conditions
  // rotate systemmatrix_ using GetLocSysTrafo()!=Teuchos::null
  LINALG::ApplyDirichlettoSystem(
    systemmatrix_,
    iterinc_,
    rhs_,
    Teuchos::null,
    zeros_,
    *CombinedDBCMap()
    );

  if(nopencond_.size())
  {
    const Teuchos::RCP<const Epetra_Map >& nopenetrationmap = nopenetration_.Map(1);
    LINALG::ApplyDirichlettoSystem(iterinc_,rhs_,cond_rhs_,*nopenetrationmap);
  }

  solver_->Solve(
             systemmatrix_->EpetraOperator(),
             iterinc_,
             rhs_,
             true,
             iter_==1
             );

#endif  // POROELASTBLOCKMATRIXMERGE
}

/*----------------------------------------------------------------------*
 | create linear solver                                   wiesner 07/11 |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::CreateLinearSolver()
{
  // get dynamic section
  const Teuchos::ParameterList& porodyn = DRT::Problem::Instance()->PoroelastDynamicParams();
  // get the solver number used for linear TSI solver
  const int linsolvernumber = porodyn.get<int>("LINEAR_SOLVER");
  // check if the TSI solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for monolithic Poroelasticity. Please set LINEAR_SOLVER in POROELASTICITY DYNAMIC to a valid number!");

  // get parameter list of structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // use solver blocks for structure
  // get the solver number used for structural solver
  const int slinsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (slinsolvernumber == (-1))
    dserror("no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

  // get parameter list of fluid dynamics
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  // use solver blocks for fluid
  // get the solver number used for fluid solver
  const int flinsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
  // check if the fluid solver has a valid solver number
  if (flinsolvernumber == (-1))
    dserror("no linear solver defined for fluid field. Please set LINEAR_SOLVER in FLUID DYNAMIC to a valid number!");

  // get solver parameter list of linear Poroelasticity solver
  const Teuchos::ParameterList& porosolverparams
    = DRT::Problem::Instance()->SolverParams(linsolvernumber);

  const int solvertype
    = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(
        porosolverparams,
        "SOLVER"
        );

  if (solvertype != INPAR::SOLVER::aztec_msr &&
      solvertype != INPAR::SOLVER::belos)
  {
    cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << " Note: the BGS2x2 preconditioner now "                  << endl;
    cout << " uses the structural solver and fluid solver blocks"  << endl;
    cout << " for building the internal inverses"                    << endl;
    cout << " Remove the old BGS PRECONDITIONER BLOCK entries "      << endl;
    cout << " in the dat files!"                                     << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << endl;
    dserror("aztec solver expected");
  }
  const int azprectype
    = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(
        porosolverparams,
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
      if (solvertype != INPAR::SOLVER::stratimikos_amesos &&
          solvertype != INPAR::SOLVER::stratimikos_aztec  &&
          solvertype != INPAR::SOLVER::stratimikos_belos)
      dserror("Teko expects a STRATIMIKOS solver object in STRUCTURE SOLVER");

      solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(flinsolvernumber), "SOLVER");
      if (solvertype != INPAR::SOLVER::stratimikos_amesos &&
          solvertype != INPAR::SOLVER::stratimikos_aztec  &&
          solvertype != INPAR::SOLVER::stratimikos_belos)
        dserror("Teko expects a STRATIMIKOS solver object in thermal solver %3d",flinsolvernumber);
#else
      dserror("Teko preconditioners only available with HAVE_TEKO flag (Trilinos >Q1/2011)");
#endif
    }
    break;
    default:
          dserror("Block Gauss-Seidel BGS2x2 preconditioner expected");
          break;
  }

  solver_ = rcp(new LINALG::Solver(
                          porosolverparams,
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
  StructureField()->Discretization()->ComputeNullSpaceIfNecessary(
                                       solver_->Params().sublist("Inverse1")
                                       );
  FluidField().Discretization()->ComputeNullSpaceIfNecessary(
                                    solver_->Params().sublist("Inverse2")
                                    );
}

/*----------------------------------------------------------------------*
 | initial guess of the displacements/velocities           vuong 01/12  |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::InitialGuess");

  // InitalGuess() is called of the single fields and results are put in TSI
  // increment vector ig
  SetupVector(*ig,
      // returns residual displacements \f$\Delta D_{n+1}^{<k>}\f$ - disi_
      StructureField()->InitialGuess(),
      // returns residual velocities or iterative fluid increment - incvel_
      FluidField().InitialGuess());
} // InitialGuess()

/*----------------------------------------------------------------------*
 | setup vector of the structure and fluid field            vuong 01/12|
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::SetupVector(Epetra_Vector &f,
                                        Teuchos::RCP<const Epetra_Vector> sv,
                                        Teuchos::RCP<const Epetra_Vector> fv)
{
  // extract dofs of the two fields
  // and put the structural/fluid field vector into the global vector f
  // noticing the block number

  Extractor().InsertVector(*sv, 0, f);
  Extractor().InsertVector(*fv, 1, f);
}

/*----------------------------------------------------------------------*
 | check convergence of Newton iteration (public)         vuong 01/12   |
 *----------------------------------------------------------------------*/
bool POROELAST::Monolithic::Converged()
{
  // check for single norms
  bool convinc = false;
  bool convfres = false;

  // residual increments
  switch (normtypeinc_)
  {
    case INPAR::POROELAST::convnorm_abs:
      convinc = norminc_ < tolinc_;
      break;
    default:
      dserror("Cannot check for convergence of residual values!");
      break;
  }

  // residual forces
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs:
    convfres = normrhs_ < tolfres_;
    break;
    default:
    dserror("Cannot check for convergence of residual forces!");
    break;
  }

  // combine increments and forces
  bool conv = false;
  if (combincfres_==INPAR::POROELAST::bop_and)
    conv = convinc and convfres;
  else if (combincfres_==INPAR::POROELAST::bop_or)
    conv = convinc or convfres;
  else
    dserror("Something went terribly wrong with binary operator!");

  // return things
  return conv;
}  // Converged()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file              |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::PrintNewtonIter()
{
  // print to standard out
  // replace myrank_ here general by Comm().MyPID()
  if ((Comm().MyPID() == 0) and printscreen_ and (Step()%printscreen_==0) and printiter_)
  {
    if (iter_ == 1)
      PrintNewtonIterHeader(stdout);
    PrintNewtonIterText(stdout);
  }

  // print to error file
  if (printerrfile_ and printiter_)
  {
    if (iter_ == 1)
      PrintNewtonIterHeader(errfile_);
    PrintNewtonIterText(errfile_);
  }

  return;
} // PrintNewtonIter()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file              |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::PrintNewtonIterHeader(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6) << "numiter";

  // different style due relative or absolute error checking
  // displacement
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs:
      oss << std::setw(18) << "abs-res";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::POROELAST::convnorm_abs :
    oss <<std::setw(18)<< "abs-inc";
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  switch ( normtypefres_ )
  {
    case INPAR::POROELAST::convnorm_abs :
    oss <<std::setw(18)<< "abs-s-res";
   // oss <<std::setw(18)<< "abs-f-res";
    oss <<std::setw(18)<< "abs-fvel-res";
    oss <<std::setw(18)<< "abs-fpres-res";
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::POROELAST::convnorm_abs :
    oss <<std::setw(18)<< "abs-s-inc";
   // oss <<std::setw(18)<< "abs-f-inc";
    oss <<std::setw(18)<< "abs-fvel-inc";
    oss <<std::setw(18)<< "abs-fpres-inc";
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  // add solution time
  oss << std::setw(14)<< "wct";

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile==NULL)
  dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;
}  // PrintNewtonIterHeader()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen                       |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::PrintNewtonIterText(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(7) << iter_;

  // different style due relative or absolute error checking
  // displacement
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << normrhs_;
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::POROELAST::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << norminc_;
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  switch ( normtypefres_ )
  {
    case INPAR::POROELAST::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normrhsstruct_;
  //  oss << std::setw(18) << std::setprecision(5) << std::scientific << normrhsfluid_;
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normrhsfluidvel_;
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normrhsfluidpres_;
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::POROELAST::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normincstruct_;
 //   oss << std::setw(18) << std::setprecision(5) << std::scientific << normincfluid_;
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normincfluidvel_;
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normincfluidpres_;
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  // add solution time
  oss << std::setw(14) << std::setprecision(2) << std::scientific << timer_.ElapsedTime();

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile==NULL)
  dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;

}  // PrintNewtonIterText

/*----------------------------------------------------------------------*
 | print statistics of converged NRI                                     |
 | orignially by bborn 08/09                                            |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::PrintNewtonConv()
{
  // somebody did the door
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate mechanical-fluid system matrix at state                    |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::ApplyStrCouplMatrix(
    Teuchos::RCP<
    LINALG::SparseOperator> k_sf //!< off-diagonal tangent matrix term
)
{
  k_sf->Zero();
  const Teuchos::ParameterList& sdynparams
  = DRT::Problem::Instance()->StructuralDynamicParams();

  // create the parameters for the discretization
  Teuchos::ParameterList sparams;

  switch (strmethodname_)
  {
  /*
   case  INPAR::STR::dyna_statics :
   {
   // continue
   break;
   }*/
  case INPAR::STR::dyna_onesteptheta:
  {
    double theta = sdynparams.sublist("ONESTEPTHETA").get<double> ("THETA");
    sparams.set("theta", theta);
    break;
  }
    // TODO: time factor for genalpha
    /*
     case INPAR::STR::dyna_genalpha :
     {
     double alphaf_ = sdynparams.sublist("GENALPHA").get<double>("ALPHA_F");
     // K_Teffdyn(T_n+1) = (1-alphaf_) . kst
     // Lin(dT_n+1-alphaf_/ dT_n+1) = (1-alphaf_)
     k_st->Scale(1.0 - alphaf_);
     }*/
  default:
  {
    dserror("Don't know what to do... only one-step theta time integration implemented");
    break;
  }
  } // end of switch(strmethodname_)

  const std::string action = "calc_struct_multidofsetcoupling";
  sparams.set("action", action);
  // other parameters that might be needed by the elements
  sparams.set("delta time", Dt());
  sparams.set("total time", Time());

  StructureField()->Discretization()->ClearState();
  StructureField()->Discretization()->SetState(0,"displacement",StructureField()->Dispnp());
  StructureField()->Discretization()->SetState(0,"velocity",StructureField()->ExtractVelnp());
  StructureField()->Discretization()->SetState(1,"fluidvel",FluidField().Velnp());

  // build specific assemble strategy for mechanical-fluid system matrix
  // from the point of view of StructureField:
  // structdofset = 0, fluiddofset = 1
  DRT::AssembleStrategy structuralstrategy(
      0,               // structdofset for row
      1,               // fluiddofset for column
      k_sf,            // mechanical-fluid coupling matrix
      Teuchos::null ,
      Teuchos::null ,
      Teuchos::null,
      Teuchos::null
  );

  // evaluate the mechancial-fluid system matrix on the structural element
  //StructureField()->Discretization()->Evaluate( sparams, structuralstrategy );
  StructureField()->Discretization()->EvaluateCondition( sparams, structuralstrategy,"PoroCoupling" );
  StructureField()->Discretization()->ClearState();

  return;
}    // ApplyStrCouplMatrix()

/*----------------------------------------------------------------------*
 |    evaluate fluid-structural system matrix at state                |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::ApplyFluidCouplMatrix(
    Teuchos::RCP< LINALG::SparseOperator> k_fs //!< off-diagonal tangent matrix term
  )
{
  k_fs->Zero();

  // create the parameters for the discretization
  Teuchos::ParameterList fparams;
  // action for elements
  fparams.set<int>("action", FLD::calc_porousflow_fluid_coupling);
  // other parameters that might be needed by the elements
  fparams.set("delta time", Dt());
  fparams.set("total time", Time());
  // create specific time integrator
  const Teuchos::ParameterList& fdyn =
      DRT::Problem::Instance()->FluidDynamicParams();
  fparams.set<int> ("time integrator", DRT::INPUT::IntegralValue<
      INPAR::FLUID::TimeIntegrationScheme>(fdyn, "TIMEINTEGR"));

  switch (DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn,
      "TIMEINTEGR"))
  {
  /*
   // Static analysis
   case INPAR::THR::dyna_statics :
   {
   break;
   }*/
  // Static analysis
  case INPAR::FLUID::timeint_one_step_theta:
  {
    double theta = fdyn.get<double> ("THETA");
    fparams.set("theta", theta);
    break;
  }
    /*case INPAR::THR::dyna_genalpha :
     {
     dserror("Genalpha not yet implemented");
     break;
     }*/
  default:
  {
    dserror("Don't know what to do...");
    break;
  }
  }

  FluidField().Discretization()->ClearState();

  // set general vector values needed by elements
  FluidField().Discretization()->SetState(0,"hist",FluidField().Hist());
  FluidField().Discretization()->SetState(0,"accam",FluidField().Accam());
  FluidField().Discretization()->SetState(0,"dispnp",FluidField().Dispnp());
  FluidField().Discretization()->SetState(0,"gridv",FluidField().GridVel());
  FluidField().Discretization()->SetState(0,"dispn",FluidField().Dispn());
  FluidField().Discretization()->SetState(0,"veln",FluidField().Veln());
  FluidField().Discretization()->SetState(0,"accnp",FluidField().Accnp());

  // set scheme-specific element parameters and vector values
  //TODO
  //if (is_genalpha_)
  //    discret_->SetState("velaf",velaf_);
  //else

  FluidField().Discretization()->SetState(0,"velaf",FluidField().Velnp());

  FluidField().Discretization()->SetState(0,"velnp",FluidField().Velnp());

  // build specific assemble strategy for the fluid-mechanical system matrix
  // from the point of view of FluidField:
  // fluiddofset = 0, structdofset = 1
  DRT::AssembleStrategy fluidstrategy(
      0,              // fluiddofset for row
      1,              // structdofset for column
      k_fs,           // fluid-mechanical matrix
      Teuchos::null,  // no other matrix or vectors
      Teuchos::null ,
      Teuchos::null,
      Teuchos::null
  );

  // evaluate the fluid-mechancial system matrix on the fluid element
  //FluidField().Discretization()->Evaluate(fparams,fluidstrategy);
  FluidField().Discretization()->EvaluateCondition( fparams, fluidstrategy,"PoroCoupling" );
  FluidField().Discretization()->ClearState();

  //apply normal flux condition on coupling part
  if(nopencond_.size())
  {
    k_fs->Complete(StructureField()->SystemMatrix()->RangeMap(), FluidField().SystemMatrix()->RangeMap());
    const Teuchos::RCP<const Epetra_Map >& nopenetrationmap = nopenetration_->Map(1);
    k_fs->ApplyDirichlet(*nopenetrationmap, false);

    cond_rhs_ = rcp(new Epetra_Vector(*DofRowMap(), true));
    cond_rhs_->PutScalar(0.0);

    EvaluateCondition(k_fs,cond_rhs_,1);
  }
}    // ApplyFluidCouplMatrix()

/*----------------------------------------------------------------------*
 |  map containing the dofs with Dirichlet BC               vuong 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> POROELAST::Monolithic::CombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map> scondmap =
      StructureField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> fcondmap =
      FluidField().GetDBCMapExtractor()->CondMap();
  Teuchos::RCP<Epetra_Map> condmap =
      LINALG::MergeMap(scondmap, fcondmap, false);
  return condmap;
} // CombinedDBCMap()

/*----------------------------------------------------------------------*
 |  check tangent stiffness matrix vie finite differences               |
 *----------------------------------------------------------------------*/
/*
void POROELAST::Monolithic::PoroFDCheck()
{
  cout << "\n******************finite difference check***************" << endl;

  int dof_struct = (StructureField()->Discretization()->NumGlobalNodes()) * 3;
  int dof_fluid = (FluidField().Discretization()->NumGlobalNodes()) * 4;

  cout << "structure field has " << dof_struct << " DOFs" << endl;
  cout << "fluid field has " << dof_fluid << " DOFs" << endl;

  Teuchos::RCP<Epetra_Vector> iterinc = Teuchos::null;
  iterinc = LINALG::CreateVector(*DofRowMap(), true);

  const int dofs = iterinc->GlobalLength();
  cout << "in total " << dofs << " DOFs" << endl;
  const double delta = 1e-8;

  iterinc->PutScalar(0.0);

  iterinc->ReplaceGlobalValue(0, 0, delta);

  Teuchos::RCP<Epetra_CrsMatrix> stiff_approx = Teuchos::null;
  stiff_approx = LINALG::CreateMatrix(*DofRowMap(), 81);

  //Teuchos::RCP<Epetra_Vector> rhs_old= null;
  Teuchos::RCP<Epetra_Vector> rhs_old = rcp(new Epetra_Vector(*DofRowMap(),
      true));
  rhs_old->Update(1.0, *rhs_, 0.0);
  Teuchos::RCP<Epetra_Vector> rhs_copy = rcp(new Epetra_Vector(*DofRowMap(),
      true));
  //rhs_old = rcp(new Epetra_Vector(*rhs_));
  //cout<<"rhs_"<<endl<<*rhs_<<endl;
  //cout<<"rhs_old"<<endl<<*rhs_old<<endl;

  Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();
  //cout<<"DBCMap:"<<endl<<*CombinedDBCMap();
  Teuchos::RCP<LINALG::SparseMatrix> sparse_copy = Teuchos::rcp(
      new LINALG::SparseMatrix(*(sparse->EpetraMatrix())));

  if (false)
  {
    cout << "iterinc_" << endl << *iterinc_ << endl;
    cout << "iterinc" << endl << *iterinc << endl;
    cout << "meshdisp: " << endl << *(FluidField().Dispnp());
    cout << "disp: " << endl << *(StructureField()->Dispnp());
    cout << "fluid vel" << endl << *(FluidField().Velnp());
    cout << "fluid acc" << endl << *(FluidField().Accnp());
    cout << "gridvel fluid" << endl << *(FluidField().GridVel());
    cout << "gridvel struct" << endl << *(StructureField()->ExtractVelnp());
  }

  int spaltenr = -1;
  int zeilennr = -1;
  for (int i = 0; i < dofs; ++i)
  {
    if (CombinedDBCMap()->MyGID(i))
    {
      iterinc->ReplaceGlobalValue(i, 0, 0.0);
    }

    if (i == spaltenr)
      cout << "\n******************" << spaltenr + 1
          << ". Spalte!!***************" << endl;

    // cout<<"iterinc anfang: "<<endl<< *iterinc<<endl;

    Evaluate(iterinc);
    SetupRHS();
    if(nopencond_.size())
    {
      const Teuchos::RCP<const Epetra_Map >& nopenetrationmap = nopenetration_.Map(1);
      LINALG::ApplyDirichlettoSystem(iterinc_,rhs_,cond_rhs_,*nopenetrationmap);
    }
    rhs_copy->Update(1.0, *rhs_, 0.0);

    iterinc_->PutScalar(0.0); // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(sparse_copy, iterinc_, rhs_copy,
        Teuchos::null, zeros_, *CombinedDBCMap());


    if (i == spaltenr)
    {

      cout << "rhs_: " << (*rhs_copy)[zeilennr] << endl;
      cout << "rhs_old: " << (*rhs_old)[zeilennr] << endl;
    }

    rhs_copy->Update(-1.0, *rhs_old, 1.0);
    rhs_copy->Scale(-1.0 / delta);

    int* index = &i;
    for (int j = 0; j < dofs; ++j)
    {
      double value = (*rhs_copy)[j];
      stiff_approx->InsertGlobalValues(j, 1, &value, index);

      if ((j == zeilennr) and (i == spaltenr))
      {
        cout << "\n******************" << zeilennr + 1
            << ". Zeile!!***************" << endl;
        cout << "iterinc_" << endl << *iterinc_ << endl;
        cout << "iterinc" << endl << *iterinc << endl;
        cout << "meshdisp: " << endl << *(FluidField().Dispnp());
        cout << "disp: " << endl << *(StructureField()->Dispnp());
        cout << "fluid vel" << endl << *(FluidField().Velnp());
        cout << "fluid acc" << endl << *(FluidField().Accnp());
        cout << "gridvel fluid" << endl << *(FluidField().GridVel());
        cout << "gridvel struct" << endl << *(StructureField()->ExtractVelnp());

        cout << "stiff_apprx(" << zeilennr << "," << spaltenr << "): "
            << (*rhs_copy)[zeilennr] << endl;

        cout << "stiff_apprx(" << zeilennr << "," << spaltenr << "): "
            << (*rhs_copy)[zeilennr] << endl;
        cout << "value(" << zeilennr << "," << spaltenr << "): " << value
            << endl;
        cout << "\n******************" << zeilennr + 1
            << ". Zeile Ende!!***************" << endl;
      }
    }

    //  cout<<"stiff_approx, column "<< i<<endl<<*stiff_approx<<endl;

    if (not CombinedDBCMap()->MyGID(i))
      iterinc->ReplaceGlobalValue(i, 0, -delta);

    iterinc->ReplaceGlobalValue(i - 1, 0, 0.0);

    if (i != dofs - 1)
      iterinc->ReplaceGlobalValue(i + 1, 0, delta);

    if (i == spaltenr)
      cout << "\n******************" << spaltenr + 1
          << ". Spalte Ende!!***************" << endl;

  }

  //cout<<"iterinc ende"<<endl<<*iterinc<<endl;
  Evaluate(iterinc);
  SetupRHS();

  //cout<<"vel ende"<<endl<<*(FluidField().Velnp());

  stiff_approx->FillComplete();

  //    cout<<"stiff_approx"<<endl<<*stiff_approx;
  //    cout<<"systemmatrix_"<<endl<<*systemmatrix_;
  //    cout<<"sparse"<<endl<<*sparse;

  Teuchos::RCP<LINALG::SparseMatrix> stiff_approx_sparse = Teuchos::null;
  stiff_approx_sparse = Teuchos::rcp(new LINALG::SparseMatrix(*stiff_approx));

  stiff_approx_sparse->Add(*sparse_copy, false, -1.0, 1.0);

  Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = sparse_copy->EpetraMatrix();

  Teuchos::RCP<Epetra_CrsMatrix> error_crs =
      stiff_approx_sparse->EpetraMatrix();

  error_crs->FillComplete();
  sparse_crs->FillComplete();

  //  cout<<"stiff_approx"<<endl<<*stiff_approx;
  //  cout<<"error_crs"<<endl<<*error_crs;
  // cout<<"sparse_crs"<<endl<<*sparse_crs;
 // cout << "DBCMap:" << endl << *CombinedDBCMap();

  bool success = true;
  double error_max = 0.0;
  for (int i = 0; i < dofs; ++i)
  {
    if (not CombinedDBCMap()->MyGID(i))
    {
      for (int j = 0; j < dofs; ++j)
      {
        if (not CombinedDBCMap()->MyGID(j))
        {
          double stiff_approx_ij = ((*stiff_approx)[i][j]);
          double sparse_ij = ((*sparse_crs)[i][j]);

          double error = 0.0;
          if (abs(stiff_approx_ij) > 1e-5)
            error = (*error_crs)[i][j] / (stiff_approx_ij);
          else if (abs(sparse_ij) > 1e-5)
            error = (*error_crs)[i][j] / (sparse_ij);

          if (abs(error) > abs(error_max))
            error_max = abs(error);

          if ((abs(error) > 1e-4))
          {
            if ((abs((*error_crs)[i][j]) > 1e-5))
            //  if( (sparse_ij>1e-1) or (stiff_approx_ij>1e-1) )
            {
              cout << "finite difference check failed entry (" << i << "," << j
                  << ")! stiff: " << sparse_ij << ", approx: "
                  << stiff_approx_ij << " ,abs. error: " << (*error_crs)[i][j]
                  << " , rel. error: " << error << endl;

              success = false;
            }
          }
        }
      }
    }
    // else cout<<"GID "<<i<<" mit DBC"<<endl;
  }

  if(success)
  {
    cout << "finite difference check successful, max. rel. error: "
        << error_max << endl;
    cout << "******************finite difference check done***************\n\n"
        << endl;
  }

  return;
}
*/

/*----------------------------------------------------------------------*
 |   evaluate poroelasticity specific constraint            vuong 03/12 |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::EvaluateCondition(Teuchos::RCP<LINALG::SparseOperator> Sysmat,
                                              Teuchos::RCP<Epetra_Vector> Cond_RHS,
                                              int coupltype)
{

  Teuchos::RCP<LINALG::SparseMatrix> ConstraintMatrix = Teuchos::null;
  ConstraintMatrix = Teuchos::rcp(new LINALG::SparseMatrix(
                        *(FluidField().Discretization()->DofRowMap(0)),
                        //*(FluidField().DofRowMap()),
                        StructureField()->Discretization()->DofRowMap()->NumGlobalElements(),
                        true, true));

  Teuchos::RCP<LINALG::SparseMatrix> StructVelConstraintMatrix = Teuchos::null;

  StructVelConstraintMatrix = Teuchos::rcp(new LINALG::SparseMatrix(
                        *(FluidField().Discretization()->DofRowMap(0)),
                        StructureField()->Discretization()->DofRowMap()->NumGlobalElements(),
                        true, true));

  ADAPTER::FluidPoro& fluidfield = dynamic_cast<ADAPTER::FluidPoro&>(FluidField());

  condIDs_ = rcp(new std::set<int>());
  condIDs_->clear();

  //evaluate condition on elements and assemble matrixes
  fluidfield.EvaluateNoPenetrationCond( Cond_RHS,
                                        ConstraintMatrix,
                                        StructVelConstraintMatrix,
                                        condIDs_,
                                        coupltype);

  if(coupltype==0)//fluid fluid part
  {
    ConstraintMatrix->Complete();
    BuidNoPenetrationMap();
  }
  else //fluid structure part
  {
    //double timescale = FluidField().TimeScaling();
    double timescale = FluidField().ResidualScaling();
    StructVelConstraintMatrix->Scale(timescale);
    StructVelConstraintMatrix->Complete(StructureField()->SystemMatrix()->RangeMap(), FluidField().SystemMatrix()->RangeMap());
    ConstraintMatrix->Add(*StructVelConstraintMatrix, false, 1.0, 1.0);
    ConstraintMatrix->Complete(StructureField()->SystemMatrix()->RangeMap(), FluidField().SystemMatrix()->RangeMap());
  }

  const Teuchos::RCP<const Epetra_Map >& nopenetrationmap = nopenetration_->Map(1);
  Sysmat->ApplyDirichlet(*nopenetrationmap, false);
  Sysmat->UnComplete();
  Sysmat->Add(*ConstraintMatrix, false, 1.0, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> POROELAST::Monolithic::StructFluidCouplingMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_sf_);
} // StructFluidCouplingMatrix()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> POROELAST::Monolithic::FluidStructCouplingMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_fs_);
} // FluidStructCouplingMatrix()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> POROELAST::Monolithic::StructFluidCouplingBlockMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(k_sf_);
} // StructFluidCouplingBlockMatrix()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> POROELAST::Monolithic::FluidStructCouplingBlockMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(k_fs_);
} // FluidStructCouplingBlockMatrix()



/*----------------------------------------------------------------------*
 | Setup Newton-Raphson iteration            vuong 01/12   |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::SetupNewton()
{

  // initialise equilibrium loop and norms
  iter_ = 1;
  normrhs_ = 0.0;
  norminc_ = 0.0;
  normrhsfluid_ = 0.0;
  normincfluid_ = 0.0;
  normrhsfluidvel_ = 0.0;
  normincfluidvel_ = 0.0;
  normrhsstruct_ = 0.0;
  normincstruct_ = 0.0;

  // incremental solution vector with length of all dofs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  iterinc_->PutScalar(0.0);

  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  zeros_->PutScalar(0.0);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::BuildCovergenceNorms()
{
  // build residual force norm
  // for now use for simplicity only L2/Euclidian norm
  rhs_->Norm2(&normrhs_);
  Teuchos::RCP<const Epetra_Vector> rhs_s;
  Teuchos::RCP<const Epetra_Vector> rhs_f;
  Teuchos::RCP<const Epetra_Vector> rhs_fvel;
  Teuchos::RCP<const Epetra_Vector> rhs_fpres;

  // process structure unknowns of the first field
  rhs_s = Extractor().ExtractVector(rhs_, 0);
  // process fluid unknowns of the second field
  rhs_f = Extractor().ExtractVector(rhs_, 1);
  rhs_fvel = FluidField().ExtractVelocityPart(rhs_f);
  rhs_fpres = FluidField().ExtractPressurePart(rhs_f);

  rhs_s->Norm2(&normrhsstruct_);
  rhs_f->Norm2(&normrhsfluid_);
  rhs_fvel->Norm2(&normrhsfluidvel_);
  rhs_fpres->Norm2(&normrhsfluidpres_);

  // build residual increment norm
  iterinc_->Norm2(&norminc_);

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> interincs;
  Teuchos::RCP<const Epetra_Vector> interincf;
  Teuchos::RCP<const Epetra_Vector> interincfvel;
  Teuchos::RCP<const Epetra_Vector> interincfpres;
  // process structure unknowns of the first field
  interincs = Extractor().ExtractVector(iterinc_, 0);
  // process fluid unknowns of the second field
  interincf = Extractor().ExtractVector(iterinc_, 1);
  interincfvel = FluidField().ExtractVelocityPart(interincf);
  interincfpres = FluidField().ExtractPressurePart(interincf);

  interincs->Norm2(&normincstruct_);
  interincf->Norm2(&normincfluid_);
  interincfvel->Norm2(&normincfluidvel_);
  interincfpres->Norm2(&normincfluidpres_);

  return;
}
