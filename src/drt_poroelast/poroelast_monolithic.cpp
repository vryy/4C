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

#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
#endif

/*----------------------------------------------------------------------*
 | headers                                                  vuong 01/12 |
 *----------------------------------------------------------------------*/
#include "poroelast_monolithic.H"
#include "poroelast_defines.H"

#include "../drt_adapter/adapter_coupling.H"
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
#include "../drt_adapter/ad_fld_poro.H"

//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB genprob;

/*----------------------------------------------------------------------*
 | constructor (public)                                    vuong 01/12  |
 *----------------------------------------------------------------------*/
POROELAST::MonolithicBase::MonolithicBase(const Epetra_Comm& comm) :
      AlgorithmBase(comm, DRT::Problem::Instance()->PoroelastDynamicParams()),
      StructureBaseAlgorithm(DRT::Problem::Instance()->PoroelastDynamicParams()),
      FluidBaseAlgorithm(DRT::Problem::Instance()->PoroelastDynamicParams(),
      true)
{
  // monolithic Poroelasticity must know the other discretization
  // build a proxy of the structure discretization for the fluid field
  Teuchos::RCP<DRT::DofSet> structdofset =
      StructureField().Discretization()->GetDofSetProxy();
  // build a proxy of the fluid discretization for the structure field
  Teuchos::RCP<DRT::DofSet> fluiddofset =
      FluidField().Discretization()->GetDofSetProxy();

  // check if FluidField has 2 discretizations, so that coupling is possible
  if (FluidField().Discretization()->AddDofSet(structdofset) != 1)
    dserror("unexpected dof sets in fluid field");
  if (StructureField().Discretization()->AddDofSet(fluiddofset)!=1)
    dserror("unexpected dof sets in structure field");

  // access the problem-specific parameter lists
  const Teuchos::ParameterList& sdyn
  = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& fdyn
  = DRT::Problem::Instance()->FluidDynamicParams();

  // check time integration algo -> currently only one-step-theta scheme supported
  INPAR::STR::DynamicType structtimealgo
  = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP");
  INPAR::FLUID::TimeIntegrationScheme fluidtimealgo
  = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn,"TIMEINTEGR");

  if ( structtimealgo != INPAR::STR::dyna_onesteptheta or
      fluidtimealgo != INPAR::FLUID::timeint_one_step_theta )
  dserror("monolithic Poroelasticity is limited in functionality (only one-step-theta scheme possible)");

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* structurenodemap = StructureField().Discretization()->NodeRowMap();

  coupfa_ = Teuchos::rcp(new ADAPTER::Coupling());
  const int ndim = DRT::Problem::Instance()->NDim();
  coupfa_->SetupCoupling(*FluidField().Discretization(),
      *StructureField().Discretization(),
      *fluidnodemap,
      *structurenodemap,
      ndim);

  FluidField().SetMeshMap(coupfa_->MasterDofMap());

  //extractor for constraints on structure phase
  //
  // when using constraints applied via Lagrange-Multipliers there is a
  // difference between StructureField().DofRowMap() and StructureField().DofRowMap(0).
  // StructureField().DofRowMap(0) returns the DofRowMap
  // known to the discretization (without lagrange multipliers)
  // while StructureField().DofRowMap() returns the DofRowMap known to
  // the constraint manager (with lagrange multipliers)
  consplitter_=LINALG::MapExtractor(*StructureField().DofRowMap(),
      StructureField().DofRowMap(0));

  FluidField().Discretization()->GetCondition("NoPenetration", nopencond_);
}

/*----------------------------------------------------------------------*
 | destructor (public)                                    vuong 01/12   |
 *----------------------------------------------------------------------*/
POROELAST::MonolithicBase::~MonolithicBase()
{
}

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 01/12  |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicBase::ReadRestart(int step)
{
  FluidField().ReadRestart(step);
  StructureField().ReadRestart(step);

  // apply current velocity and pressures to structure
  StructureField().ApplyVelAndPress(FluidField().Velnp());

  Teuchos::RCP<Epetra_Vector> dispn;
  if (StructureField().HaveConstraint())
  {
    //displacment vector without lagrange-multipliers
    dispn = consplitter_.ExtractCondVector(StructureField().Dispnp());
  }
  else
    dispn = StructureField().ExtractDispnp();

  // transfer the current structure displacement to the fluid field
  Teuchos::RCP<Epetra_Vector> structdisp = StructureToFluidField(dispn);
  FluidField().ApplyMeshDisplacement(structdisp);

  // transfer the current structure velocity to the fluid field
  Teuchos::RCP<Epetra_Vector> structvel = StructureToFluidField(
      StructureField().ExtractVelnp());
  FluidField().ApplyMeshVelocity(structvel);

  // second ReadRestart needed due to the coupling variables
  FluidField().ReadRestart(step);
  StructureField().ReadRestart(step);

  SetTimeStep(FluidField().Time(), step);

  return;
}

/*----------------------------------------------------------------------*
 | prepare time step (public)                         vuong 01/12       |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicBase::PrepareTimeStep()
{
  // counter and print header
  IncrementTimeAndStep();
  PrintHeader();

  // call the predictor
  StructureField().PrepareTimeStep();
  FluidField().PrepareTimeStep();
}

/*----------------------------------------------------------------------*
 | update (protected)                                     vuong 01/12   |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicBase::Update()
{
  StructureField().Update();
  FluidField().Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> POROELAST::MonolithicBase::StructureToFluidField(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> POROELAST::MonolithicBase::FluidToStructureField(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::MonolithicBase::BuidNoPenetrationMap()
{
  std::vector<int> condIDs;
  std::set<int>::iterator it;
  for(it=condIDs_->begin();it!=condIDs_->end();it++)
  {
    condIDs.push_back(*it);
  }
  Teuchos::RCP<Epetra_Map> nopendofmap = rcp(new Epetra_Map(-1, condIDs.size(), &condIDs[0], 0, FluidField().Discretization()->Comm()));

  nopenetration_ = LINALG::MapExtractor(*FluidField().DofRowMap(), nopendofmap);

  return;
}

/*----------------------------------------------------------------------*
 | monolithic                                              vuong 01/12  |
 *----------------------------------------------------------------------*/
POROELAST::Monolithic::Monolithic(const Epetra_Comm& comm,
    const Teuchos::ParameterList& sdynparams) :
    MonolithicBase(comm),
    solveradapttol_(DRT::INPUT::IntegralValue<int>(sdynparams, "ADAPTCONV") == 1),
    solveradaptolbetter_(sdynparams.get<double> ("ADAPTCONV_BETTER")),
    printscreen_(true), // ADD INPUT PARAMETER
    printiter_(true), // ADD INPUT PARAMETER
    printerrfile_(true and errfile_), // ADD INPUT PARAMETER FOR 'true'
    errfile_(NULL),
    zeros_(Teuchos::null),
    strmethodname_(DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams,"DYNAMICTYP")),
    timer_(comm),
    veln_(Teuchos::null),
    dispn_(Teuchos::null)
{
  // add extra parameters (a kind of work-around)
  Teuchos::RCP<Teuchos::ParameterList> xparams = Teuchos::rcp(
      new Teuchos::ParameterList());
  xparams->set<FILE*> ("err file",
      DRT::Problem::Instance()->ErrorFile()->Handle());
  errfile_ = xparams->get<FILE*> ("err file");

  // velocities V_{n+1} at t_{n+1}
  veln_ = LINALG::CreateVector(*(StructureField().DofRowMap()), true);
  veln_->PutScalar(0.0);

  // displacements V_{n+1} at t_{n+1}
  dispn_ = LINALG::CreateVector(*(StructureField().DofRowMap()), true);
  dispn_->PutScalar(0.0);

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

    solver_ = rcp(new LINALG::Solver(
            solverparams,
            Comm(),
            DRT::Problem::Instance()->ErrorFile()->Handle()
        )
    );

#else
        dserror("implicit solver not implemented");
#endif

}

/*----------------------------------------------------------------------*
 | output                                                vuong 01/12    |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  FluidField().Output();
  StructureField().Output();
} // Monolithic::Output()

/*----------------------------------------------------------------------*
 | time loop of the monolithic system                    vuong 01/12     |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::TimeLoop(const Teuchos::ParameterList& sdynparams)
{
  // time loop
  while (NotFinished())
  {
    // counter and print header
    // predict solution of both field (call the adapter)
    PrepareTimeStep();

    // Newton-Raphson iteration
    NewtonFull(sdynparams);

    // calculate stresses, strains, energies
    PrepareOutput();

    // update all single field solvers
    Update();

    // write output to screen and files
    Output();

  } // NotFinished
} // TimeLoop

/*----------------------------------------------------------------------*
 | solution with full Newton-Raphson iteration            vuong 01/12   |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::NewtonFull(const Teuchos::ParameterList& sdynparams)
{

  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #rhs_ is the positive force residuum
  // --> On #systemmatrix_ is the effective dynamic tangent matrix

  // time parameters
  // call the Poroelasticity parameter list
  const Teuchos::ParameterList& poroelastdyn =
      DRT::Problem::Instance()->PoroelastDynamicParams();
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

  // initialise equilibrium loop
  iter_ = 1;
  normrhs_ = 0.0;
  norminc_ = 0.0;
  normrhsfluid_ = 0.0;
  normincfluid_ = 0.0;
  normrhsstruct_ = 0.0;
  normincstruct_ = 0.0;
  Epetra_Time timerporoelast(Comm());
  timerporoelast.ResetStartTime();

  // incremental solution vector with length of all Poroelasticity dofs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  // iterinc_ = LINALG::CreateVector(*(StructureField().DofRowMap(0)), true);
  iterinc_->PutScalar(0.0);

  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  //zeros_ = LINALG::CreateVector(*(StructureField().DofRowMap(0)), true);
  zeros_->PutScalar(0.0);

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
    SetupSystemMatrix(sdynparams);

    // check whether we have a sanely filled tangent matrix
    if (not systemmatrix_->Filled())
    {
      dserror("Effective tangent matrix must be filled here");
    }

    // create full monolithic rhs vector
    SetupRHS();

    //#ifdef POROFDCHECK
    // PoroFDCheck();
    //#endif

    // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
    // is done in PrepareSystemForNewtonSolve() within Evaluate(iterinc_)
    LinearSolve();

    // reset solver tolerance
    solver_->ResetTolerance();

    // build residual force norm
    // for now use for simplicity only L2/Euclidian norm
    rhs_->Norm2(&normrhs_);
    Teuchos::RCP<const Epetra_Vector> rhs_s;
    Teuchos::RCP<const Epetra_Vector> rhs_f;
    ExtractFieldVectors(rhs_,rhs_s,rhs_f);
    rhs_s->Norm2(&normrhsstruct_);
    rhs_f->Norm2(&normrhsfluid_);

    // build residual increment norm
    iterinc_->Norm2(&norminc_);

    // displacement and fluid velocity & pressure incremental vector
    Teuchos::RCP<const Epetra_Vector> interincs;
    Teuchos::RCP<const Epetra_Vector> interincf;
    ExtractFieldVectors(iterinc_,interincs,interincf);

    interincs->Norm2(&normincstruct_);
    interincf->Norm2(&normincfluid_);

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
    // extract displacement sx and temperature tx incremental vector of global
    // unknown incremental vector x
    ExtractFieldVectors(x, sx, fx);
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
  StructureField().ApplyVelAndPress(FluidField().Velnp());

  // Monolithic Poroelasticity accesses the linearised structure problem:
  //   UpdaterIterIncrementally(sx),
  //   EvaluateForceStiffResidual()
  //   PrepareSystemForNewtonSolve()
  StructureField().Evaluate(sx);
  //cout << "  structure time for calling Evaluate: " << timerstructure.ElapsedTime() << "\n";

  /// fluid field

  // fluid Evaluate
  // (builds tangent, residual and applies DBC and recent coupling values)
  //Epetra_Time timerfluid(Comm());

  /*
   // apply current displacements and velocities to the fluid field
   if (strmethodname_==INPAR::STR::dyna_statics)
   {
   // calculate velocity V_n+1^k = (D_n+1^k-D_n)/Dt()
   veln_ = CalcVelocity(StructureField().Dispnp());
   }
   else
   {
   veln_ = StructureField().ExtractVelnp();
   }*/

  if (StructureField().HaveConstraint())
  {
    //displacement vector without lagrange-multipliers
    dispn_ = consplitter_.ExtractCondVector(StructureField().Dispnp());
  }
  else
  {
    dispn_ = StructureField().ExtractDispnp();
  }

  veln_ = StructureField().ExtractVelnp();

  // transfer the current structure displacement to the fluid field
  Teuchos::RCP<Epetra_Vector> structdisp = StructureToFluidField(dispn_);
  FluidField().ApplyMeshDisplacement(structdisp);

  // transfer the current structure velocity to the fluid field
  Teuchos::RCP<Epetra_Vector> structvel = StructureToFluidField(veln_);
  FluidField().ApplyMeshVelocity(structvel);

  // monolithic Poroelasticity accesses the linearised thermo problem
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
    Teuchos::RCP<const Epetra_Vector>& fx)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::ExtractFieldVectors");

  // process structure unknowns of the first field
  sx = Extractor().ExtractVector(x, 0);

  // process thermo unknowns of the second field
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
  vel = rcp(new Epetra_Vector(*(StructureField().ExtractDispn())));
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
  // difference between StructureField().DofRowMap() and StructureField().DofRowMap(0).
  // StructureField().DofRowMap(0) returns the DofRowMap
  // known to the discretization (without lagrange multipliers)
  // while StructureField().DofRowMap() returns the DofRowMap known to
  // the constraint manager (with lagrange multipliers)
  // In the constrained case we want the "whole" RowDofMap,
  // otherwise both calls are equivalent

  vecSpaces.push_back(StructureField().DofRowMap());
  vecSpaces.push_back(FluidField().DofRowMap(0));

  if (vecSpaces[0]->NumGlobalElements() == 0)
    dserror("No structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements()==0)
    dserror("No fluid equation. Panic.");

    SetDofRowMaps(vecSpaces);

  } // SetupSystem()

/*----------------------------------------------------------------------*
 | put the single maps to one full                                      |
 | Poroelasticity map together                              vuong 01/12 |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::SetDofRowMaps(const std::vector<Teuchos::RCP<
    const Epetra_Map> >& maps)
{
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);

  // full Poroelasticity-blockmap
  blockrowdofmap_.Setup(*fullmap, maps);
}

/*----------------------------------------------------------------------*
 | setup system matrix of poroelasticity                   vuong 01/12  |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::SetupSystemMatrix(
    const Teuchos::ParameterList& sdynparams)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::SetupSystemMatrix");

  /*----------------------------------------------------------------------*/
  // initialize Poroelasticity-systemmatrix_
  systemmatrix_ = rcp(new LINALG::BlockSparseMatrix<
      LINALG::DefaultBlockMatrixStrategy>(Extractor(), Extractor(), 81, false,
      true));

  /*----------------------------------------------------------------------*/
  // pure structural part k_ss (3nx3n)

  // build pure structural block k_ss
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  Teuchos::RCP<LINALG::SparseMatrix> k_ss = StructureField().SystemMatrix();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  k_ss->UnComplete();

  // assign structure part to the TSI matrix
  systemmatrix_->Assign(0, 0, View, *k_ss);

  /*----------------------------------------------------------------------*/
  // structural part k_sf (3nxn)
  // build mechanical-fluid block

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> k_sf = Teuchos::null;
  k_sf = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(StructureField().DofRowMap()), 81, true, true));

  // call the element and calculate the matrix block
  ApplyStrCouplMatrix(k_sf, sdynparams);

  // Uncomplete mechanical-fluid matrix to be able to deal with slightly
  // defective interface meshes.
  k_sf->UnComplete();

  // assign fluid part to the Poroelasticity matrix
  systemmatrix_->Assign(0, 1, View, *(k_sf));

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
  systemmatrix_->Assign(1, 1, View, *(k_ff));

  /*----------------------------------------------------------------------*/
  // fluid part k_fs ( (3n+1)x3n )
  // build fluid-mechanical block

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> k_fs = Teuchos::null;
  k_fs = Teuchos::rcp(new LINALG::SparseMatrix(
                        *(FluidField().Discretization()->DofRowMap(0)),
                        //*(FluidField().DofRowMap()),
                        81, true, true));

  // call the element and calculate the matrix block
  ApplyFluidCouplMatrix(k_fs, sdynparams);

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  k_fs->UnComplete();

  // assign fluid part to the Poroelasticity matrix
  systemmatrix_->Assign(1, 0, View, *(k_fs));

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  systemmatrix_->Complete();
} // SetupSystemMatrix

/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                                 vuong 01/12  |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::SetupRHS()
{
  //  cout << " POROELAST::Monolithic::SetupRHS()" << endl;
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::SetupRHS");

  // create full monolithic rhs vector
  rhs_ = rcp(new Epetra_Vector(*DofRowMap(), true));

  // fill the Poroelasticity rhs vector rhs_ with the single field rhss
  SetupVector(*rhs_, StructureField().RHS(), FluidField().RHS());

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
  LINALG::ApplyDirichlettoSystem(sparse, iterinc_, rhs_, Teuchos::null, zeros_,
      *CombinedDBCMap());
  //  if ( Comm().MyPID()==0 ) { cout << " DBC applied to system" << endl; }

  if(nopencond_.size())
  {
    const Teuchos::RCP<const Epetra_Map >& nopenetrationmap = nopenetration_.Map(1);
    //Teuchos::RCP<Epetra_Vector> iterinc = Teuchos::null;
    //iterinc = LINALG::CreateVector(*DofRowMap(), true);
    LINALG::ApplyDirichlettoSystem(iterinc_,rhs_,cond_rhs_,*nopenetrationmap);
  }

  // standard solver call
  solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, true, iter_ == 1);
  //  if ( Comm().MyPID()==0 ) { cout << " Solved" << endl; }

#else // use bgs2x2_operator
  dserror("implicit solver not implemented");

#endif
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
      StructureField().InitialGuess(),
      // returns residual temperatures or iterative thermal increment - tempi_
      FluidField().InitialGuess());
} // InitialGuess()

/*----------------------------------------------------------------------*
 | setup vector of the structure and fluid field            vuong 01/12|
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::SetupVector(Epetra_Vector &f, Teuchos::RCP<
    const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv)
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
    }

    // residual forces
    switch (normtypefres_)
    {
      case INPAR::POROELAST::convnorm_abs:
      convfres = normrhs_ < tolfres_;
      break;
      default:
      dserror("Cannot check for convergence of residual forces!");
    }

    // combine temperature-like and force-like residuals
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
      oss << std::setw(18) << "abs-res-norm";
      break;
    default:
      dserror("You should not turn up here.");
    }

    switch ( normtypeinc_ )
    {
      case INPAR::POROELAST::convnorm_abs :
      oss <<std::setw(18)<< "abs-inc-norm";
      break;
      default:
      dserror("You should not turn up here.");
    }

    switch ( normtypefres_ )
    {
      case INPAR::POROELAST::convnorm_abs :
      oss <<std::setw(18)<< "abs-s-res-norm";
      oss <<std::setw(18)<< "abs-f-res-norm";
      break;
      default:
      dserror("You should not turn up here.");
    }

    switch ( normtypeinc_ )
    {
      case INPAR::POROELAST::convnorm_abs :
      oss <<std::setw(18)<< "abs-s-inc-norm";
      oss <<std::setw(18)<< "abs-f-inc-norm";
      break;
      default:
      dserror("You should not turn up here.");
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
 | print Newton-Raphson iteration to screen                  			|
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
    }

    switch ( normtypeinc_ )
    {
      case INPAR::POROELAST::convnorm_abs :
      oss << std::setw(18) << std::setprecision(5) << std::scientific << norminc_;
      break;
      default:
      dserror("You should not turn up here.");
    }

    switch ( normtypefres_ )
    {
      case INPAR::POROELAST::convnorm_abs :
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normrhsstruct_;
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normrhsfluid_;
      break;
      default:
      dserror("You should not turn up here.");
    }

    switch ( normtypeinc_ )
    {
      case INPAR::POROELAST::convnorm_abs :
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normincstruct_;
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normincfluid_;
      break;
      default:
      dserror("You should not turn up here.");
    }

    //Epetra_Time timerporoelast(Comm());
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
 | print statistics of converged NRI                          			|
 | orignially by bborn 08/09                                            |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::PrintNewtonConv()
{
  // somebody did the door
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate mechanical-fluid system matrix at state      				 |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::ApplyStrCouplMatrix(Teuchos::RCP<
    LINALG::SparseMatrix> k_sf, //!< off-diagonal tangent matrix term
    const Teuchos::ParameterList& sdynparams)
{

  // if ( Comm().MyPID()==0 )
  //  cout << " POROELAST::Monolithic::ApplyStrCouplMatrix()" << endl;

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

  const Teuchos::ParameterList& porodyn
  = DRT::Problem::Instance()->PoroelastDynamicParams();
  sparams.set<double>("initporosity", porodyn.get<double>("INITPOROSITY"));

  StructureField().Discretization()->ClearState();
  StructureField().Discretization()->SetState(0,"displacement",dispn_);
  StructureField().Discretization()->SetState(0,"velocity",veln_);
  StructureField().Discretization()->SetState(1,"fluidvel",FluidField().Velnp());

  // build specific assemble strategy for mechanical-thermal system matrix
  // from the point of view of StructureField:
  // structdofset = 0, fluiddofset = 1
  DRT::AssembleStrategy structuralstrategy(
      0, // structdofset for row
      1, // fluiddofset for column
      k_sf, // build mechanical-fluid coupling matrix (static part)
      Teuchos::null ,//k_sf_rea,  // build mechanical-fluid coupling matrix (transient part)
      Teuchos::null ,//rhs_sf_,
      Teuchos::null, //rhs_sf_rea,
      Teuchos::null
  );

  // evaluate the mechancial-thermal system matrix on the structural element
  StructureField().Discretization()->Evaluate( sparams, structuralstrategy );
  StructureField().Discretization()->ClearState();

  return;
}    // ApplyStrCouplMatrix()

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
/***********************************************************************/
void POROELAST::Monolithic::ApplyFluidCouplMatrix(Teuchos::RCP<
    LINALG::SparseMatrix> k_fs, //!< off-diagonal tangent matrix term
    const Teuchos::ParameterList& sdynparams)
{
  // create the parameters for the discretization
  Teuchos::ParameterList fparams;
  // action for elements
  const std::string action = "calc_porousflow_fluid_coupling";
  fparams.set("action", action);
  // other parameters that might be needed by the elements
  fparams.set("delta time", Dt());
  fparams.set("total time", Time());
  // create specific time integrator
  const Teuchos::ParameterList& fdyn =
      DRT::Problem::Instance()->FluidDynamicParams();
  fparams.set<int> ("time integrator", DRT::INPUT::IntegralValue<
      INPAR::FLUID::TimeIntegrationScheme>(fdyn, "TIMEINTEGR"));
  const Teuchos::ParameterList& porodyn =
      DRT::Problem::Instance()->PoroelastDynamicParams();
  fparams.set<double> ("initporosity", porodyn.get<double> ("INITPOROSITY"));

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
      0, // fluiddofset for row
      1, // structdofset for column
      k_fs, // fluid-mechancial matrix
      Teuchos::null, // no other matrix or vectors
      Teuchos::null ,//rhs_fs_,
      Teuchos::null,
      Teuchos::null
  );

  // evaluate the fluid-mechancial system matrix on the fluid element
  FluidField().Discretization()->Evaluate(fparams,fluidstrategy);
  FluidField().Discretization()->ClearState();

  //apply normal flux condition on coupling part
  if(nopencond_.size())
  {
    k_fs->Complete(StructureField().SystemMatrix()->RangeMap(), FluidField().SystemMatrix()->RangeMap());
    const Teuchos::RCP<const Epetra_Map >& nopenetrationmap = nopenetration_.Map(1);
    k_fs->ApplyDirichlet(*nopenetrationmap, false);

    cond_rhs_ = rcp(new Epetra_Vector(*DofRowMap(), true));
    cond_rhs_->PutScalar(0.0);

    EvaluateCondition(k_fs,cond_rhs_,1);
    //cond_rhs_->Scale(-1.0);
  }
}    // ApplyFluidCouplMatrix()

/*----------------------------------------------------------------------*
 |  map containing the dofs with Dirichlet BC               vuong 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> POROELAST::Monolithic::CombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map> scondmap =
      StructureField().GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> fcondmap =
      FluidField().GetDBCMapExtractor()->CondMap();
  Teuchos::RCP<Epetra_Map> condmap =
      LINALG::MergeMap(scondmap, fcondmap, false);
  return condmap;

} // CombinedDBCMap()

/*----------------------------------------------------------------------*
 |  check tangent stiffness matrix vie finite differences               |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::PoroFDCheck()
{
  cout << "\n******************finite difference check***************" << endl;

  int dof_struct = (StructureField().Discretization()->NumGlobalNodes()) * 3;
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

  if (true)
  {
    cout << "iterinc_" << endl << *iterinc_ << endl;
    cout << "iterinc" << endl << *iterinc << endl;
    cout << "meshdisp: " << endl << *(FluidField().Dispnp());
    cout << "disp: " << endl << *(StructureField().Dispnp());
    cout << "fluid vel" << endl << *(FluidField().Velnp());
    cout << "fluid acc" << endl << *(FluidField().Accnp());
    cout << "gridvel fluid" << endl << *(FluidField().GridVel());
    cout << "gridvel struct" << endl << *(StructureField().ExtractVelnp());
  }

  int spaltenr = 1;
  int zeilennr = 48;
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
        cout << "disp: " << endl << *(StructureField().Dispnp());
        cout << "fluid vel" << endl << *(FluidField().Velnp());
        cout << "fluid acc" << endl << *(FluidField().Accnp());
        cout << "gridvel fluid" << endl << *(FluidField().GridVel());
        cout << "gridvel struct" << endl << *(StructureField().ExtractVelnp());

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

  /*
   iterinc_->PutScalar(0.0);  // Useful? depends on solver and more
   LINALG::ApplyDirichlettoSystem(
   sparse,
   iterinc_,
   rhs_,
   Teuchos::null,
   zeros_,
   *CombinedDBCMap()
   );*/

  //cout<<"vel ende"<<endl<<*(FluidField().Velnp());

  stiff_approx->FillComplete();

  //	  cout<<"stiff_approx"<<endl<<*stiff_approx;
  //	  cout<<"systemmatrix_"<<endl<<*systemmatrix_;
  //	  cout<<"sparse"<<endl<<*sparse;

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
  cout << "DBCMap:" << endl << *CombinedDBCMap();

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

  // if(success)
  {
    cout << "finite difference check successful, max. rel. error: "
        << error_max << endl;
    cout << "******************finite difference check done***************\n\n"
        << endl;
  }
  /*
   else
   {
   cout<<"stiff_approx: "<<endl<<(*stiff_approx)<<endl;
   cout<<"sparse_crs: "<<endl<<(*sparse_crs)<<endl;
   dserror("finite difference check failed!");
   }*/
  return;
}

/*----------------------------------------------------------------------*
 |   evaluate poroelasticity specific constraint            vuong 03/12 |
 *----------------------------------------------------------------------*/
void POROELAST::Monolithic::EvaluateCondition(Teuchos::RCP<LINALG::SparseMatrix> Sysmat,
                                              Teuchos::RCP<Epetra_Vector> Cond_RHS,
                                              int coupltype)
{

  Teuchos::RCP<LINALG::SparseMatrix> ConstraintMatrix = Teuchos::null;
  ConstraintMatrix = Teuchos::rcp(new LINALG::SparseMatrix(
                        *(FluidField().Discretization()->DofRowMap(0)),
                        //*(FluidField().DofRowMap()),
                        StructureField().Discretization()->DofRowMap()->NumGlobalElements(),
                        true, true));

  Teuchos::RCP<LINALG::SparseMatrix> StructVelConstraintMatrix = Teuchos::null;

  StructVelConstraintMatrix = Teuchos::rcp(new LINALG::SparseMatrix(
                        *(FluidField().Discretization()->DofRowMap(0)),
                        StructureField().Discretization()->DofRowMap()->NumGlobalElements(),
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

  if(coupltype==0)
  {
    ConstraintMatrix->Complete();
    BuidNoPenetrationMap();
  }
  else
  {
    const Teuchos::ParameterList& fdyn =
        DRT::Problem::Instance()->FluidDynamicParams();
    double timefactor = -1.0;
    double dt=Dt();

    switch (DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn,
        "TIMEINTEGR"))
    {
      case INPAR::FLUID::timeint_one_step_theta:
      {
        timefactor = fdyn.get<double> ("THETA");
        break;
      }
      default:
      {
        dserror("No penetration condition implemented for one-step-theta time integration scheme only");
        break;
      }//TODO other time integration schemes
    }

    StructVelConstraintMatrix->Scale(-1/(dt*timefactor));
    StructVelConstraintMatrix->Complete(StructureField().SystemMatrix()->RangeMap(), FluidField().SystemMatrix()->RangeMap());
    ConstraintMatrix->Add(*StructVelConstraintMatrix, false, 1.0, 1.0);
    ConstraintMatrix->Complete(StructureField().SystemMatrix()->RangeMap(), FluidField().SystemMatrix()->RangeMap());
  }

  const Teuchos::RCP<const Epetra_Map >& nopenetrationmap = nopenetration_.Map(1);
  Sysmat->ApplyDirichlet(*nopenetrationmap, false);
  Sysmat->UnComplete();
  Sysmat->Add(*ConstraintMatrix, false, 1.0, 1.0);

  return;
}

#endif  // CCADISCRET
