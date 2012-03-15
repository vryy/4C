
#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "../drt_lib/standardtypes_cpp.H"

#include "fsi_debugwriter.H"
#include "fsi_monolithic_nonox.H"
#include "fsi_statustest.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_colors.H"

#include "../drt_io/io_control.H"
#include "../drt_structure/stru_aux.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

/*----------------------------------------------------------------------*/
// constructor (public)
/*----------------------------------------------------------------------*/
FSI::MonolithicNoNOX::MonolithicNoNOX(const Epetra_Comm& comm,
                            const Teuchos::ParameterList& timeparams)
  : MonolithicBase(comm,timeparams),
    cout0_(StructureField().Discretization()->Comm(), std::cout),
    zeros_(Teuchos::null)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  // enable debugging
  if (DRT::INPUT::IntegralValue<int>(fsidyn,"DEBUGOUTPUT")==1)
  {
    sdbg_ = Teuchos::rcp(new UTILS::DebugWriter(StructureField().Discretization()));
    //fdbg_ = Teuchos::rcp(new UTILS::DebugWriter(FluidField().Discretization()));
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

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::Timeloop()
{
  while (NotFinished())
  {
    PrepareTimeStep();
    Newton();
    PrepareOutput();
    Update();
    Output();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::Newton()
{
  // initialise equilibrium loop
  iter_ = 1;
  normrhs_ = 0.0;
  normstrrhs_ = 0.0;
  normflrhs_ = 0.0;
  normalerhs_ = 0.0;
  norminc_ = 0.0;

  // get length of the structural, fluid and ale vector
  ns_ = (*(StructureField().RHS())).GlobalLength();
  nf_ = (*(FluidField().RHS())).GlobalLength();
  na_ = (*(AleField().RHS())).GlobalLength();

  x_sum_ = LINALG::CreateVector(*DofRowMap(),true);
  x_sum_->PutScalar(0.0);

  // incremental solution vector with length of all FSI dofs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  iterinc_->PutScalar(0.0);

  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  zeros_->PutScalar(0.0);

  // residual vector with length of all FSI dofs
  rhs_ = LINALG::CreateVector(*DofRowMap(), true);
  rhs_->PutScalar(0.0);
  nall_ = (*rhs_).GlobalLength();

  firstcall_ = true;

  // equilibrium iteration loop (loop over k)
  while ( ((not Converged()) and (iter_ <= itermax_)) or (iter_ ==  1) )
  {
    // compute residual forces #rhs_ and tangent #tang_
    // build linear system stiffness matrix and rhs/force
    // residual for each field

    Evaluate(iterinc_);

    // create the linear system
    // J(x_i) \Delta x_i = - R(x_i)
    // create the systemmatrix
    SetupSystemMatrix();

    // check whether we have a sanely filled tangent matrix
    if (not systemmatrix_->Filled())
    {
      dserror("Effective tangent matrix must be filled here");
    }

    SetupRHS(*rhs_,firstcall_);

    LinearSolve();

    // reset solver tolerance
    solver_->ResetTolerance();

    // build residual force norm
    // for now use for simplicity only L2/Euclidian norm
    rhs_->Norm2(&normrhs_);
    StructureField().RHS()->Norm2(&normstrrhs_);
    FluidField().RHS()->Norm2(&normflrhs_);
    AleField().RHS()->Norm2(&normalerhs_);

    // build residual increment norm
    iterinc_->Norm2(&norminc_);

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;

    firstcall_ = false;

  }// end while loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ( (Converged()) and (Comm().MyPID()==0) )
  {
    cout << endl;
    cout << endl;
    cout << BLUE_LIGHT << "  Newton Converged! " <<  END_COLOR<<  endl;
  }
  else if (iter_ >= itermax_)
  {
    cout << endl;
    cout << endl;
    cout << RED_LIGHT << " Newton unconverged in "<< iter_ << " iterations " << END_COLOR<<  endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicNoNOX::Converged()
{
  // check for single norms
  bool convinc = false;
  bool convfres = false;

  // residual increments
  switch (normtypeinc_)
  {
    case INPAR::FSI::convnorm_abs:
      convinc = norminc_ < tolinc_;
      break;
      //case INPAR::FSI::convnorm_rel:

      //break;
    default:
      dserror("Cannot check for convergence of residual values!");
  }

  // structural, fluid and ale residual forces
  switch (normtypefres_)
  {
  case INPAR::FSI::convnorm_abs:
    convfres = normrhs_ < tolfres_;
    //convfres = (normrhs_/nall_) < tolfres_;
    break;
  case INPAR::FSI::convnorm_rel:
    convfres = ( ((normstrrhs_/ns_) < tolfres_) and ((normflrhs_/nf_) < tolfres_) and ((normalerhs_/na_) < tolfres_));
    //convfres = ( (normstrrhs_) and (normflrhs_)  and (normalerhs_));
    break;
  case INPAR::FSI::convnorm_mix:
    convfres = ( (normstrrhs_ < tolfres_) and (normflrhs_ < tolfres_) and (normalerhs_ < tolfres_) );
    break;
  default:
    dserror("Cannot check for convergence of residual forces!");
  }

  // combined
  bool conv = false;
  if (combincfres_==INPAR::FSI::bop_and)
     conv = convinc and convfres;
   else
     dserror("Something went terribly wrong with binary operator!");

  return conv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::LinearSolve()
{
  // merge blockmatrix to SparseMatrix and solve
  Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

  // apply Dirichlet BCs to system of equations
  if(firstcall_)
    InitialGuess(iterinc_);
  else
    iterinc_->PutScalar(0.0);

  LINALG::ApplyDirichlettoSystem(
    sparse,
    iterinc_,
    rhs_,
    Teuchos::null,
    zeros_,
    *CombinedDBCMap()
    );

  // get UMFPACK...

  // TAW: ???????????? This is nonsense ??????????????????
  /*Teuchos::ParameterList solverparams = DRT::Problem::Instance()->FluidSolverParams();

  solver_ = rcp(new LINALG::Solver(solverparams, Comm(),
                                   DRT::Problem::Instance()->ErrorFile()->Handle()));

  // standard solver call
  solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, true, iter_==1);*/
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
   TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::Evaluate");

   Teuchos::RCP<const Epetra_Vector> sx;
   Teuchos::RCP<const Epetra_Vector> fx;
   Teuchos::RCP<const Epetra_Vector> ax;

   if (x!=Teuchos::null)
   {
     // structure, ale and fluid fields expects the step increment. So
     // we add all of the increments together to build the step
     // increment.
     //
     // The update of the latest increment with iteration increments:
     // x^n+1_i+1 = x^n+1_i + iterinc
     //
     // The update of the latest increment with step increment:
     // x^n+1_i+1 = x^n     + stepinc

     x_sum_->Update(1.0,*x,1.0);

     ExtractFieldVectors(x_sum_,sx,fx,ax);
     if (sdbg_!=Teuchos::null)
     {
       sdbg_->NewIteration();
       sdbg_->WriteVector("x",*StructureField().Interface().ExtractFSICondVector(sx));
     }
   }

   // Call all fileds evaluate method and assemble rhs and matrices

   {
     Epetra_Time ts(Comm());
     StructureField().Evaluate(sx);
     //cout0_  << "structure time: " << ts.ElapsedTime() << endl;
   }

  {
     // ALE field expects the sum of all increments and not the
     // latest increment. It adds the sum of all increments to the
     // displacement of the last time step. So we need to build the
     // sum of all increments and give it to ALE.

     Epetra_Time ta(Comm());
     AleField().Evaluate(ax);
  }

   // transfer the current ale mesh positions to the fluid field
   Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluid(AleField().ExtractDisplacement());
   FluidField().ApplyMeshDisplacement(fluiddisp);

   {
     Epetra_Time tf(Comm());
     FluidField().Evaluate(fx);
     //cout0_ << "fluid time : " << tf.ElapsedTime() << endl;
   }

   //cout0_ << endl;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps)
{
  fullmap_ = LINALG::MultiMapExtractor::MergeMaps(maps);
  blockrowdofmap_.Setup(*fullmap_,maps);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::SetDefaultParameters(const Teuchos::ParameterList& fsidyn,
                                           Teuchos::ParameterList& list)
{
  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = list;

  nlParams.set<std::string>("Nonlinear Solver", "Line Search Based");
  //nlParams.set("Preconditioner", "None");
  //nlParams.set("Norm abs F", fsidyn.get<double>("CONVTOL"));

  nlParams.set("Max Iterations", fsidyn.get<int>("ITEMAX"));
  //nlParams.set("Max Iterations", 1);

  nlParams.set("Norm abs pres", fsidyn.get<double>("CONVTOL"));
  nlParams.set("Norm abs vel",  fsidyn.get<double>("CONVTOL"));
  nlParams.set("Norm abs disp", fsidyn.get<double>("CONVTOL"));

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  //Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");



#if 0 // NonlinearCG

  dirParams.set("Method","NonlinearCG");
  lineSearchParams.set("Method","NonlinearCG");
  Teuchos::ParameterList& nlcgParams = dirParams.sublist("Nonlinear CG");
  ParameterList& lsParams = nlcgParams.sublist("Linear Solver");
  nlcgParams.set("Restart Frequency",10);
  nlcgParams.set("Precondition", "On");
  nlcgParams.set("Orthogonalize", "Fletcher-Reeves");
  lsParams.set("Preconditioning", "User Supplied Preconditioner");
  lsParams.set("Preconditioner","User Defined");

#endif // NonlinearCG

#if 0 // Steepest Decent

  dirParams.set("Method","Steepest Descent");
  //lineSearchParams.set("Method","Full Step");
  Teuchos::ParameterList& SDParams = dirParams.sublist("Steepest Descent");
  ParameterList& lsParams = SDParams.sublist("Linear Solver");
  dirParams.set("Precondition","On");
  SDParams.set("Precondition","On"); // not recognized
  // Quadratic Model Min / 2-Norm / F 2-Norm / None
  SDParams.set("Scaling Type","2-Norm");
  lsParams.set("Preconditioning","User Supplied Preconditioner");
  lsParams.set("Preconditioner","User Defined");

#endif // Steepest Decent


#if 1 // Newton

  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  dirParams.set<std::string>("Method","User Defined");
//   Teuchos::RCP<NOX::Direction::UserDefinedFactory> newtonfactory = Teuchos::rcp(this,false);
//   dirParams.set("User Defined Direction Factory",newtonfactory);


#if 0 // polynominal line search
  lineSearchParams.set("Method","Polynomial");
  Teuchos::ParameterList& polyparams = lineSearchParams.sublist("Polynomial");
  polyparams.set("Default Step",1.0);
  polyparams.set("Max Iters",5);
  polyparams.set("Minimum Step",0.01);
  polyparams.set("Maximum Step",1.5);
  polyparams.set("Recovery Step Type","Constant");
  polyparams.set("Recovery Step",0.1);
  //polyparams.set("Interpolation Type","Cubic");
  polyparams.set("Interpolation Type","Quadratic");
  //polyparams.set("Interpolation Type","Quadratic3");
  //polyparams.set("Sufficient Decrease Condition","Armijo-Goldstein");
  polyparams.set("Sufficient Decrease Condition","Ared/Pred");
  //polyparams.set("Sufficient Decrease Condition","None");
  //polyparams.set("Sufficient Decrease",0.001);
  polyparams.set("Alpha Factor",0.1);
  polyparams.set("Use Counters",false);
  polyparams.set("Maximum Iteration for Increase",0);
  polyparams.set("Optimize Slope Calculation",false);
  polyparams.set("Allowed Relative Increase",1.0);
  polyparams.set("Min Bounds Factor",0.1);
  polyparams.set("Max Bounds Factor",0.95);
#endif

#if 0 // backtracking
  lineSearchParams.set("Method","Backtrack");
  Teuchos::ParameterList& backtrack = lineSearchParams.sublist("Backtrack");
  backtrack.set("Default Step",1.0);
  backtrack.set("Minimum Step",0.01);
  backtrack.set("Recovery Step",0.1);
  backtrack.set("Reduction Factor",0.5);
  backtrack.set("Max Iters",10);
  backtrack.set("Maximum Iteration for Increase",0);
  backtrack.set("Optimize Slope Calculation",false);
  polyparams.set("Min Bounds Factor",0.1);
  polyparams.set("Max Bounds Factor",1.5);
#endif

#if 0 // More'-Thuente
  lineSearchParams.set("Method","More'-Thuente");
  Teuchos::ParameterList& mt = lineSearchParams.sublist("More'-Thuente");
  mt.set("Default Step",1.0);
  mt.set("Minimum Step",0.01);
  mt.set("Maximum Step",1.5);
  mt.set("Max Iters",5);
  mt.set("Recovery Step",0.1);
  //mt.set("Sufficient Decrease Condition","Armijo-Goldstein");
  mt.set("Sufficient Decrease Condition","Ared/Pred");
  mt.set("Alpha Factor",0.1);
  mt.set("Optimize Slope Calculation",false);
  mt.set("Maximum Iteration for Increase",0);
  mt.set("Curvature Condition",0.9999);
  mt.set("Min Bounds Factor",0.1);
  mt.set("Max Bounds Factor",1.1);
#endif


#endif // Newton

  // status tests are expensive, but instructive
  //solverOptions.set<std::string>("Status Test Check Type","Minimal");
  solverOptions.set<std::string>("Status Test Check Type","Complete");

  // be explicit about linear solver parameters
  lsParams.set<std::string>("Aztec Solver","GMRES");
  //lsParams.set<std::string>("BiCGStab","GMRES");
  lsParams.set<std::string>("Orthogonalization","Modified");

  // "r0", "rhs", "norm", "no scaling", "sol"
  lsParams.set<std::string>("Convergence Test","r0");

  lsParams.set<int>("Size of Krylov Subspace",50);
  lsParams.set<int>("Max Iterations",1000);
  lsParams.set<std::string>("Preconditioner","User Defined");
  lsParams.set<int>("Output Frequency",10);
  lsParams.set<bool>("Output Solver Details",true);

  // adaptive tolerance settings
  lsParams.set<double>("base tolerance",fsidyn.get<double>("BASETOL"));
  lsParams.set<double>("adaptive distance",fsidyn.get<double>("ADAPTIVEDIST"));

#if 0
  // add Aitken relaxation to Newton step
  // there is nothing to be gained...
  Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> aitkenfactory =
    Teuchos::rcp(new NOX::FSI::AitkenFactory());
  lineSearchParams.set("Method","User Defined");
  lineSearchParams.set("User Defined Line Search Factory", aitkenfactory);

  //lineSearchParams.sublist("Aitken").set("max step size",
  //fsidyn.get<double>("MAXOMEGA"));
#endif
}

/*----------------------------------------------------------------------*/
/*  print Newton-Raphson iteration to screen and error file             */
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::PrintNewtonIter()
{
  // print to standard out
  // replace myrank_ here general by Comm().MyPID()
  if ( Comm().MyPID()==0 )
  {
    if (iter_== 1)
      PrintNewtonIterHeader(stdout);
    PrintNewtonIterText(stdout);
  }

//   // print to error file
//   if ( printerrfile_ and printiter_ )
//   {
//     if (iter_== 1)
//       PrintNewtonIterHeader(errfile_);
//     PrintNewtonIterText(errfile_);
//   }
}

/*----------------------------------------------------------------------*/
/* print Newton-Raphson iteration to screen and error file              */
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::PrintNewtonIterHeader(FILE* ofile)
{
  cout << "CONVTOL: " << tolfres_ << endl;
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6)<< "numiter  |";

  // different style due relative or absolute error checking
  // displacement
  switch ( normtypefres_ )
  {
  case INPAR::FSI::convnorm_abs :
    oss <<std::setw(18)<< "abs-res-norm  |";
    break;
  case INPAR::FSI::convnorm_rel :
    oss <<std::setw(18)<< "rel-res-norm  |";
    break;
  case INPAR::FSI::convnorm_mix :
    oss <<std::setw(18)<< "mix-res-norm";
    break;
  default:
    dserror("You should not turn up here.");
  }

  switch ( normtypeinc_ )
  {
  case INPAR::FSI::convnorm_abs :
    oss <<std::setw(18)<< "abs-inc-norm |";
    break;
  default:
    dserror("You should not turn up here.");
  }

  // add solution time
  oss << std::setw(12)<< "wct    |";
  cout << "=========================================================="<< endl;

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile==NULL)
    dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);
  cout << "==========================================================";
}

/*---------------------------------------------------------------------*/
/*  print Newton-Raphson iteration to screen                           */
/*---------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::PrintNewtonIterText(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(7)<< iter_ << "/" << itermax_;

  // different style due relative or absolute error checking
  // displacement
  switch ( normtypefres_ )
  {
  case INPAR::FSI::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << (normrhs_/nall_);
    break;
  case INPAR::FSI::convnorm_rel :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << "str " << (normstrrhs_/ns_) << "fl " << (normflrhs_/nf_) << "ale "<< (normalerhs_/na_) ;
    break;
  case INPAR::FSI::convnorm_mix :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << "str " << (normstrrhs_/ns_) << "fl " << (normflrhs_/nf_) << "ale "<< (normalerhs_/na_) ;
    break;
  default:
    dserror("You should not turn up here.");
 }

  switch ( normtypeinc_ )
  {
  case INPAR::FSI::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << norminc_;
    printf("\n");
    break;
  default:
    dserror("You should not turn up here.");
  }

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile==NULL)
    dserror("no ofile available");
  fprintf(ofile, "%s", oss.str().c_str());

  // print it, now
  fflush(ofile);
}

#endif
