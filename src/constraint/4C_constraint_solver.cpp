/*----------------------------------------------------------------------*/
/*! \file
\brief Class containing Uzawa algorithm to solve linear system.
\level 2


*----------------------------------------------------------------------*/


#include "4C_constraint_solver.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method_linalg.hpp"

#include <Teuchos_ParameterList.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 11/07|
 *----------------------------------------------------------------------*/
CONSTRAINTS::ConstraintSolver::ConstraintSolver(Teuchos::RCP<Core::FE::Discretization> discr,
    Core::LinAlg::Solver& solver, Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps,
    Teuchos::ParameterList params)
    : actdisc_(discr),
      max_iter_(params.get<int>("UZAWAMAXITER", 50)),
      dirichtoggle_(Teuchos::null),
      dbcmaps_(dbcmaps)
{
  setup(discr, solver, dbcmaps, params);
}

/*----------------------------------------------------------------------*
 |  set-up (public)                                             tk 11/07|
 *----------------------------------------------------------------------*/
void CONSTRAINTS::ConstraintSolver::setup(Teuchos::RCP<Core::FE::Discretization> discr,
    Core::LinAlg::Solver& solver, Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps,
    Teuchos::ParameterList params)
{
  solver_ = Teuchos::rcp(&solver, false);

  algochoice_ = Core::UTILS::IntegralValue<Inpar::STR::ConSolveAlgo>(params, "UZAWAALGO");

  // different setup for #adapttol_
  isadapttol_ = true;
  isadapttol_ = (Core::UTILS::IntegralValue<int>(params, "ADAPTCONV") == 1);

  // simple parameters
  adaptolbetter_ = params.get<double>("ADAPTCONV_BETTER", 0.01);
  iterationparam_ = params.get<double>("UZAWAPARAM", 1);
  minparam_ = iterationparam_ * 1E-3;
  iterationtol_ = params.get<double>("UZAWATOL", 1E-8);


  counter_ = 0;
  return;
}



/*----------------------------------------------------------------------*
|(public)                                                               |
|Solve linear constrained system                                        |
*-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstraintSolver::Solve(Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> constr,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> constrT, Teuchos::RCP<Epetra_Vector> dispinc,
    Teuchos::RCP<Epetra_Vector> lagrinc, const Teuchos::RCP<Epetra_Vector> rhsstand,
    const Teuchos::RCP<Epetra_Vector> rhsconstr)
{
  switch (algochoice_)
  {
    case Inpar::STR::consolve_uzawa:
      solve_uzawa(stiff, constr, constrT, dispinc, lagrinc, rhsstand, rhsconstr);
      break;
    case Inpar::STR::consolve_direct:
      solve_direct(stiff, constr, constrT, dispinc, lagrinc, rhsstand, rhsconstr);
      break;
    case Inpar::STR::consolve_simple:
      solve_simple(stiff, constr, constrT, dispinc, lagrinc, rhsstand, rhsconstr);
      break;
    default:
      FOUR_C_THROW("Unknown constraint solution technique!");
  }
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                               |
|Solve linear constrained system by iterative Uzawa algorithm           |
*-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstraintSolver::solve_uzawa(Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> constr,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> constrT, Teuchos::RCP<Epetra_Vector> dispinc,
    Teuchos::RCP<Epetra_Vector> lagrinc, const Teuchos::RCP<Epetra_Vector> rhsstand,
    const Teuchos::RCP<Epetra_Vector> rhsconstr)
{
  const int myrank = (actdisc_->Comm().MyPID());
  // For every iteration step an uzawa algorithm is used to solve the linear system.
  // Preparation of uzawa method to solve the linear system.
  double norm_uzawa;
  double norm_uzawa_old;
  double quotient;
  double norm_constr_uzawa;
  int numiter_uzawa = 0;
  // counter used for adaptivity
  const int adaptstep = 2;
  const int minstep = 1;
  int count_paramadapt = 1;

  const double computol = 1E-8;

  Teuchos::RCP<Epetra_Vector> constrTLagrInc = Teuchos::rcp(new Epetra_Vector(rhsstand->Map()));
  Teuchos::RCP<Epetra_Vector> constrTDispInc = Teuchos::rcp(new Epetra_Vector(rhsconstr->Map()));
  // Core::LinAlg::SparseMatrix constrT =
  // *(Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(constr));

  // ONLY compatability
  // dirichtoggle_ changed and we need to rebuild associated DBC maps
  if (dirichtoggle_ != Teuchos::null)
    dbcmaps_ = Core::LinAlg::ConvertDirichletToggleVectorToMaps(dirichtoggle_);

  Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(rhsstand->Map(), true));
  Teuchos::RCP<Epetra_Vector> dirichzeros = dbcmaps_->extract_cond_vector(zeros);

  // Compute residual of the uzawa algorithm
  Teuchos::RCP<Epetra_Vector> fresmcopy = Teuchos::rcp(new Epetra_Vector(*rhsstand));
  Epetra_Vector uzawa_res(*fresmcopy);
  (*stiff).Multiply(false, *dispinc, uzawa_res);
  uzawa_res.Update(1.0, *fresmcopy, -1.0);

  // blank residual DOFs which are on Dirichlet BC
  dbcmaps_->insert_cond_vector(dirichzeros, Teuchos::rcp(&uzawa_res, false));

  uzawa_res.Norm2(&norm_uzawa);
  Epetra_Vector constr_res(lagrinc->Map());

  constr_res.Update(1.0, *(rhsconstr), 0.0);
  constr_res.Norm2(&norm_constr_uzawa);
  quotient = 1;
  // Solve one iteration step with augmented lagrange
  // Since we calculate displacement norm as well, at least one step has to be taken
  while (((norm_uzawa > iterationtol_ or norm_constr_uzawa > iterationtol_) and
             numiter_uzawa < max_iter_) or
         numiter_uzawa < minstep)
  {
    // solve for disi
    // Solve K . IncD = -R  ===>  IncD_{n+1}
    Core::LinAlg::SolverParams solver_params;
    if (isadapttol_ && counter_ && numiter_uzawa)
    {
      solver_params.nonlin_tolerance = tolres_ / 10.0;
      solver_params.nonlin_residual = norm_uzawa;
      solver_params.lin_tol_better = adaptolbetter_;
    }

    solver_params.refactor = true;
    solver_params.reset = numiter_uzawa == 0 && counter_ == 0;
    solver_->Solve(stiff->EpetraMatrix(), dispinc, fresmcopy, solver_params);
    solver_->ResetTolerance();

    // compute Lagrange multiplier increment
    constrTDispInc->PutScalar(0.0);
    constrT->Multiply(true, *dispinc, *constrTDispInc);
    lagrinc->Update(iterationparam_, *constrTDispInc, iterationparam_, *rhsconstr, 1.0);

    // Compute residual of the uzawa algorithm
    constr->Multiply(false, *lagrinc, *constrTLagrInc);

    fresmcopy->Update(-1.0, *constrTLagrInc, 1.0, *rhsstand, 0.0);
    Epetra_Vector uzawa_res(*fresmcopy);
    (*stiff).Multiply(false, *dispinc, uzawa_res);
    uzawa_res.Update(1.0, *fresmcopy, -1.0);

    // blank residual DOFs which are on Dirichlet BC
    dbcmaps_->insert_cond_vector(dirichzeros, Teuchos::rcp(&uzawa_res, false));
    norm_uzawa_old = norm_uzawa;
    uzawa_res.Norm2(&norm_uzawa);
    Epetra_Vector constr_res(lagrinc->Map());

    constr_res.Update(1.0, *constrTDispInc, 1.0, *rhsconstr, 0.0);
    constr_res.Norm2(&norm_constr_uzawa);
    //-------------Adapt Uzawa parameter--------------
    // For a constant parameter the quotient of two successive residual norms
    // stays nearly constant during the computation. So this quotient seems to be a good
    // measure for the parameter choice
    // Adaptivity only takes place every second step. Otherwise the quotient is not significant.
    if (count_paramadapt >= adaptstep)
    {
      double quotient_new = norm_uzawa / norm_uzawa_old;
      // In case of divergence the parameter must be too high
      if (quotient_new > (1. + computol))
      {
        if (iterationparam_ > 2. * minparam_) iterationparam_ = iterationparam_ / 2.;
        quotient = 1;
      }
      else
      {
        // In case the newly computed quotient is better than the one obtained from the
        // previous parameter, the parameter is increased by a factor (1+quotient_new)
        if (quotient >= quotient_new)
        {
          iterationparam_ = iterationparam_ * (1. + quotient_new);
          quotient = quotient_new;
        }
        // In case the newly computed quotient is worse than the one obtained from the
        // previous parameter, the parameter is decreased by a factor 1/(1+quotient_new)
        else
        {
          if (iterationparam_ > 2. * minparam_)
            iterationparam_ = iterationparam_ / (1. + quotient_new);
          quotient = quotient_new;
        }
      }

      if (iterationparam_ <= minparam_)
      {
        if (!myrank)
          std::cout << "leaving uzawa loop since Uzawa parameter is too low" << std::endl;
        iterationparam_ *= 1E2;
        break;
      }
      count_paramadapt = 0;
    }
    count_paramadapt++;
    numiter_uzawa++;
  }  // Uzawa loop

  if (!myrank)
  {
    std::cout << "Uzawa steps " << numiter_uzawa << ", Uzawa parameter: " << iterationparam_;
    std::cout << ", residual norms for linear system: " << norm_constr_uzawa << " and "
              << norm_uzawa << std::endl;
  }
  counter_++;
  return;
}

/*----------------------------------------------------------------------*
|(public)                                                               |
|Solve linear constrained system by iterative Uzawa algorithm           |
*-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstraintSolver::solve_direct(Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> constr,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> constrT, Teuchos::RCP<Epetra_Vector> dispinc,
    Teuchos::RCP<Epetra_Vector> lagrinc, const Teuchos::RCP<Epetra_Vector> rhsstand,
    const Teuchos::RCP<Epetra_Vector> rhsconstr)
{
  // define maps of standard dofs and additional lagrange multipliers
  Teuchos::RCP<Epetra_Map> standrowmap = Teuchos::rcp(new Epetra_Map(stiff->RowMap()));
  Teuchos::RCP<Epetra_Map> conrowmap = Teuchos::rcp(new Epetra_Map(constr->DomainMap()));
  // merge maps to one large map
  Teuchos::RCP<Epetra_Map> mergedmap = Core::LinAlg::MergeMap(standrowmap, conrowmap, false);
  // define MapExtractor
  Core::LinAlg::MapExtractor mapext(*mergedmap, standrowmap, conrowmap);

  // initialize large Sparse Matrix and Epetra_Vectors
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mergedmatrix =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mergedmap, 81));
  Teuchos::RCP<Epetra_Vector> mergedrhs = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  Teuchos::RCP<Epetra_Vector> mergedsol = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  // ONLY compatability
  // dirichtoggle_ changed and we need to rebuild associated DBC maps
  if (dirichtoggle_ != Teuchos::null)
    dbcmaps_ = Core::LinAlg::ConvertDirichletToggleVectorToMaps(dirichtoggle_);
  // fill merged matrix using Add
  mergedmatrix->Add(*stiff, false, 1.0, 1.0);
  mergedmatrix->Add(*constr, false, 1.0, 1.0);
  mergedmatrix->Add(*constrT, true, 1.0, 1.0);
  mergedmatrix->Complete(*mergedmap, *mergedmap);
  // fill merged vectors using Export
  Core::LinAlg::Export(*rhsconstr, *mergedrhs);
  mergedrhs->Scale(-1.0);
  Core::LinAlg::Export(*rhsstand, *mergedrhs);

  // solve
  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = counter_ == 0;
  solver_->Solve(mergedmatrix->EpetraMatrix(), mergedsol, mergedrhs, solver_params);
  solver_->ResetTolerance();
  // store results in smaller vectors
  mapext.extract_cond_vector(mergedsol, dispinc);
  mapext.extract_other_vector(mergedsol, lagrinc);

  counter_++;
  return;
}

void CONSTRAINTS::ConstraintSolver::solve_simple(Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> constr,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> constrT, Teuchos::RCP<Epetra_Vector> dispinc,
    Teuchos::RCP<Epetra_Vector> lagrinc, const Teuchos::RCP<Epetra_Vector> rhsstand,
    const Teuchos::RCP<Epetra_Vector> rhsconstr)
{
  // row maps (assumed to equal to range map) and extractor
  Teuchos::RCP<Epetra_Map> standrowmap = Teuchos::rcp(new Epetra_Map(stiff->RowMap()));
  Teuchos::RCP<Epetra_Map> conrowmap = Teuchos::rcp(new Epetra_Map(constr->DomainMap()));
  Teuchos::RCP<Epetra_Map> mergedrowmap = Core::LinAlg::MergeMap(standrowmap, conrowmap, false);
  Core::LinAlg::MapExtractor rowmapext(*mergedrowmap, conrowmap, standrowmap);

  // domain maps and extractor
  Teuchos::RCP<Epetra_Map> standdommap = Teuchos::rcp(new Epetra_Map(stiff->DomainMap()));
  Teuchos::RCP<Epetra_Map> condommap = Teuchos::rcp(new Epetra_Map(constr->DomainMap()));
  Teuchos::RCP<Epetra_Map> mergeddommap = Core::LinAlg::MergeMap(standdommap, condommap, false);
  Core::LinAlg::MapExtractor dommapext(*mergeddommap, condommap, standdommap);

  // cast constraint operators to matrices and save transpose of constraint matrix
  Core::LinAlg::SparseMatrix constrTrans(*conrowmap, 81, false, true);
  constrTrans.Add(*constrT, true, 1.0, 0.0);
  constrTrans.Complete(constrT->RangeMap(), constrT->DomainMap());

  // ONLY compatability
  // dirichtoggle_ changed and we need to rebuild associated DBC maps
  if (dirichtoggle_ != Teuchos::null)
    dbcmaps_ = Core::LinAlg::ConvertDirichletToggleVectorToMaps(dirichtoggle_);

  // stuff needed for Dirichlet BCs
  Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(rhsstand->Map(), true));
  Teuchos::RCP<Epetra_Vector> dirichzeros = dbcmaps_->extract_cond_vector(zeros);
  Teuchos::RCP<Epetra_Vector> rhscopy = Teuchos::rcp(new Epetra_Vector(*rhsstand));

  // FIXME: The solver should not be taken from the contact dynamic section here,
  // but must be specified somewhere else instead (popp 11/2012)

  /*
  //make solver CheapSIMPLE-ready
  // meshtying/contact for structure
  const Teuchos::ParameterList& mcparams = Global::Problem::Instance()->contact_dynamic_params();
  // get the solver number used for meshtying/contact problems
  const int linsolvernumber = mcparams.get<int>("LINEAR_SOLVER");
  // check if the meshtying/contact solver has a valid solver number
  if (linsolvernumber == (-1))
   FOUR_C_THROW("no linear solver defined for meshtying/contact problem. Please set LINEAR_SOLVER in
  CONTACT DYNAMIC to a valid number!");
   */

  Teuchos::ParameterList sfparams =
      solver_->Params();  // save copy of original solver parameter list
  const Teuchos::ParameterList& mcparams = Global::Problem::Instance()->contact_dynamic_params();
  const int linsolvernumber = mcparams.get<int>("LINEAR_SOLVER");
  solver_->Params() = Core::LinAlg::Solver::translate_solver_parameters(
      Global::Problem::Instance()->SolverParams(linsolvernumber),
      Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY"));

  // Teuchos::ParameterList sfparams = solver_->Params();  // save copy of original solver parameter
  // list solver_->Params() =
  // Core::LinAlg::Solver::translate_solver_parameters(Global::Problem::Instance()->SolverParams(linsolvernumber));
  if (!solver_->Params().isSublist("Belos Parameters")) FOUR_C_THROW("Iterative solver expected!");

  solver_->Params().set<bool>(
      "CONSTRAINT", true);  // handling of constraint null space within Simple type preconditioners
  solver_->Params().sublist(
      "CheapSIMPLE Parameters");  // this automatically sets preconditioner to CheapSIMPLE!
  /*solver_->Params().sublist("Inverse1") = sfparams;
  // get the solver number used for meshtying/contact problems
  const int simplersolvernumber = mcparams.get<int>("SIMPLER_SOLVER");
  // check if the SIMPLER solver has a valid solver number
  if (simplersolvernumber == (-1))
    FOUR_C_THROW("no linear solver defined for Lagrange multipliers. Please set SIMPLER_SOLVER in
  CONTACT DYNAMIC to a valid number!"); solver_->put_solver_params_to_sub_params("Inverse2",
      Global::Problem::Instance()->SolverParams(simplersolvernumber));
  */

  // build block matrix for SIMPLE
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>> mat =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          dommapext, rowmapext, 81, false, false));
  mat->Assign(0, 0, Core::LinAlg::View, *stiff);
  mat->Assign(0, 1, Core::LinAlg::View, *constr);
  mat->Assign(1, 0, Core::LinAlg::View, constrTrans);
  mat->Complete();

  // merged rhs using Export
  Teuchos::RCP<Epetra_Vector> mergedrhs = Teuchos::rcp(new Epetra_Vector(*mergedrowmap));
  Core::LinAlg::Export(*rhsconstr, *mergedrhs);
  mergedrhs->Scale(-1.0);
  Core::LinAlg::Export(*rhscopy, *mergedrhs);

  // solution vector
  Teuchos::RCP<Epetra_Vector> mergedsol = Teuchos::rcp(new Epetra_Vector(*mergedrowmap));

  // solve
  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = counter_ == 0;
  solver_->Solve(mat->EpetraOperator(), mergedsol, mergedrhs, solver_params);
  solver_->ResetTolerance();
  solver_->Params() = sfparams;  // store back original parameter list

  // store results in smaller vectors
  rowmapext.extract_cond_vector(mergedsol, lagrinc);
  rowmapext.extract_other_vector(mergedsol, dispinc);

  counter_++;
  return;
}

FOUR_C_NAMESPACE_CLOSE
