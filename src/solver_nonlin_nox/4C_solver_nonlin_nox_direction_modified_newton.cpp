/*----------------------------------------------------------------------------*/
/*! \file

\brief Implementation of a modified Newton approach



\level 3

*/
/*----------------------------------------------------------------------------*/

#include "4C_solver_nonlin_nox_direction_modified_newton.hpp"

#include "4C_solver_nonlin_nox_direction_defaultsteptest.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <NOX_GlobalData.H>
#include <NOX_Solver_Generic.H>
#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Direction::ModifiedNewton::ModifiedNewton(
    const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& p)
    : NOX::Nln::Direction::Newton(gd, p),
      init_primal_diag_corr_(p.sublist("Newton", true)
                                 .sublist("Modified", true)
                                 .get<double>("Initial Primal Diagonal Correction")),
      min_primal_diag_corr_(p.sublist("Newton", true)
                                .sublist("Modified", true)
                                .get<double>("Minimal Primal Diagonal Correction")),
      max_primal_diag_corr_(p.sublist("Newton", true)
                                .sublist("Modified", true)
                                .get<double>("Maximal Primal Diagonal Correction")),
      primal_red_fac_(p.sublist("Newton", true)
                          .sublist("Modified", true)
                          .get<double>("Primal Reduction Factor")),
      primal_acc_fac_(p.sublist("Newton", true)
                          .sublist("Modified", true)
                          .get<double>("Primal Accretion Factor")),
      primal_high_acc_fac_(p.sublist("Newton", true)
                               .sublist("Modified", true)
                               .get<double>("Primal High Accretion Factor")),
      utils_(gd->getUtils()),
      params_(&p)
{
  Teuchos::ParameterList& pmod = p.sublist("Newton", true).sublist("Modified", true);
  fill_default_step_tests(pmod);

  if (p.sublist("Newton").get("Forcing Term Method", "Constant") == "Constant")
    use_adjustable_forcing_term_ = false;
  else
    use_adjustable_forcing_term_ = true;

  fp_except_.shall_be_caught_ = pmod.get<bool>("Catch Floating Point Exceptions");

  int correction_counter = 0;
  pmod.set<int>("Number of Corrections", correction_counter);
  corr_counter_ = &pmod.get<int>("Number of Corrections");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Direction::ModifiedNewton::fill_default_step_tests(
    Teuchos::ParameterList& pmodnewton)
{
  enum DefaultStepTest dstest_type =
      Teuchos::getIntegralValue<DefaultStepTest>(pmodnewton, "Default Step Test");

  switch (dstest_type)
  {
    case NOX::Nln::Direction::DefaultStepTest::none:
    {
      // do nothing
      break;
    }
    case NOX::Nln::Direction::DefaultStepTest::volume_change_control:
    {
      Teuchos::RCP<Test::Generic> dstest = Teuchos::rcp(new Test::VolumeChange(utils_));
      dstests_.push_back(dstest);
      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::ModifiedNewton::compute_correction_direction(::NOX::Abstract::Vector& dir,
    ::NOX::Abstract::Group& grp, const ::NOX::Solver::Generic& solver,
    NOX::Nln::CorrectionType corr_type)
{
  switch (corr_type)
  {
    /* The system matrix has been newly evaluated and, thus, the modification
     * must be repeated. However, no new default step tests are performed,
     * instead the correction of the previous step is used. */
    case NOX::Nln::CorrectionType::soc_full:
    {
      Teuchos::RCP<Epetra_Vector> diagonal = get_diagonal(grp);
      if (not use_unmodified_system()) modify_system(grp, diagonal.get(), primal_diag_corr_last_);
      return solve_modified_system(dir, grp, solver);
    }
    /* In this case the system matrix did not change and still contains the
     * previous modification. Therefore, just a solve call is performed. */
    case NOX::Nln::CorrectionType::soc_cheap:
    {
      return solve_modified_system(dir, grp, solver);
    }
    default:
    {
      FOUR_C_THROW(
          "The NOX::Nln::CorrectionType \"%s\" is not yet supported "
          "by the NOX::Nln::Direction::ModifiedNewton object!",
          NOX::Nln::CorrectionType2String(corr_type).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::ModifiedNewton::compute(
    ::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp, const ::NOX::Solver::Generic& solver)
{
  *corr_counter_ = 0;
  stagnation_counter_ = 0;

  // special treatment for correction steps
  NOX::Nln::Group& nln_grp = dynamic_cast<NOX::Nln::Group&>(grp);
  NOX::Nln::CorrectionType corr_type = nln_grp.get_correction_type();
  if (corr_type != NOX::Nln::CorrectionType::vague)
  {
    const bool corr_status = compute_correction_direction(dir, grp, solver, corr_type);
    print(utils_->out(), &corr_type);

    return corr_status;
  }

  const bool newton_status = solve_unmodified_system(dir, grp, solver);
  bool mod_newton_status = false;

  // Modify the system if a default Newton step failed
  if (not newton_status) mod_newton_status = compute_modified_newton(dir, grp, solver);

  Teuchos::RCP<Epetra_Vector> diagonal = Teuchos::null;
  bool dst_success = test_default_step_quality(dir, grp, diagonal, true);

  while (not dst_success)
  {
    mod_newton_status = compute_modified_newton(dir, grp, solver, diagonal.get());
    if (not mod_newton_status) break;

    dst_success = test_default_step_quality(dir, grp, diagonal);
  }

  // update the successive reduction counter
  update_successive_reduction_counter();

  // print info to screen
  if (*corr_counter_ > 0) print(utils_->out());

  /* If the system has been modified and the modification has been successful,
   * the correction factor is stored for the future. In this way, the history
   * information can be used to avoid unnecessary linear solver calls. */
  if (mod_newton_status) store_correction_factor();

  return (*corr_counter_ == 0 ? newton_status : mod_newton_status);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::ModifiedNewton::use_unmodified_system() const
{
  return (successive_red_counter_ > 2);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Direction::ModifiedNewton::update_successive_reduction_counter()
{
  if (*corr_counter_ <= 1)
    ++successive_red_counter_;
  else
    successive_red_counter_ = 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> NOX::Nln::Direction::ModifiedNewton::get_diagonal(
    const ::NOX::Abstract::Group& grp) const
{
  // fill the diagonal
  Teuchos::RCP<Epetra_Vector> diagonal = Teuchos::null;

  return diagonal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::ModifiedNewton::test_default_step_quality(::NOX::Abstract::Vector& dir,
    ::NOX::Abstract::Group& grp, Teuchos::RCP<Epetra_Vector>& diagonal, bool first_test)
{
  bool status = true;
  // perform the test
  for (auto& dstest : dstests_)
  {
    if (first_test)
      status = dstest->init_and_check_test(dir, grp);
    else
      status = dstest->check_test(dir, grp);
    if (not status) break;
  }

  // fill the diagonal
  diagonal = Teuchos::null;
  diagonal = get_diagonal(grp);

  return status;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::ModifiedNewton::compute_modified_newton(::NOX::Abstract::Vector& dir,
    ::NOX::Abstract::Group& grp, const ::NOX::Solver::Generic& solver, Epetra_Vector* diagonal)
{
  // Compute F and Jacobian at current solution.
  NOX::Nln::Group& nln_grp = dynamic_cast<NOX::Nln::Group&>(grp);
  ::NOX::Abstract::Group::ReturnType status = nln_grp.compute_f_and_jacobian();
  if (status != ::NOX::Abstract::Group::Ok)
    throw_error(__LINE__, "compute", "Unable to compute F and/or Jacobian");

  if (not modify_system(grp, diagonal)) return false;

  // try to solve the modified system
  bool success = solve_modified_system(dir, grp, solver);

  // recursive call if the solving attempt fails
  return (success ? true : compute_modified_newton(dir, grp, solver));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::ModifiedNewton::solve_modified_system(
    ::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp, const ::NOX::Solver::Generic& solver)
{
  // sanity check
  if (not grp.isF()) FOUR_C_THROW("At this point the right hand side must be fully evaluated!");
  if (not grp.isJacobian()) FOUR_C_THROW("At this point the jacobian must be fully evaluated!");

  // Reset the linear solver tolerance.
  if (use_adjustable_forcing_term_)
  {
    resetForcingTerm(grp, solver.getPreviousSolutionGroup(), solver.getNumIterations(), solver);
  }

  // reset the IsValid flag of the previous Newton solver attempt
  dynamic_cast<NOX::Nln::Group&>(grp).reset_is_valid_newton();

  Teuchos::ParameterList& plinsolver =
      params_->sublist("Newton", true).sublist("Linear Solver", true);

  fp_except_.precompute();
  ::NOX::Abstract::Group::ReturnType status = grp.computeNewton(plinsolver);
  if (fp_except_.postcompute(utils_->out(::NOX::Utils::Warning)))
    status = ::NOX::Abstract::Group::ReturnType::Failed;

  dir = grp.getNewton();
  if (status == ::NOX::Abstract::Group::Ok) set_stagnation_counter(dir);

  return (status == ::NOX::Abstract::Group::Ok);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::ModifiedNewton::solve_unmodified_system(
    ::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp, const ::NOX::Solver::Generic& solver)
{
  if (not use_unmodified_system()) return false;

  try
  {
    fp_except_.precompute();
    const bool success = Newton::compute(dir, grp, solver);
    if (fp_except_.postcompute(utils_->out(::NOX::Utils::Warning))) return false;

    return success;
  }
  catch (const char* e)
  {
    utils_->out(::NOX::Utils::Warning)
        << ::NOX::Utils::fill(40, '*') << "\n"
        << "WARNING: Error caught during the linear solver call = " << e << "\n"
        << ::NOX::Utils::fill(40, '*') << "\n";

    return false;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::ModifiedNewton::modify_system(
    ::NOX::Abstract::Group& grp, Epetra_Vector* diagonal)
{
  primal_diag_corr_ = get_primal_diag_correction(*corr_counter_ == 0);

  // the maximal correction factor is reached without success
  if (primal_diag_corr_ > max_primal_diag_corr_) return false;

  return modify_system(grp, diagonal, primal_diag_corr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Direction::ModifiedNewton::set_stagnation_counter(::NOX::Abstract::Vector& dir)
{
  double dir_nrm2 = dir.norm(::NOX::Abstract::Vector::TwoNorm);
  static double prev_dir_nrm2 = 0.0;

  if (std::abs(prev_dir_nrm2 - dir_nrm2) < 1.0e-3 * prev_dir_nrm2) ++stagnation_counter_;

  if (stagnation_counter_ > 3)
    FOUR_C_THROW(
        "There has been %d stagnation detected. This indicates that the "
        "modified Newton does no longer help to improve the solvability "
        "of the problem. Thus, the problem might be badly posed.",
        stagnation_counter_);

  prev_dir_nrm2 = dir_nrm2;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::ModifiedNewton::modify_system(
    ::NOX::Abstract::Group& grp, Epetra_Vector* diagonal, const double primal_diag_corr)
{
  NOX::Nln::Group& nln_grp = dynamic_cast<NOX::Nln::Group&>(grp);

  if (*corr_counter_ == 0)
  {
    // store the copy of the original diagonal
    original_diag_ptr_ = nln_grp.get_diagonal_of_jacobian(0);
  }

  Teuchos::RCP<Epetra_Vector> mod_diag_ptr = Teuchos::null;
  if (diagonal)
  {
    if (not diagonal->Map().PointSameAs(original_diag_ptr_->Map()))
      FOUR_C_THROW("The provided diagonal vector has the wrong map layout!");

    mod_diag_ptr = Teuchos::rcpFromRef(*diagonal);
  }
  else
  {
    const Epetra_BlockMap& map = original_diag_ptr_->Map();
    mod_diag_ptr = Teuchos::rcp(new Epetra_Vector(map, false));
    mod_diag_ptr->PutScalar(1.0);
  }

  mod_diag_ptr->Scale(primal_diag_corr);
  mod_diag_ptr->Update(1.0, *original_diag_ptr_, 1.0);

  nln_grp.replace_diagonal_of_jacobian(*mod_diag_ptr, 0);

  ++(*corr_counter_);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::Direction::ModifiedNewton::get_primal_diag_correction(const bool first) const
{
  if (first) return get_first_primal_diag_correction();

  return get_primal_diag_correction();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::Direction::ModifiedNewton::get_first_primal_diag_correction() const
{
  if (primal_diag_corr_last_ == 0.0)
    return init_primal_diag_corr_;
  else
    return std::max(primal_red_fac_ * primal_diag_corr_last_, min_primal_diag_corr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::Direction::ModifiedNewton::get_primal_diag_correction() const
{
  if (primal_diag_corr_last_ == 0.0)
    return primal_diag_corr_ * primal_high_acc_fac_;
  else
    return primal_diag_corr_ * primal_acc_fac_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Direction::ModifiedNewton::store_correction_factor()
{
  primal_diag_corr_last_ = primal_diag_corr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Direction::ModifiedNewton::print(
    std::ostream& os, const NOX::Nln::CorrectionType* corr_type) const
{
  os << ::NOX::Utils::fill(80, '=') << "\nNOX::Nln::Direction::ModifiedNewton\n"
     << "The system matrix has been " << *corr_counter_ << " time(s) corrected"
     << " before a reliable solution \nof the linear system could be obtained.\n"
     << ::NOX::Utils::fill(80, '=');

  if (corr_type)
    os << "\nCorrection Type = " << NOX::Nln::CorrectionType2String(*corr_type) << "\n"
       << ::NOX::Utils::fill(80, '-');

  os << "\nStagnation detected                       = " << stagnation_counter_
     << "\nSuccessive reduction counter              = " << successive_red_counter_ << "\n"
     << ::NOX::Utils::fill(80, '-')
     << "\nLast primal diagonal correction factor    = " << primal_diag_corr_last_
     << "\nCurrent primal diagonal correction factor = " << primal_diag_corr_ << "\n"
     << ::NOX::Utils::fill(80, '-')
     << "\nInitial primal diagonal correction = " << init_primal_diag_corr_
     << "\nMinimal primal diagonal correction = " << min_primal_diag_corr_
     << "\nMaximal primal diagonal correction = " << max_primal_diag_corr_
     << "\nPrimal reduction factor            = " << primal_red_fac_
     << "\nPrimal acceration factor           = " << primal_acc_fac_
     << "\nPrimal high acceration factor      = " << primal_high_acc_fac_ << "\n"
     << ::NOX::Utils::fill(80, '=') << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Direction::ModifiedNewton::throw_error(
    const int line, const std::string& functionName, const std::string& errorMsg) const
{
  if (utils_->isPrintType(::NOX::Utils::Error))
    utils_->err() << line << " -- "
                  << "NOX::Nln::Direction::ModifiedNewton::" << functionName << " - " << errorMsg
                  << std::endl;
  throw "NOX Error";
}

FOUR_C_NAMESPACE_CLOSE
