/*----------------------------------------------------------------------------*/
/*! \file

\brief Implementation of a modified Newton approach

\maintainer Anh-Tu Vuong


\level 3

*/
/*----------------------------------------------------------------------------*/

#include "nox_nln_direction_modified_newton.H"
#include "nox_nln_group.H"
#include "nox_nln_direction_defaultsteptest.H"

#include "../drt_lib/drt_dserror.H"

#include <NOX_GlobalData.H>
#include <NOX_Solver_Generic.H>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include <Teuchos_StandardParameterEntryValidators.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Direction::ModifiedNewton::ModifiedNewton(
    const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& p)
    : NOX::NLN::Direction::Newton(gd, p),
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
  fillDefaultStepTests(pmod);

  if (p.sublist("Newton").get("Forcing Term Method", "Constant") == "Constant")
    useAdjustableForcingTerm_ = false;
  else
    useAdjustableForcingTerm_ = true;

  fp_except_.shall_be_caught_ = pmod.get<bool>("Catch Floating Point Exceptions");

  int correction_counter = 0;
  pmod.set<int>("Number of Corrections", correction_counter);
  corr_counter_ = &pmod.get<int>("Number of Corrections");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Direction::ModifiedNewton::fillDefaultStepTests(Teuchos::ParameterList& pmodnewton)
{
  enum DefaultStepTest dstest_type =
      Teuchos::getIntegralValue<DefaultStepTest>(pmodnewton, "Default Step Test");

  switch (dstest_type)
  {
    case NOX::NLN::Direction::DefaultStepTest::none:
    {
      // do nothing
      break;
    }
    case NOX::NLN::Direction::DefaultStepTest::volume_change_control:
    {
      Teuchos::RCP<Test::Generic> dstest = Teuchos::rcp(new Test::VolumeChange(utils_));
      dstests_.push_back(dstest);
      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::ModifiedNewton::computeCorrectionDirection(NOX::Abstract::Vector& dir,
    NOX::Abstract::Group& grp, const NOX::Solver::Generic& solver,
    NOX::NLN::CorrectionType corr_type)
{
  switch (corr_type)
  {
    /* The system matrix has been newly evaluated and, thus, the modification
     * must be repeated. However, no new default step tests are performed,
     * instead the correction of the previous step is used. */
    case NOX::NLN::CorrectionType::soc_full:
    {
      Teuchos::RCP<Epetra_Vector> diagonal = getDiagonal(grp);
      if (not useUnmodifiedSystem()) modifySystem(grp, diagonal.get(), primal_diag_corr_last_);
      return solveModifiedSystem(dir, grp, solver);
    }
    /* In this case the system matrix did not change and still contains the
     * previous modification. Therefore, just a solve call is performed. */
    case NOX::NLN::CorrectionType::soc_cheap:
    {
      return solveModifiedSystem(dir, grp, solver);
    }
    default:
    {
      dserror(
          "The NOX::NLN::CorrectionType \"%s\" is not yet supported "
          "by the NOX::NLN::Direction::ModifiedNewton object!",
          NOX::NLN::CorrectionType2String(corr_type).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::ModifiedNewton::compute(
    NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, const NOX::Solver::Generic& solver)
{
  *corr_counter_ = 0;
  stagnation_counter_ = 0;

  // special treatment for correction steps
  NOX::NLN::Group& nln_grp = dynamic_cast<NOX::NLN::Group&>(grp);
  NOX::NLN::CorrectionType corr_type = nln_grp.GetCorrectionType();
  if (corr_type != NOX::NLN::CorrectionType::vague)
  {
    const bool corr_status = computeCorrectionDirection(dir, grp, solver, corr_type);
    print(utils_->out(), &corr_type);

    return corr_status;
  }

  const bool newton_status = solveUnmodifiedSystem(dir, grp, solver);
  bool mod_newton_status = false;

  // Modify the system if a default Newton step failed
  if (not newton_status) mod_newton_status = computeModifiedNewton(dir, grp, solver);

  Teuchos::RCP<Epetra_Vector> diagonal = Teuchos::null;
  bool dst_success = testDefaultStepQuality(dir, grp, diagonal, true);

  while (not dst_success)
  {
    mod_newton_status = computeModifiedNewton(dir, grp, solver, diagonal.get());
    if (not mod_newton_status) break;

    dst_success = testDefaultStepQuality(dir, grp, diagonal);
  }

  // update the successive reduction counter
  updateSuccessiveReductionCounter();

  // print info to screen
  if (*corr_counter_ > 0) print(utils_->out());

  /* If the system has been modified and the modification has been successful,
   * the correction factor is stored for the future. In this way, the history
   * information can be used to avoid unnecessary linear solver calls. */
  if (mod_newton_status) storeCorrectionFactor();

  return (*corr_counter_ == 0 ? newton_status : mod_newton_status);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::ModifiedNewton::useUnmodifiedSystem() const
{
  return (successive_red_counter_ > 2);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Direction::ModifiedNewton::updateSuccessiveReductionCounter()
{
  if (*corr_counter_ <= 1)
    ++successive_red_counter_;
  else
    successive_red_counter_ = 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> NOX::NLN::Direction::ModifiedNewton::getDiagonal(
    const NOX::Abstract::Group& grp) const
{
  // fill the diagonal
  Teuchos::RCP<Epetra_Vector> diagonal = Teuchos::null;
  //  for ( auto& dstest : dstests_ )
  //  {
  //    if ( diagonal.is_null() )
  //      diagonal = dstest->getCurrentDiagonal( grp );
  //    else
  //      dstest->fillDiagonal( *diagonal );
  //  }
  //
  //  // debugging
  //  if ( not diagonal.is_null() )
  //  {
  //    double num_non_zeros = 0;
  //    diagonal->Norm2( &num_non_zeros );
  //    utils_->out() << "The diagonal of the jacobian is going to be modified at "
  //        << num_non_zeros*num_non_zeros << " entries!\n";
  //  }

  return diagonal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::ModifiedNewton::testDefaultStepQuality(NOX::Abstract::Vector& dir,
    NOX::Abstract::Group& grp, Teuchos::RCP<Epetra_Vector>& diagonal, bool first_test)
{
  bool status = true;
  // perform the test
  for (auto& dstest : dstests_)
  {
    if (first_test)
      status = dstest->initAndCheckTest(dir, grp);
    else
      status = dstest->checkTest(dir, grp);
    if (not status) break;
  }

  // fill the diagonal
  diagonal = Teuchos::null;
  diagonal = getDiagonal(grp);

  return status;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::ModifiedNewton::computeModifiedNewton(NOX::Abstract::Vector& dir,
    NOX::Abstract::Group& grp, const NOX::Solver::Generic& solver, Epetra_Vector* diagonal)
{
  // Compute F and Jacobian at current solution.
  NOX::NLN::Group& nln_grp = dynamic_cast<NOX::NLN::Group&>(grp);
  NOX::Abstract::Group::ReturnType status = nln_grp.computeFandJacobian();
  if (status != NOX::Abstract::Group::Ok)
    throwError(__LINE__, "compute", "Unable to compute F and/or Jacobian");

  if (not modifySystem(grp, diagonal)) return false;

  // try to solve the modified system
  bool success = solveModifiedSystem(dir, grp, solver);

  // recursive call if the solving attempt fails
  return (success ? true : computeModifiedNewton(dir, grp, solver));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::ModifiedNewton::solveModifiedSystem(
    NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, const NOX::Solver::Generic& solver)
{
  // sanity check
  if (not grp.isF()) dserror("At this point the right hand side must be fully evaluated!");
  if (not grp.isJacobian()) dserror("At this point the jacobian must be fully evaluated!");

  // Reset the linear solver tolerance.
  if (useAdjustableForcingTerm_)
  {
    resetForcingTerm(grp, solver.getPreviousSolutionGroup(), solver.getNumIterations(), solver);
  }

  // reset the IsValid flag of the previous Newton solver attempt
  dynamic_cast<NOX::NLN::Group&>(grp).resetIsValidNewton();

  Teuchos::ParameterList& plinsolver =
      params_->sublist("Newton", true).sublist("Linear Solver", true);

  fp_except_.precompute();
  NOX::Abstract::Group::ReturnType status = grp.computeNewton(plinsolver);
  if (fp_except_.postcompute(utils_->out(NOX::Utils::Warning)))
    status = NOX::Abstract::Group::ReturnType::Failed;

  dir = grp.getNewton();
  if (status == NOX::Abstract::Group::Ok) setStagnationCounter(dir);

  return (status == NOX::Abstract::Group::Ok);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::ModifiedNewton::solveUnmodifiedSystem(
    NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, const NOX::Solver::Generic& solver)
{
  if (not useUnmodifiedSystem()) return false;

  try
  {
    fp_except_.precompute();
    const bool success = Newton::compute(dir, grp, solver);
    if (fp_except_.postcompute(utils_->out(NOX::Utils::Warning))) return false;

    return success;
  }
  catch (const char* e)
  {
    utils_->out(NOX::Utils::Warning)
        << NOX::Utils::fill(40, '*') << "\n"
        << "WARNING: Error caught during the linear solver call = " << e << "\n"
        << NOX::Utils::fill(40, '*') << "\n";

    return false;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::ModifiedNewton::modifySystem(
    NOX::Abstract::Group& grp, Epetra_Vector* diagonal)
{
  primal_diag_corr_ = getPrimalDiagCorrection(*corr_counter_ == 0);

  // the maximal correction factor is reached without success
  if (primal_diag_corr_ > max_primal_diag_corr_) return false;

  return modifySystem(grp, diagonal, primal_diag_corr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Direction::ModifiedNewton::setStagnationCounter(NOX::Abstract::Vector& dir)
{
  double dir_nrm2 = dir.norm(NOX::Abstract::Vector::TwoNorm);
  static double prev_dir_nrm2 = 0.0;

  if (std::abs(prev_dir_nrm2 - dir_nrm2) < 1.0e-3 * prev_dir_nrm2) ++stagnation_counter_;

  if (stagnation_counter_ > 3)
    dserror(
        "There has been %d stagnation detected. This indicates that the "
        "modified Newton does no longer help to improve the solvability "
        "of the problem. Thus, the problem might be badly posed.",
        stagnation_counter_);

  prev_dir_nrm2 = dir_nrm2;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::ModifiedNewton::modifySystem(
    NOX::Abstract::Group& grp, Epetra_Vector* diagonal, const double primal_diag_corr)
{
  NOX::NLN::Group& nln_grp = dynamic_cast<NOX::NLN::Group&>(grp);

  if (*corr_counter_ == 0)
  {
    // store the copy of the original diagonal
    original_diag_ptr_ = nln_grp.getDiagonalOfJacobian(0);
  }

  Teuchos::RCP<Epetra_Vector> mod_diag_ptr = Teuchos::null;
  if (diagonal)
  {
    if (not diagonal->Map().PointSameAs(original_diag_ptr_->Map()))
      dserror("The provided diagonal vector has the wrong map layout!");

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

  nln_grp.replaceDiagonalOfJacobian(*mod_diag_ptr, 0);

  ++(*corr_counter_);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::Direction::ModifiedNewton::getPrimalDiagCorrection(const bool first) const
{
  if (first) return getFirstPrimalDiagCorrection();

  return getPrimalDiagCorrection();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::Direction::ModifiedNewton::getFirstPrimalDiagCorrection() const
{
  if (primal_diag_corr_last_ == 0.0)
    return init_primal_diag_corr_;
  else
    return std::max(primal_red_fac_ * primal_diag_corr_last_, min_primal_diag_corr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::Direction::ModifiedNewton::getPrimalDiagCorrection() const
{
  if (primal_diag_corr_last_ == 0.0)
    return primal_diag_corr_ * primal_high_acc_fac_;
  else
    return primal_diag_corr_ * primal_acc_fac_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Direction::ModifiedNewton::storeCorrectionFactor()
{
  primal_diag_corr_last_ = primal_diag_corr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Direction::ModifiedNewton::print(
    std::ostream& os, const NOX::NLN::CorrectionType* corr_type) const
{
  os << NOX::Utils::fill(80, '=') << "\nNOX::NLN::Direction::ModifiedNewton\n"
     << "The system matrix has been " << *corr_counter_ << " time(s) corrected"
     << " before a reliable solution \nof the linear system could be obtained.\n"
     << NOX::Utils::fill(80, '=');

  if (corr_type)
    os << "\nCorrection Type = " << NOX::NLN::CorrectionType2String(*corr_type) << "\n"
       << NOX::Utils::fill(80, '-');

  os << "\nStagnation detected                       = " << stagnation_counter_
     << "\nSuccessive reduction counter              = " << successive_red_counter_ << "\n"
     << NOX::Utils::fill(80, '-')
     << "\nLast primal diagonal correction factor    = " << primal_diag_corr_last_
     << "\nCurrent primal diagonal correction factor = " << primal_diag_corr_ << "\n"
     << NOX::Utils::fill(80, '-')
     << "\nInitial primal diagonal correction = " << init_primal_diag_corr_
     << "\nMinimal primal diagonal correction = " << min_primal_diag_corr_
     << "\nMaximal primal diagonal correction = " << max_primal_diag_corr_
     << "\nPrimal reduction factor            = " << primal_red_fac_
     << "\nPrimal acceration factor           = " << primal_acc_fac_
     << "\nPrimal high acceration factor      = " << primal_high_acc_fac_ << "\n"
     << NOX::Utils::fill(80, '=') << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Direction::ModifiedNewton::throwError(
    const int line, const std::string& functionName, const std::string& errorMsg) const
{
  if (utils_->isPrintType(NOX::Utils::Error))
    utils_->err() << line << " -- "
                  << "NOX::NLN::Direction::ModifiedNewton::" << functionName << " - " << errorMsg
                  << std::endl;
  throw "NOX Error";
}
