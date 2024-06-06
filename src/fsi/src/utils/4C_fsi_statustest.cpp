/*----------------------------------------------------------------------*/
/*! \file

\brief Test routines for monolithic FSI convergence test

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_fsi_statustest.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fsi_nox_newton.hpp"
#include "4C_utils_exceptions.hpp"

#include <NOX_Abstract_Group.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Solver_Generic.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::FSI::GenericNormF::GenericNormF(
    std::string name, double tolerance, ::NOX::Abstract::Vector::NormType normType, ScaleType stype)
    : status_(::NOX::StatusTest::Unevaluated),
      norm_type_(normType),
      scale_type_(stype),
      specified_tolerance_(tolerance),
      initial_tolerance_(1.0),
      true_tolerance_(tolerance),
      norm_f_(0.0),
      name_(name)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double NOX::FSI::GenericNormF::compute_norm(const Epetra_Vector& v)
{
  int n = v.GlobalLength();
  double norm;
  int err;

  switch (norm_type_)
  {
    case ::NOX::Abstract::Vector::TwoNorm:
      err = v.Norm2(&norm);
      if (err != 0) FOUR_C_THROW("norm failed");
      if (scale_type_ == Scaled) norm /= sqrt(1.0 * n);
      break;

    case ::NOX::Abstract::Vector::OneNorm:
      err = v.Norm1(&norm);
      if (err != 0) FOUR_C_THROW("norm failed");
      if (scale_type_ == Scaled) norm /= n;
      break;

    case ::NOX::Abstract::Vector::MaxNorm:
      err = v.NormInf(&norm);
      if (err != 0) FOUR_C_THROW("norm failed");
      if (scale_type_ == Scaled)
        FOUR_C_THROW("It does not make sense to scale a MaxNorm by the vector length.");
      break;

    default:
      FOUR_C_THROW("norm type confusion");
      break;
  }

  return norm;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
::NOX::StatusTest::StatusType NOX::FSI::GenericNormF::checkStatus(
    const ::NOX::Solver::Generic& problem, ::NOX::StatusTest::CheckType checkType)
{
  if (checkType == ::NOX::StatusTest::None)
  {
    norm_f_ = 0.0;
    status_ = ::NOX::StatusTest::Unevaluated;
  }
  else
  {
    norm_f_ = compute_norm(problem.getSolutionGroup());
    if ((norm_f_ != -1) and (norm_f_ < true_tolerance_))
    {
      status_ = ::NOX::StatusTest::Converged;
    }
    else
    {
      status_ = ::NOX::StatusTest::Unconverged;
    }
  }

  return status_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
::NOX::StatusTest::StatusType NOX::FSI::GenericNormF::getStatus() const { return status_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& NOX::FSI::GenericNormF::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j++) stream << ' ';

  stream << status_  // test status
         << name_    // what is tested?
         << " ";

  // check which norm is used and print its name
  if (norm_type_ == ::NOX::Abstract::Vector::TwoNorm)
    stream << "Two-Norm";
  else if (norm_type_ == ::NOX::Abstract::Vector::OneNorm)
    stream << "One-Norm";
  else if (norm_type_ == ::NOX::Abstract::Vector::MaxNorm)
    stream << "Max-Norm";

  // print current value of norm and given tolerance
  stream << " = " << ::NOX::Utils::sciformat(norm_f_, 3) << " < "
         << ::NOX::Utils::sciformat(true_tolerance_, 3) << "\n";

  // Note: All norms are hard-coded absolute norms. So, we do not need to print
  // this.                                                      mayt.mt 01/2012

  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double NOX::FSI::GenericNormF::getNormF() const { return norm_f_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double NOX::FSI::GenericNormF::getTrueTolerance() const { return true_tolerance_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double NOX::FSI::GenericNormF::get_specified_tolerance() const { return specified_tolerance_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double NOX::FSI::GenericNormF::getInitialTolerance() const { return initial_tolerance_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::FSI::PartialNormF::PartialNormF(std::string name,
    const Core::LinAlg::MultiMapExtractor& extractor, int blocknum, double tolerance,
    ::NOX::Abstract::Vector::NormType normType, ScaleType stype)
    : AdaptiveNewtonNormF(name, tolerance, normType, stype),
      extractor_(extractor),
      blocknum_(blocknum)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double NOX::FSI::PartialNormF::compute_norm(const ::NOX::Abstract::Group& grp)
{
  if (!grp.isF()) return -1.0;

  // extract the block epetra vector

  const ::NOX::Abstract::Vector& abstract_f = grp.getF();
  const ::NOX::Epetra::Vector& f = Teuchos::dyn_cast<const ::NOX::Epetra::Vector>(abstract_f);

  // extract the inner vector elements we are interested in

  Teuchos::RCP<Epetra_Vector> v = extractor_.ExtractVector(f.getEpetraVector(), blocknum_);

  double norm = FSI::GenericNormF::compute_norm(*v);

  if (newton() != Teuchos::null)
  {
    newton()->Residual(norm, tolerance());
  }

  return norm;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::FSI::PartialSumNormF::PartialSumNormF(std::string name,
    const Core::LinAlg::MapExtractor& extractor1, double scale1,
    const Core::LinAlg::MapExtractor& extractor2, double scale2,
    Teuchos::RCP<Core::Adapter::CouplingConverter> converter, double tolerance, ScaleType stype)
    : AdaptiveNewtonNormF(name, tolerance, ::NOX::Abstract::Vector::TwoNorm, stype),
      extractor1_(extractor1),
      extractor2_(extractor2),
      scale1_(scale1),
      scale2_(scale2),
      converter_(converter)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double NOX::FSI::PartialSumNormF::compute_norm(const ::NOX::Abstract::Group& grp)
{
  if (!grp.isF()) return -1.0;

  // extract the block epetra vector

  const ::NOX::Abstract::Vector& abstract_f = grp.getF();
  const ::NOX::Epetra::Vector& f = Teuchos::dyn_cast<const ::NOX::Epetra::Vector>(abstract_f);

  // extract the inner vector elements we are interested in

  Teuchos::RCP<Epetra_Vector> v1 = extractor1_.ExtractCondVector(f.getEpetraVector());
  Teuchos::RCP<Epetra_Vector> v2 = extractor2_.ExtractCondVector(f.getEpetraVector());

  Teuchos::RCP<Epetra_Vector> v = converter_->SrcToDst(v2);
  v->Update(scale1_, *v1, scale2_);

  double norm = FSI::GenericNormF::compute_norm(*v);

  if (newton() != Teuchos::null)
  {
    newton()->Residual(norm, tolerance());
  }

  return norm;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::FSI::GenericNormUpdate::GenericNormUpdate(
    std::string name, double tol, ::NOX::Abstract::Vector::NormType ntype, ScaleType stype)
    : status_(::NOX::StatusTest::Unevaluated),
      norm_type_(ntype),
      scale_type_(stype),
      tolerance_(tol),
      norm_update_(0.0),
      name_(name)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::FSI::GenericNormUpdate::GenericNormUpdate(std::string name, double tol, ScaleType stype)
    : status_(::NOX::StatusTest::Unevaluated),
      norm_type_(::NOX::Abstract::Vector::TwoNorm),
      scale_type_(stype),
      tolerance_(tol),
      norm_update_(0.0),
      name_(name)
{
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
::NOX::StatusTest::StatusType NOX::FSI::GenericNormUpdate::checkStatus(
    const ::NOX::Solver::Generic& problem, ::NOX::StatusTest::CheckType checkType)
{
  if (checkType == ::NOX::StatusTest::None)
  {
    status_ = ::NOX::StatusTest::Unevaluated;
    norm_update_ = -1.0;
    return status_;
  }

  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0)
  {
    status_ = ::NOX::StatusTest::Unconverged;
    norm_update_ = -1.0;
    return status_;
  }

  // Check that F exists!
  if (!problem.getSolutionGroup().isF())
  {
    status_ = ::NOX::StatusTest::Unconverged;
    norm_update_ = -1.0;
    return status_;
  }

  const ::NOX::Abstract::Vector& oldSoln = problem.getPreviousSolutionGroup().getX();
  const ::NOX::Abstract::Vector& curSoln = problem.getSolutionGroup().getX();

  if (Teuchos::is_null(update_vector_ptr_)) update_vector_ptr_ = curSoln.clone();

  update_vector_ptr_->update(1.0, curSoln, -1.0, oldSoln, 0.0);

  compute_norm(
      Teuchos::rcp_dynamic_cast<::NOX::Epetra::Vector>(update_vector_ptr_)->getEpetraVector());

  status_ =
      (norm_update_ < tolerance_) ? ::NOX::StatusTest::Converged : ::NOX::StatusTest::Unconverged;
  return status_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double NOX::FSI::GenericNormUpdate::compute_norm(const Epetra_Vector& v)
{
  ::NOX::Epetra::Vector vec(v);
  int n = (scale_type_ == Scaled) ? vec.length() : 0;

  switch (norm_type_)
  {
    case ::NOX::Abstract::Vector::TwoNorm:
      norm_update_ = vec.norm();
      if (scale_type_ == Scaled) norm_update_ /= sqrt(1.0 * n);
      break;

    case ::NOX::Abstract::Vector::OneNorm:
      norm_update_ = vec.norm(norm_type_);
      if (scale_type_ == Scaled) norm_update_ /= n;
      break;

    case ::NOX::Abstract::Vector::MaxNorm:
      norm_update_ = vec.norm(norm_type_);
      if (scale_type_ == Scaled)
        FOUR_C_THROW("It does not make sense to scale a MaxNorm by the vector length.");
      break;

    default:
      FOUR_C_THROW("norm type confusion");
      break;
  }

  return norm_update_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
::NOX::StatusTest::StatusType NOX::FSI::GenericNormUpdate::getStatus() const { return status_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& NOX::FSI::GenericNormUpdate::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j++) stream << ' ';

  stream << status_  // test status
         << name_    // what is tested?
         << " ";

  // check which norm is used and print its name
  if (norm_type_ == ::NOX::Abstract::Vector::TwoNorm)
    stream << "Two-Norm";
  else if (norm_type_ == ::NOX::Abstract::Vector::OneNorm)
    stream << "One-Norm";
  else if (norm_type_ == ::NOX::Abstract::Vector::MaxNorm)
    stream << "Max-Norm";

  // print current value of norm and given tolerance
  stream << " = " << ::NOX::Utils::sciformat(norm_update_, 3) << " < "
         << ::NOX::Utils::sciformat(tolerance_, 3) << "\n";

  // Note: All norms are hard-coded absolute norms. So, we do not neet to print
  // this.                                                      mayt.mt 01/2012

  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double NOX::FSI::GenericNormUpdate::getNormUpdate() const { return norm_update_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double NOX::FSI::GenericNormUpdate::getTolerance() const { return tolerance_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::FSI::PartialNormUpdate::PartialNormUpdate(std::string name,
    const Core::LinAlg::MultiMapExtractor& extractor, int blocknum, double tolerance,
    ::NOX::Abstract::Vector::NormType ntype, ScaleType stype)
    : GenericNormUpdate(name, tolerance, ntype, stype), extractor_(extractor), blocknum_(blocknum)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::FSI::PartialNormUpdate::PartialNormUpdate(std::string name,
    const Core::LinAlg::MultiMapExtractor& extractor, int blocknum, double tolerance,
    ScaleType stype)
    : GenericNormUpdate(name, tolerance, stype), extractor_(extractor), blocknum_(blocknum)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double NOX::FSI::PartialNormUpdate::compute_norm(const Epetra_Vector& v)
{
  return FSI::GenericNormUpdate::compute_norm(*extractor_.ExtractVector(v, blocknum_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::FSI::MinIters::MinIters(int minIterations, const ::NOX::Utils* u)
    : miniters_(minIterations), niters_(0), status_(::NOX::StatusTest::Unevaluated)
{
  if (u != nullptr) utils_ = *u;

  if (miniters_ < 1)
  {
    utils_.err() << "::NOX::StatusTest::MinIters - must choose a number greater than zero"
                 << std::endl;
    throw "NOX Error";
  }
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
::NOX::StatusTest::StatusType NOX::FSI::MinIters::checkStatus(
    const ::NOX::Solver::Generic& problem, ::NOX::StatusTest::CheckType checkType)
{
  switch (checkType)
  {
    case ::NOX::StatusTest::Complete:
    case ::NOX::StatusTest::Minimal:
      niters_ = problem.getNumIterations();
      status_ =
          (niters_ < miniters_) ? ::NOX::StatusTest::Unconverged : ::NOX::StatusTest::Converged;
      break;

    case ::NOX::StatusTest::None:
    default:
      niters_ = -1;
      status_ = ::NOX::StatusTest::Unevaluated;
      break;
  }

  return status_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
::NOX::StatusTest::StatusType NOX::FSI::MinIters::getStatus() const { return status_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& NOX::FSI::MinIters::print(std::ostream& stream, int indent) const { return stream; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int NOX::FSI::MinIters::getMinIters() const { return miniters_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int NOX::FSI::MinIters::getNumIters() const { return niters_; }

FOUR_C_NAMESPACE_CLOSE
