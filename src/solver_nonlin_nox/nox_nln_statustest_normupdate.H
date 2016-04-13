/*-----------------------------------------------------------*/
/*!
\file nox_nln_statustest_normupdate.H

\maintainer Michael Hiermeier

\date Sep 17, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_STATUSTEST_NORMUPDATE_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_STATUSTEST_NORMUPDATE_H_


#include <NOX_StatusTest_Generic.H>

#include <NOX_Abstract_Vector.H>
#include <NOX_Utils.H>
#include "nox_nln_statustest_factory.H"

namespace NOX {
namespace Abstract {
  class Group;
} // namespace Abstract
namespace NLN {
namespace StatusTest {
  //! Absolute and relative convergence tests based on the norm of the increment.
  /*!
   I recommend to use the weighted root mean square test (normwrms), instead or additionally
   (AND combination)!

   Nevertheless, I want to give you the possibility to check the absolute and relative norm
   of the solution increment corresponding to the different quantities as well. This class is
   leaned on the NOX::NLN::StatusTest::NormF and NOX::NLN::StatusTest::NormWRMS test and has
   an equivalent structure. It is capable to handle a std::vector of increment norms. In this
   way it's easily possible to check the norms of different quantities at once in a very
   general way.

   One difference to the NormF case is the calculation of the reference value. We use
   the following equation:

             $\frac{\| x_i^{(k+1)} - x_i^{(k)} \|}{\|x_i^{(k)}\|} \leq RTOL$

   */
  class NormUpdate : public NOX::StatusTest::Generic
  {
  public:
    //! Type that determines whether to scale the norm by the problem size.
    enum ScaleType {
      //! No norm scaling
      Unscaled,
      //! Scale the norm by the length of the vector
      Scaled
    };

    //! Type that determines whether the check is absolute or relative
    enum ToleranceType {
      //! Relative convergence check
      Relative,
      //! Absolute
      Absolute
    };

  public:
    //! Constructor.
    NormUpdate(
        const std::vector<NOX::NLN::StatusTest::QuantityType>& checkList,
        const std::vector<NormUpdate::ToleranceType>& toltype,
        const std::vector<double>& tolerance,
        const std::vector<NOX::Abstract::Vector::NormType>& ntype,
        const std::vector<NormUpdate::ScaleType>& stype,
        const double& alpha,
        const double& beta,
        const NOX::Utils* u = NULL);

    //! Destructor.
    virtual ~NormUpdate() {};

    //! derived
    virtual NOX::StatusTest::StatusType checkStatus(
        const NOX::Solver::Generic& problem,
        NOX::StatusTest::CheckType checkType);

    /*! \brief Returns the norm type as <int> of the desired quantity
     *
     *  If the given quantity cannot be found a default value of -100 is returned. */
    int GetNormType(const NOX::NLN::StatusTest::QuantityType& qType) const;

    //! Check for the given quantity
    bool IsQuantity(const NOX::NLN::StatusTest::QuantityType& qType) const;

    //! derived
    //! NOTE: returns the global status of all normF tests
    virtual NOX::StatusTest::StatusType getStatus() const;

    //! derived
    virtual std::ostream& print(std::ostream& stream, int indent = 0) const;
  protected:
    virtual void ComputeNorm(
        const NOX::Abstract::Group& grp,
        const NOX::Solver::Generic& problem);

  private:
    //! throws an NOX error
    void throwError(
        const std::string& functionName,
        const std::string& errorMsg) const;

  protected:
    //! number of status tests
    const std::size_t nChecks_;

    //! computed step size
    double computedStepSize_;

    //! Minimum step size allowed during a line search for increment norm to be flagged as converged.
    double alpha_;

    //! Actual tolerance achieved by the linear solver during the last linear solve.
    double achievedTol_;

    //! Maximum linear solve tolerance allowed for incr norm to be flagged as converged.
    double beta_;

    //! enums of the quantities we want check
    std::vector<NOX::NLN::StatusTest::QuantityType> checkList_;

    //! global status
    NOX::StatusTest::StatusType gStatus_;

    //! Status
    std::vector<NOX::StatusTest::StatusType> status_;

    //! Type of norm to use
    std::vector<NOX::Abstract::Vector::NormType> normType_;

    //! Scaling to use
    std::vector<NormUpdate::ScaleType> scaleType_;

    //! Tolerance type (i.e., relative or absolute)
    std::vector<NormUpdate::ToleranceType> toleranceType_;

    //! Tolerance required for convergence.
    std::vector<double> specifiedTolerance_;

    //! reference norm of the solution vector (relative only)
    Teuchos::RCP<const std::vector<double> > normRefSol_;

    //! True tolerance value, i.e., specifiedTolerance * initialTolerance
    std::vector<double> trueTolerance_;

    //! Norm of F to be compared to trueTolerance
    Teuchos::RCP<const std::vector<double> > normUpdate_;

    //! Ostream used to print errors
    NOX::Utils utils_;

    //! Flag that tells the print method whether to print the criteria 2 information.
    bool printCriteria2Info_;

    //! Flag that tells the print method whether to print the criteria 3 information.
    bool printCriteria3Info_;


  }; // class NormUpdate
} // namespace StatusTest
} // namespace NLN
} // namespace NOX

#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_STATUSTEST_NORMUPDATE_H_ */
