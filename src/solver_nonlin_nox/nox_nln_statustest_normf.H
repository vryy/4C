/*-----------------------------------------------------------*/
/*!
\file nox_nln_statustest_normf.H

\maintainer Michael Hiermeier

\date Aug 3, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_STATUSTEST_NORMF_H_
#define NOX_NLN_STATUSTEST_NORMF_H_

#include <NOX_StatusTest_NormF.H>   // (derived) base class

#include "nox_nln_statustest_factory.H"

namespace NOX {
namespace NLN {
namespace StatusTest {
  //! Various convergence tests based on the norm of the residual.
  /*!

   This Status Test class is closely related to the NOX_StatusTest_NormF class.
   Because of some incompatibilities and necessary generalizations, we restate
   this class here. The main difference is, that the new class is capable to
   handle a std_vector of right-hand-side norms. In this way it's easily possible
   to check the norm of different quantities at once in a very general way.

   The norms are calculated in the corresponding required interfaces.
   */

  class NormF : public NOX::StatusTest::Generic
  {
  public:

    //! Constructor
    NormF(
      const std::vector<NOX::NLN::StatusTest::QuantityType>& checkList,
      const std::vector<NOX::StatusTest::NormF::ToleranceType>& toltype,
      const std::vector<double>& tolerance,
      const std::vector<NOX::Abstract::Vector::NormType>& ntype,
      const std::vector<NOX::StatusTest::NormF::ScaleType>& stype,
      const NOX::Utils* u = NULL);

    //! Destructor.
    virtual ~NormF() {};

    /* @name Accessor Functions
       Used to query current values of variables in the status test.
     */
    //@{
    //! Returns the value of the F-norm of the corresponding quantity computed in the last call to checkStatus.
    virtual double GetNormF(const NOX::NLN::StatusTest::QuantityType& qType) const;

    //! Returns the true tolerance of the corresponding quantity.
    virtual double GetTrueTolerance(const NOX::NLN::StatusTest::QuantityType& qType) const;

    //! Returns the specified tolerance set in the constructor for the corresponding quantity.
    virtual double GetSpecifiedTolerance(const NOX::NLN::StatusTest::QuantityType& qType) const;

    //! Returns the initial tolerance of the corresponding quantity.
    virtual double GetInitialTolerance(const NOX::NLN::StatusTest::QuantityType& qType) const;

    /*! \brief Returns the norm type as <int> of the desired quantity
     *
     *  If the given quantity cannot be found a default value of -100 is returned. */
    int GetNormType(const NOX::NLN::StatusTest::QuantityType& qType) const;

    //! Check for the given quantity
    bool IsQuantity(const NOX::NLN::StatusTest::QuantityType& qType) const;
    //@}

    //! derived
    virtual NOX::StatusTest::StatusType checkStatus(
        const NOX::Solver::Generic& problem,
        NOX::StatusTest::CheckType checkType);

    //! derived
    //! NOTE: returns the global status of all normF tests
    virtual NOX::StatusTest::StatusType getStatus() const;

    //! derived
    virtual std::ostream& print(std::ostream& stream, int indent = 0) const;

  protected:
    /*! In the case of a relative norm calculation, initializes
      \c trueTolerance based on the F-value at the initial guess.*/
    void relativeSetup(Teuchos::RCP<const NOX::Abstract::Group>& initialGuess);

    /*! \brief Calculate the norm of the specified quantities
      of the rhs for the given group according to the scaling
      type, norm type, and tolerance type.

      \note Returns Teuchos::null if F(x) has not been calculated for the given
      grp (i.e., grp.isF() is false).
    */
    Teuchos::RCP<const std::vector<double> > ComputeNorm(Teuchos::RCP<const NOX::Abstract::Group>& initialGuess);

  protected:
    //! number of status tests
    const std::size_t nChecks_;

    //! enums of the quantities we want to check
    std::vector<NOX::NLN::StatusTest::QuantityType> checkList_;

    //! global status
    NOX::StatusTest::StatusType gStatus_;

    //! local/entry-wise status
    std::vector<NOX::StatusTest::StatusType> status_;

    //! Type of norm to use
    std::vector<NOX::Abstract::Vector::NormType> normType_;

    //! Scaling to use
    std::vector<NOX::StatusTest::NormF::ScaleType> scaleType_;

    //! Tolerance type (i.e., relative or absolute)
    std::vector<NOX::StatusTest::NormF::ToleranceType> toleranceType_;

    //! Tolerance required for convergence.
    std::vector<double> specifiedTolerance_;

    //! Initial tolerance
    Teuchos::RCP<const std::vector<double> > initialTolerance_;

    //! True tolerance value, i.e., specifiedTolerance * initialTolerance
    std::vector<double> trueTolerance_;

    //! Norm of F to be compared to trueTolerance
    Teuchos::RCP<const std::vector<double> > normF_;

    //! Ostream used to print errors
    NOX::Utils utils_;
  };  // class NormF
} // namespace StatusTest
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_STATUSTEST_NORMF_H_ */
