/*-----------------------------------------------------------*/
/*!
\file nox_nln_group.H

\maintainer Michael Hiermeier

\date Jun 29, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_GROUP_H_
#define NOX_NLN_GROUP_H_

/*----------------------------------------------------------------------------*/
/* headers */
#include <NOX_Epetra_Group.H>   // base class
#include <NOX_StatusTest_NormF.H>

#include "nox_nln_statustest_factory.H"

// forward declarations
namespace NOX {
namespace NLN {
namespace Interface {
  class Required;
} // namespace Interface

  class Group : public NOX::Epetra::Group
  {
  public:

    //! Standard Constructor
    Group
    (
      Teuchos::ParameterList& printParams,  //!< printing parameters
      const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,  //!< basically the NOXified user interface
      const NOX::Epetra::Vector& x,  //!< current solution vector
      const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys  //!< linear system, matrix and RHS etc.
    );

    /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
      valid shared linear system. */
    Group(const NOX::NLN::Group& source,
      NOX::CopyType type = NOX::DeepCopy);

    //! Destructor
    virtual ~Group(){};

    //! generate a clone of the given object concerning the given \c CopyType
    virtual Teuchos::RCP<NOX::Abstract::Group> clone(CopyType type = DeepCopy) const;

    //! compute/update the current state variables
    virtual void computeX(
        const NOX::NLN::Group& grp,
        const NOX::Epetra::Vector& d,
        double step);
    virtual void computeX(const NOX::Abstract::Group& grp,
              const NOX::Abstract::Vector& d,
              double step);

    //! Compute and store \f$F(x)\f$ and the jacobian \f$\frac{\partial F(x)}{\partial x}\f$ at the same time.
    //! This can result in a huge performance gain in some special cases, e.g. contact problems.
    virtual NOX::Abstract::Group::ReturnType computeFandJacobian();

    //! set right hand side
    NOX::Abstract::Group::ReturnType setF(const Teuchos::RCP<const NOX::Epetra::Vector> Fptr);

    //! set jacobian operator
    virtual NOX::Abstract::Group::ReturnType setJacobianOperator(const Teuchos::RCP<const Epetra_Operator> jacOperator);

    //! returns the nox_nln_interface_required pointer
    Teuchos::RCP<const NOX::NLN::Interface::Required> GetNlnReqInterfacePtr() const;

    //! returns the primary rhs norms
    virtual Teuchos::RCP<const std::vector<double> > GetRHSNorms(
        const std::vector<NOX::Abstract::Vector::NormType>& type,
        const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
        Teuchos::RCP<const std::vector<NOX::StatusTest::NormF::ScaleType> > scale = Teuchos::null) const;

    //! returns the Root Mean Squares (abbr.: RMS) of the primary solution updates
    virtual Teuchos::RCP<std::vector<double> > GetSolutionUpdateRMS(
        const std::vector<double>& aTol,
        const std::vector<double>& rTol,
        const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
        const std::vector<bool>& disable_implicit_weighting) const;

    //! compute and return some energy representative
    virtual const double GetObjectiveModelValue(const std::string& name) const;

  private:
    //! Throw an NOX_error
    void throwError(const std::string& functionName,const std::string& errorMsg) const;
  };
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_GROUP_H_ */
