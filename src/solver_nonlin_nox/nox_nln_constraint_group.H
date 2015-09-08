/*-----------------------------------------------------------*/
/*!
\file nox_nln_constraint_group.H

\maintainer Michael Hiermeier

\date Jun 29, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_CONSTRAINT_GROUP_H_
#define NOX_NLN_CONSTRAINT_GROUP_H_

#include "nox_nln_group.H"
#include "nox_nln_constraint_interface_required.H"
#include "../drt_lib/drt_dserror.H"

// forward declaration
namespace STR {
  class TimIntImpl;
}

namespace NOX {
namespace NLN {
namespace CONSTRAINT {

  class Group : public NOX::NLN::Group
  {
  public:
    //! Standard constructor
    Group
    (
      Teuchos::ParameterList& printParams,  //!< printing parameters
      const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,  //!< basically the NOXified time integrator
      const NOX::Epetra::Vector& x,  //!< current solution vector
      const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys,  //!< linear system, matrix and RHS etc.
      const Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required>& iConstr  //!< constraint interface
    );

    /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
      valid shared linear system. */
    Group(const NOX::NLN::CONSTRAINT::Group& source, NOX::CopyType type = NOX::DeepCopy);

    //! destructor
    virtual ~Group() {};

    //! generate a clone of the given object concerning the given \c CopyType
    Teuchos::RCP<NOX::Abstract::Group> clone(CopyType type) const;

    //! derived
    virtual void computeX(
        const NOX::NLN::CONSTRAINT::Group& grp,
        const NOX::Epetra::Vector& d,
        double step);
    virtual void computeX(
        const NOX::Abstract::Group& grp,
        const NOX::Abstract::Vector& d,
        double step);

    //! Returns the constraint interface
    virtual const NOX::NLN::CONSTRAINT::Interface::Required& GetConstr() const
    {
      if (userConstraintInterfacePtr_.is_null())
        dserror("The userConstraintInterfacePtr_ is not defined!");

      return *userConstraintInterfacePtr_;
    }

    //! Returns the constraint interface pointer
    virtual Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> GetConstrPtr() const
    {
      return userConstraintInterfacePtr_;
    }

    // @name "Get" functions
    //@{
    //! Returns the right-hand-side norms of the primary and constraint quantities
    virtual Teuchos::RCP<const std::vector<double> > GetRHSNorms(
        const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
        const std::vector<NOX::Abstract::Vector::NormType>& type,
        const Teuchos::RCP<const std::vector<NOX::StatusTest::NormF::ScaleType> > scale = Teuchos::null) const;

    //! Returns the root mean square norm of the primary and Lagrange multiplier updates
    virtual Teuchos::RCP<std::vector<double> > GetSolutionUpdateRMS(
        const std::vector<double>& aTol,
        const std::vector<double>& rTol,
        const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
        const std::vector<bool>& disable_implicit_weighting) const;

    //@}

  private:
    //! throw Nox error
    void throwError(const std::string& functionName,
        const std::string& errorMsg) const;

  private:
    // constraint interface pointer
    Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> userConstraintInterfacePtr_;

  };
} // end namespace CONSTRAINT
} // end namespace  STR
} // end namespace  NOX


#endif /* NOX_NLN_CONSTRAINT_GROUP_H_ */
