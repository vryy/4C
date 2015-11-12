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
#include "nox_nln_globaldata.H"
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
    Group(Teuchos::ParameterList& printParams,  //!< printing parameters
      Teuchos::ParameterList& grpOptionParams,
      const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,  //!< basically the NOXified time integrator
      const NOX::Epetra::Vector& x,  //!< current solution vector
      const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys,  //!< linear system, matrix and RHS etc.
      const std::map<enum NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >& iConstr  //!< constraint interfaces
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


    //! Returns the interface map
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >&
        GetConstrInterfaces() const;

    //! Returns a pointer to the given soltype. If the solution type is not found an error is thrown.
    Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> GetConstraintInterfacePtr(
        const NOX::NLN::SolutionType& soltype) const;

    //! If the \c errflag is set to true, a error is thrown as soon as we cannot find the corresponding entry in the stl_map. Otherwise a Teuchos::null pointer is returned.
    Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> GetConstraintInterfacePtr(
        const NOX::NLN::SolutionType& soltype, const bool& errflag) const;

    // @name "Get" functions
    //@{
    //! Returns the right-hand-side norms of the primary and constraint quantities
    virtual Teuchos::RCP<const std::vector<double> > GetRHSNorms(
        const std::vector<NOX::Abstract::Vector::NormType>& type,
        const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
        const Teuchos::RCP<const std::vector<NOX::StatusTest::NormF::ScaleType> > scale = Teuchos::null) const;

    //! Returns the root mean square norm of the primary and Lagrange multiplier updates
    virtual Teuchos::RCP<std::vector<double> > GetSolutionUpdateRMS(
        const std::vector<double>& aTol,
        const std::vector<double>& rTol,
        const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
        const std::vector<bool>& disable_implicit_weighting) const;

    //! Returns the status of the active set (only needed for inequality constraint problems)
    enum NOX::StatusTest::StatusType GetActiveSetStatus(
        const enum NOX::NLN::StatusTest::QuantityType& qtype,
        int& activesetsize,
        int& cyclesize) const;

    //@}

  private:
    //! throw Nox error
    void throwError(const std::string& functionName,
        const std::string& errorMsg) const;

  private:
    // constraint interface map
    std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> > userConstraintInterfaces_;
  }; // class Group
} // end namespace CONSTRAINT
} // end namespace  STR
} // end namespace  NOX


#endif /* NOX_NLN_CONSTRAINT_GROUP_H_ */
