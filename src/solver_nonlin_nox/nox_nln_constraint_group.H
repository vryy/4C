/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN implementation of a %NOX::Epetra::Group
       to handle constrained problems.

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_CONSTRAINT_GROUP_H_
#define NOX_NLN_CONSTRAINT_GROUP_H_

#include "nox_nln_group.H"
#include "nox_nln_constraint_interface_required.H"
#include "nox_nln_enum_lists.H"
#include "../drt_lib/drt_dserror.H"

// forward declaration
namespace STR
{
  class TimIntImpl;
}

namespace NOX
{
  namespace NLN
  {
    namespace CONSTRAINT
    {
      class Group : public virtual NOX::NLN::Group
      {
       public:
        //! Standard constructor
        Group(Teuchos::ParameterList& printParams,  //!< printing parameters
            Teuchos::ParameterList& grpOptionParams,
            const Teuchos::RCP<NOX::Epetra::Interface::Required>&
                i,                         //!< basically the NOXified time integrator
            const NOX::Epetra::Vector& x,  //!< current solution vector
            const Teuchos::RCP<NOX::Epetra::LinearSystem>&
                linSys,  //!< linear system, matrix and RHS etc.
            const std::map<enum NOX::NLN::SolutionType,
                Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required>>&
                iConstr  //!< constraint interfaces
        );

        /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
          valid shared linear system. */
        Group(const NOX::NLN::CONSTRAINT::Group& source, NOX::CopyType type = NOX::DeepCopy);

        //! destructor
        virtual ~Group(){};

        //! generate a clone of the given object concerning the given \c CopyType
        Teuchos::RCP<NOX::Abstract::Group> clone(CopyType type) const;

        NOX::Abstract::Group& operator=(const NOX::Epetra::Group& source) override;

        //! Returns the interface map
        const ReqInterfaceMap& GetConstrInterfaces() const;

        //! Returns a pointer to the given soltype. If the solution type is not found an error is
        //! thrown.
        Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> GetConstraintInterfacePtr(
            const NOX::NLN::SolutionType& soltype) const;

        //! If the \c errflag is set to true, a error is thrown as soon as we cannot find the
        //! corresponding entry in the stl_map. Otherwise a Teuchos::null pointer is returned.
        Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> GetConstraintInterfacePtr(
            const NOX::NLN::SolutionType& soltype, const bool& errflag) const;

        // @name "Get" functions
        //@{

        /// derived
        double GetModelValue(
            const enum NOX::NLN::MeritFunction::MeritFctName merit_func_type) const override;

        /// derived
        double GetLinearizedModelTerms(const NOX::Abstract::Vector& dir,
            const enum NOX::NLN::MeritFunction::MeritFctName mf_type,
            const enum NOX::NLN::MeritFunction::LinOrder linorder,
            const enum NOX::NLN::MeritFunction::LinType lintype) const override;

        //! Returns the right-hand-side norms of the primary and constraint quantities
        Teuchos::RCP<const std::vector<double>> GetRHSNorms(
            const std::vector<NOX::Abstract::Vector::NormType>& type,
            const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
            const Teuchos::RCP<const std::vector<NOX::StatusTest::NormF::ScaleType>> scale =
                Teuchos::null) const override;

        //! Returns the root mean square norm of the primary and Lagrange multiplier updates
        Teuchos::RCP<std::vector<double>> GetSolutionUpdateRMS(const NOX::Abstract::Vector& xOld,
            const std::vector<double>& aTol, const std::vector<double>& rTol,
            const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
            const std::vector<bool>& disable_implicit_weighting) const override;

        //! Returns the desired norm of the primary solution updates and Lagrange multiplier updates
        Teuchos::RCP<std::vector<double>> GetSolutionUpdateNorms(const NOX::Abstract::Vector& xOld,
            const std::vector<NOX::Abstract::Vector::NormType>& type,
            const std::vector<StatusTest::QuantityType>& chQ,
            Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale)
            const override;

        //! Returns the desired norm of the previous primary solution and Lagrange multiplier
        //! solution
        Teuchos::RCP<std::vector<double>> GetPreviousSolutionNorms(
            const NOX::Abstract::Vector& xOld,
            const std::vector<NOX::Abstract::Vector::NormType>& type,
            const std::vector<StatusTest::QuantityType>& chQ,
            Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale)
            const override;
        //! @}

        //! @name Handle active set strategies
        //! @{
        //! Returns the current active set map (only needed for inequality constraint problems)
        Teuchos::RCP<const Epetra_Map> GetCurrentActiveSetMap(
            const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

        //! Returns the active set map of the previous Newton step (only needed for inequality
        //! constraint problems)
        Teuchos::RCP<const Epetra_Map> GetOldActiveSetMap(
            const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

        //! Returns basic information about the active set status (no Epetra_Maps needed!)
        enum NOX::StatusTest::StatusType GetActiveSetInfo(
            const enum NOX::NLN::StatusTest::QuantityType& qtype, int& activeset_size) const;

        //@}

       private:
        //! throw Nox error
        void throwError(const std::string& functionName, const std::string& errorMsg) const;

       private:
        // constraint interface map
        ReqInterfaceMap userConstraintInterfaces_;
      };  // class Group
    }     // end namespace CONSTRAINT
  }       // namespace NLN
}  // end namespace  NOX


#endif /* NOX_NLN_CONSTRAINT_GROUP_H_ */
