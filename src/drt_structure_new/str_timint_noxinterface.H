/*-----------------------------------------------------------*/
/*!
\file str_timint_noxinterface.H

\brief Concrete implementation of the Jacobian, Required and
       Preconditioner %NOX::NLN interfaces.

\maintainer Michael Hiermeier

\date Nov 27, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_TIMINT_NOXINTERFACE_H_
#define SRC_DRT_STRUCTURE_NEW_STR_TIMINT_NOXINTERFACE_H_

#include "../solver_nonlin_nox/nox_nln_interface_required.H"  // (1) base class: rhs, status tests and more
#include "../solver_nonlin_nox/nox_nln_interface_jacobian.H"  // (2) base class: jacobian
#include <NOX_Epetra_Interface_Preconditioner.H>  // (3) base class: preconditioner stuff

#include "../solver_nonlin_nox/nox_nln_enum_lists.H"
// forward declaration ...
namespace NOX
{
  namespace NLN
  {
    enum class CorrectionType : int;
  }  // namespace NLN
}  // namespace NOX
namespace LINALG
{
  class SparseOperator;
  class SparseMatrix;
}  // namespace LINALG
namespace INPAR
{
  namespace STR
  {
    enum ModelType : int;
  }  // namespace STR
}  // namespace INPAR
namespace STR
{
  class Dbc;
  namespace IMPLICIT
  {
    class Generic;
  }  // namespace IMPLICIT
  namespace TIMINT
  {
    class Base;
    class BaseDataGlobalState;
    class NoxInterface : virtual public NOX::NLN::Interface::Required,
                         virtual public NOX::NLN::Interface::Jacobian,
                         virtual public NOX::Epetra::Interface::Preconditioner
    {
     public:
      //! constructor
      NoxInterface();

      //! destructor
      virtual ~NoxInterface(){};

      //! Init function
      virtual void Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
          const Teuchos::RCP<STR::IMPLICIT::Generic>& implint_ptr,
          const Teuchos::RCP<STR::Dbc>& dbc_ptr,
          const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr);

      virtual void Setup();

      //!@{
      /*! compute the right hand side entries
       *  (derived from NOX::Epetra::Interface::Required) */
      virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag);

      /*! compute jacobian
       *  ( derived from NOX::Epetra::Inteface::Jacobian) */
      virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

      /*! compute right hand side and jacobian
       *  (derived from NOX::NLN::Interface::Jacobian) */
      virtual bool computeFandJacobian(
          const Epetra_Vector& x, Epetra_Vector& rhs, Epetra_Operator& jac);

      virtual bool computeCorrectionSystem(const enum NOX::NLN::CorrectionType type,
          const NOX::Abstract::Group& grp, const Epetra_Vector& x, Epetra_Vector& rhs,
          Epetra_Operator& jac);

      /*! compute preconditioner
       *  (derived from NOX::Epetra::Interface::Preconditioner) */
      virtual bool computePreconditioner(
          const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams = 0);

      /*! Get the norm of right hand side rows/entries related to
       *  primary DoFs (derived from NOX::NLN::Interface::Required) */
      virtual double GetPrimaryRHSNorms(const Epetra_Vector& F,
          const NOX::NLN::StatusTest::QuantityType& checkquantity,
          const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
          const bool& isscaled = false) const;

      /*! Get the root mean square of the solution update (vector) entries
       *  (derived from NOX::NLN::Interface::Required) */
      virtual double GetPrimarySolutionUpdateRMS(const Epetra_Vector& xnew,
          const Epetra_Vector& xold, const double& aTol, const double& rTol,
          const NOX::NLN::StatusTest::QuantityType& checkQuantity,
          const bool& disable_implicit_weighting = false) const;

      /*! Returns the desired norm of the solution update (vector) entries
       *  (derived from NOX::NLN::Interface::Required) */
      virtual double GetPrimarySolutionUpdateNorms(const Epetra_Vector& xnew,
          const Epetra_Vector& xold, const NOX::NLN::StatusTest::QuantityType& checkquantity,
          const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
          const bool& isscaled = false) const;

      /*! Returns the previous solution norm of primary DoF fields
       *  (derived from NOX::NLN::Interface::Required) */
      virtual double GetPreviousPrimarySolutionNorms(const Epetra_Vector& xold,
          const NOX::NLN::StatusTest::QuantityType& checkquantity,
          const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
          const bool& isscaled = false) const;

      /*! Compute and return some energy representative or any other scalar value
       *  which is capable to describe the solution path progress
       *  (derived from NOX::NLN::Interface::Required) */
      virtual double GetModelValue(const Epetra_Vector& x, const Epetra_Vector& F,
          const enum NOX::NLN::MeritFunction::MeritFctName merit_func_type) const;

      virtual double GetLinearizedModelTerms(const NOX::Abstract::Group* group,
          const Epetra_Vector& dir, const enum NOX::NLN::MeritFunction::MeritFctName mf_type,
          const enum NOX::NLN::MeritFunction::LinOrder linorder,
          const enum NOX::NLN::MeritFunction::LinType lintype) const;

      /*! \brief calculate characteristic/reference norms for forces
       *
       *  Necessary for the LinearSystem objects.
       *  (derived from NOX::NLN::Interface::Required) */
      virtual double CalcRefNormForce();

      virtual void CreateBackupState(const Epetra_Vector& dir);

      virtual void RecoverFromBackupState();

      //!@}

      // Get element based scaling operator
      virtual Teuchos::RCP<LINALG::SparseMatrix> CalcJacobianContributionsFromElementLevelForPTC();

      //! Access the implicit integrator
      STR::IMPLICIT::Generic& ImplInt();

     protected:
      //! Returns the init state
      inline const bool& IsInit() const { return isinit_; };

      //! Returns the setup state
      inline const bool& IsSetup() const { return issetup_; };

      //! check if init has been called
      void CheckInit() const;

      //! check if init and setup have been called
      void CheckInitSetup() const;

      double GetLinearizedEnergyModelTerms(const NOX::Abstract::Group* group,
          const Epetra_Vector& dir, const enum NOX::NLN::MeritFunction::LinOrder linorder,
          const enum NOX::NLN::MeritFunction::LinType lintype) const;

      void FindConstraintModels(const NOX::Abstract::Group* grp,
          std::vector<INPAR::STR::ModelType>& constraint_models) const;

      //! calculate norm in Get*Norms functions
      double CalculateNorm(Teuchos::RCP<Epetra_Vector> quantity,
          const NOX::Abstract::Vector::NormType type, const bool isscaled) const;

     protected:
      //! init flag
      bool isinit_;

      //! setup flag
      bool issetup_;

     private:
      //! global state data container
      Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

      Teuchos::RCP<const STR::TIMINT::Base> timint_ptr_;

      Teuchos::RCP<STR::IMPLICIT::Generic> implint_ptr_;

      Teuchos::RCP<STR::Dbc> dbc_ptr_;
    };  // class NoxInterface
  }     // namespace TIMINT
}  // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_TIMINT_NOXINTERFACE_H_ */
