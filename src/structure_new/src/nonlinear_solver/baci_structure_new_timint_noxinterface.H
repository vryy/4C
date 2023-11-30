/*-----------------------------------------------------------*/
/*! \file

\brief Concrete implementation of the Jacobian, Required and
       Preconditioner %NOX::NLN interfaces.


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BACI_STRUCTURE_NEW_TIMINT_NOXINTERFACE_H
#define BACI_STRUCTURE_NEW_TIMINT_NOXINTERFACE_H

#include "baci_solver_nonlin_nox_enum_lists.H"
#include "baci_solver_nonlin_nox_interface_jacobian.H"  // (2) base class: jacobian
#include "baci_solver_nonlin_nox_interface_required.H"  // (1) base class: rhs, status tests and more

#include <NOX_Epetra_Interface_Preconditioner.H>  // (3) base class: preconditioner stuff
// forward declaration ...
namespace NOX
{
  namespace NLN
  {
    enum class CorrectionType : int;
  }  // namespace NLN
}  // namespace NOX
namespace CORE::LINALG
{
  class SparseOperator;
  class SparseMatrix;
}  // namespace CORE::LINALG
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
  class Integrator;
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


      //! Init function
      virtual void Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
          const Teuchos::RCP<STR::Integrator>& int_ptr, const Teuchos::RCP<STR::Dbc>& dbc_ptr,
          const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr);

      virtual void Setup();

      //!@{
      /*! compute the right hand side entries
       *  (derived from NOX::Epetra::Interface::Required) */
      bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag) override;

      /*! compute jacobian
       *  ( derived from NOX::Epetra::Inteface::Jacobian) */
      bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) override;

      /*! compute right hand side and jacobian
       *  (derived from NOX::NLN::Interface::Jacobian) */
      bool computeFandJacobian(
          const Epetra_Vector& x, Epetra_Vector& rhs, Epetra_Operator& jac) override;

      bool computeCorrectionSystem(const enum NOX::NLN::CorrectionType type,
          const NOX::Abstract::Group& grp, const Epetra_Vector& x, Epetra_Vector& rhs,
          Epetra_Operator& jac) override;

      /*! compute preconditioner
       *  (derived from NOX::Epetra::Interface::Preconditioner) */
      bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M,
          Teuchos::ParameterList* precParams = nullptr) override;

      /*! Get the norm of right hand side rows/entries related to
       *  primary DoFs (derived from NOX::NLN::Interface::Required) */
      double GetPrimaryRHSNorms(const Epetra_Vector& F,
          const NOX::NLN::StatusTest::QuantityType& checkquantity,
          const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
          const bool& isscaled = false) const override;

      /*! Get the root mean square of the solution update (vector) entries
       *  (derived from NOX::NLN::Interface::Required) */
      double GetPrimarySolutionUpdateRMS(const Epetra_Vector& xnew, const Epetra_Vector& xold,
          const double& aTol, const double& rTol,
          const NOX::NLN::StatusTest::QuantityType& checkQuantity,
          const bool& disable_implicit_weighting = false) const override;

      /*! Returns the desired norm of the solution update (vector) entries
       *  (derived from NOX::NLN::Interface::Required) */
      double GetPrimarySolutionUpdateNorms(const Epetra_Vector& xnew, const Epetra_Vector& xold,
          const NOX::NLN::StatusTest::QuantityType& checkquantity,
          const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
          const bool& isscaled = false) const override;

      /*! Returns the previous solution norm of primary DoF fields
       *  (derived from NOX::NLN::Interface::Required) */
      double GetPreviousPrimarySolutionNorms(const Epetra_Vector& xold,
          const NOX::NLN::StatusTest::QuantityType& checkquantity,
          const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
          const bool& isscaled = false) const override;

      /*! Compute and return some energy representative or any other scalar value
       *  which is capable to describe the solution path progress
       *  (derived from NOX::NLN::Interface::Required) */
      double GetModelValue(const Epetra_Vector& x, const Epetra_Vector& F,
          const enum NOX::NLN::MeritFunction::MeritFctName merit_func_type) const override;

      double GetLinearizedModelTerms(const NOX::Abstract::Group* group, const Epetra_Vector& dir,
          const enum NOX::NLN::MeritFunction::MeritFctName mf_type,
          const enum NOX::NLN::MeritFunction::LinOrder linorder,
          const enum NOX::NLN::MeritFunction::LinType lintype) const override;

      /*! \brief calculate characteristic/reference norms for forces
       *
       *  Necessary for the LinearSystem objects.
       *  (derived from NOX::NLN::Interface::Required) */
      double CalcRefNormForce() override;

      /// create back-up state of condensed solution variables (e.g. EAS)
      void CreateBackupState(const Epetra_Vector& dir) override;

      /// recover from back-up
      void RecoverFromBackupState() override;

      /// compute the current volumes for all elements
      bool computeElementVolumes(
          const Epetra_Vector& x, Teuchos::RCP<Epetra_Vector>& ele_vols) const override;

      /// fill the sets with DOFs of the desired elements
      void getDofsFromElements(
          const std::vector<int>& my_ele_gids, std::set<int>& my_ele_dofs) const override;

      //!@}

      // Get element based scaling operator
      Teuchos::RCP<CORE::LINALG::SparseMatrix> CalcJacobianContributionsFromElementLevelForPTC()
          override;

      //! Access the implicit integrator
      STR::Integrator& ImplInt();

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

      Teuchos::RCP<STR::Integrator> int_ptr_;

      Teuchos::RCP<STR::Dbc> dbc_ptr_;
    };  // class NoxInterface
  }     // namespace TIMINT
}  // namespace STR


#endif  // STRUCTURE_NEW_TIMINT_NOXINTERFACE_H
