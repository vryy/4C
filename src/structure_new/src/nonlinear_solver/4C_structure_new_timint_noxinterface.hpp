/*-----------------------------------------------------------*/
/*! \file

\brief Concrete implementation of the Jacobian, Required and
       Preconditioner %NOX::NLN interfaces.


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_NOXINTERFACE_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_NOXINTERFACE_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_solver_nonlin_nox_interface_jacobian.hpp"  // (2) base class: jacobian
#include "4C_solver_nonlin_nox_interface_required.hpp"  // (1) base class: rhs, status tests and more

#include <NOX_Epetra_Interface_Preconditioner.H>  // (3) base class: preconditioner stuff

FOUR_C_NAMESPACE_OPEN
// forward declaration ...
namespace NOX
{
  namespace Nln
  {
    enum class CorrectionType : int;
  }  // namespace Nln
}  // namespace NOX
namespace Core::LinAlg
{
  class SparseOperator;
  class SparseMatrix;
}  // namespace Core::LinAlg
namespace Inpar
{
  namespace STR
  {
    enum ModelType : int;
  }  // namespace STR
}  // namespace Inpar
namespace STR
{
  class Dbc;
  class Integrator;
  namespace TimeInt
  {
    class Base;
    class BaseDataGlobalState;
    class NoxInterface : virtual public NOX::Nln::Interface::Required,
                         virtual public NOX::Nln::Interface::Jacobian,
                         virtual public ::NOX::Epetra::Interface::Preconditioner
    {
     public:
      //! constructor
      NoxInterface();

      //! Init function
      virtual void Init(const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& gstate_ptr,
          const Teuchos::RCP<STR::Integrator>& int_ptr, const Teuchos::RCP<STR::Dbc>& dbc_ptr,
          const Teuchos::RCP<const STR::TimeInt::Base>& timint_ptr);

      virtual void Setup();

      //!@{
      /*! compute the right hand side entries
       *  (derived from ::NOX::Epetra::Interface::Required) */
      bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag) override;

      /*! compute jacobian
       *  ( derived from ::NOX::Epetra::Inteface::Jacobian) */
      bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) override;

      /*! compute right hand side and jacobian
       *  (derived from NOX::Nln::Interface::Jacobian) */
      bool computeFandJacobian(
          const Epetra_Vector& x, Epetra_Vector& rhs, Epetra_Operator& jac) override;

      bool compute_correction_system(const enum NOX::Nln::CorrectionType type,
          const ::NOX::Abstract::Group& grp, const Epetra_Vector& x, Epetra_Vector& rhs,
          Epetra_Operator& jac) override;

      /*! compute preconditioner
       *  (derived from ::NOX::Epetra::Interface::Preconditioner) */
      bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M,
          Teuchos::ParameterList* precParams = nullptr) override;

      /*! Get the norm of right hand side rows/entries related to
       *  primary DoFs (derived from NOX::Nln::Interface::Required) */
      double get_primary_rhs_norms(const Epetra_Vector& F,
          const NOX::Nln::StatusTest::QuantityType& checkquantity,
          const ::NOX::Abstract::Vector::NormType& type = ::NOX::Abstract::Vector::TwoNorm,
          const bool& isscaled = false) const override;

      /*! Get the root mean square of the solution update (vector) entries
       *  (derived from NOX::Nln::Interface::Required) */
      double get_primary_solution_update_rms(const Epetra_Vector& xnew, const Epetra_Vector& xold,
          const double& aTol, const double& rTol,
          const NOX::Nln::StatusTest::QuantityType& checkQuantity,
          const bool& disable_implicit_weighting = false) const override;

      /*! Returns the desired norm of the solution update (vector) entries
       *  (derived from NOX::Nln::Interface::Required) */
      double get_primary_solution_update_norms(const Epetra_Vector& xnew, const Epetra_Vector& xold,
          const NOX::Nln::StatusTest::QuantityType& checkquantity,
          const ::NOX::Abstract::Vector::NormType& type = ::NOX::Abstract::Vector::TwoNorm,
          const bool& isscaled = false) const override;

      /*! Returns the previous solution norm of primary DoF fields
       *  (derived from NOX::Nln::Interface::Required) */
      double get_previous_primary_solution_norms(const Epetra_Vector& xold,
          const NOX::Nln::StatusTest::QuantityType& checkquantity,
          const ::NOX::Abstract::Vector::NormType& type = ::NOX::Abstract::Vector::TwoNorm,
          const bool& isscaled = false) const override;

      /*! Compute and return some energy representative or any other scalar value
       *  which is capable to describe the solution path progress
       *  (derived from NOX::Nln::Interface::Required) */
      double get_model_value(const Epetra_Vector& x, const Epetra_Vector& F,
          const NOX::Nln::MeritFunction::MeritFctName merit_func_type) const override;

      double get_linearized_model_terms(const ::NOX::Abstract::Group* group,
          const Epetra_Vector& dir, const enum NOX::Nln::MeritFunction::MeritFctName mf_type,
          const enum NOX::Nln::MeritFunction::LinOrder linorder,
          const enum NOX::Nln::MeritFunction::LinType lintype) const override;

      /*! \brief calculate characteristic/reference norms for forces
       *
       *  Necessary for the LinearSystem objects.
       *  (derived from NOX::Nln::Interface::Required) */
      double calc_ref_norm_force() override;

      /// create back-up state of condensed solution variables (e.g. EAS)
      void create_backup_state(const Epetra_Vector& dir) override;

      /// recover from back-up
      void recover_from_backup_state() override;

      /// compute the current volumes for all elements
      bool compute_element_volumes(
          const Epetra_Vector& x, Teuchos::RCP<Epetra_Vector>& ele_vols) const override;

      /// fill the sets with DOFs of the desired elements
      void get_dofs_from_elements(
          const std::vector<int>& my_ele_gids, std::set<int>& my_ele_dofs) const override;

      //!@}

      // Get element based scaling operator
      Teuchos::RCP<Core::LinAlg::SparseMatrix>
      calc_jacobian_contributions_from_element_level_for_ptc() override;

      //! Access the implicit integrator
      STR::Integrator& ImplInt();

     protected:
      //! Returns the init state
      inline const bool& is_init() const { return isinit_; };

      //! Returns the setup state
      inline const bool& is_setup() const { return issetup_; };

      //! check if init has been called
      void check_init() const;

      //! check if init and setup have been called
      void check_init_setup() const;

      double get_linearized_energy_model_terms(const ::NOX::Abstract::Group* group,
          const Epetra_Vector& dir, const enum NOX::Nln::MeritFunction::LinOrder linorder,
          const enum NOX::Nln::MeritFunction::LinType lintype) const;

      void find_constraint_models(const ::NOX::Abstract::Group* grp,
          std::vector<Inpar::STR::ModelType>& constraint_models) const;

      //! calculate norm in Get*Norms functions
      double calculate_norm(Teuchos::RCP<Epetra_Vector> quantity,
          const ::NOX::Abstract::Vector::NormType type, const bool isscaled) const;

     protected:
      //! init flag
      bool isinit_;

      //! setup flag
      bool issetup_;

     private:
      //! global state data container
      Teuchos::RCP<STR::TimeInt::BaseDataGlobalState> gstate_ptr_;

      Teuchos::RCP<const STR::TimeInt::Base> timint_ptr_;

      Teuchos::RCP<STR::Integrator> int_ptr_;

      Teuchos::RCP<STR::Dbc> dbc_ptr_;
    };  // class nox_interface
  }     // namespace TimeInt
}  // namespace STR


FOUR_C_NAMESPACE_CLOSE

#endif
