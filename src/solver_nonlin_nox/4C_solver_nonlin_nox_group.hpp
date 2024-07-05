/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN implementation of a %::NOX::Epetra::Group
       to handle unconstrained problems.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_GROUP_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_GROUP_HPP

/*----------------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_linalg_serialdensevector.hpp"
#include "4C_solver_nonlin_nox_forward_decl.hpp"
#include "4C_solver_nonlin_nox_statustest_normupdate.hpp"

#include <NOX_Epetra_Group.H>  // base class
#include <NOX_StatusTest_NormF.H>

#include <set>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::LinAlg
{
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace NOX
{
  namespace Nln
  {
    namespace Solver
    {
      class PseudoTransient;
    }  // namespace Solver
    namespace Interface
    {
      class Required;
    }  // namespace Interface
    namespace GROUP
    {
      class PrePostOperator;
    }  // namespace GROUP

    class Group : public virtual ::NOX::Epetra::Group
    {
     public:
      //! Standard Constructor
      Group(Teuchos::ParameterList& printParams,    //!< printing parameters
          Teuchos::ParameterList& grpOptionParams,  //!< group option parameters
          const Teuchos::RCP<::NOX::Epetra::Interface::Required>&
              i,                           //!< basically the NOXified user interface
          const ::NOX::Epetra::Vector& x,  //!< current solution vector
          const Teuchos::RCP<::NOX::Epetra::LinearSystem>&
              linSys  //!< linear system, matrix and RHS etc.
      );

      /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
        valid shared linear system. */
      Group(const NOX::Nln::Group& source, ::NOX::CopyType type = ::NOX::DeepCopy);

      /// assign operator
      ::NOX::Abstract::Group& operator=(const ::NOX::Abstract::Group& source) override;
      ::NOX::Abstract::Group& operator=(const ::NOX::Epetra::Group& source) override;

      Teuchos::RCP<::NOX::Abstract::Group> clone(
          ::NOX::CopyType type = ::NOX::DeepCopy) const override;

      //! compute/update the current state variables
      void computeX(const NOX::Nln::Group& grp, const ::NOX::Epetra::Vector& d, double step);
      void computeX(const ::NOX::Abstract::Group& grp, const ::NOX::Abstract::Vector& d,
          double step) override;

      ::NOX::Abstract::Group::ReturnType computeF() override;

      ::NOX::Abstract::Group::ReturnType applyJacobianInverse(Teuchos::ParameterList& p,
          const ::NOX::Epetra::Vector& input, ::NOX::Epetra::Vector& result) const override;

      /// apply one block of the jacobian to a vector
      ::NOX::Abstract::Group::ReturnType apply_jacobian_block(const ::NOX::Epetra::Vector& input,
          Teuchos::RCP<::NOX::Epetra::Vector>& result, unsigned rbid, unsigned cbid) const;

      //! Compute and store \f$F(x)\f$ and the jacobian \f$\frac{\partial F(x)}{\partial x}\f$ at
      //! the same time. This can result in a huge performance gain in some special cases, e.g.
      //! contact problems.
      virtual ::NOX::Abstract::Group::ReturnType compute_f_and_jacobian();

      //! ToDo Move this into an extra interface
      //! @{

      /// compute element volumes
      ::NOX::Abstract::Group::ReturnType compute_element_volumes(
          Teuchos::RCP<Epetra_Vector>& ele_vols) const;

      /*! get the nodal dofs from the elements corresponding to the provided
       *  global element ids */
      void get_dofs_from_elements(
          const std::vector<int>& my_ele_gids, std::set<int>& my_ele_dofs) const;

      /// compute trial element volumes
      ::NOX::Abstract::Group::ReturnType compute_trial_element_volumes(
          Teuchos::RCP<Epetra_Vector>& ele_vols, const ::NOX::Abstract::Vector& dir, double step);

      //! @}

      /// compute the correction system of equations (e.g. in case of a second order correction
      /// step)
      ::NOX::Abstract::Group::ReturnType compute_correction_system(
          const enum NOX::Nln::CorrectionType type);

      //! set right hand side
      ::NOX::Abstract::Group::ReturnType set_f(Teuchos::RCP<::NOX::Epetra::Vector> Fptr);

      //! set jacobian operator
      ::NOX::Abstract::Group::ReturnType set_jacobian_operator(
          const Teuchos::RCP<const Epetra_Operator> jacOperator);

      //! set the solution vector to zero
      void reset_x();

      //! set flag whether update of x vector should be skipped (because it has already be done in
      //! preComputeX)
      void set_skip_update_x(bool skipUpdateX);

      /* Check the isValidJacobian flag and the ownership of the linear system
       * separately and get the ownership, if necessary. This prevents unnecessary
       * evaluation calls of the expensive computeJacobian() routines! Afterwards
       * the base class function is called.                   hiermeier 03/2016 */
      bool isJacobian() const override;

      /// are the eigenvalues valid?
      virtual bool is_eigenvalues() const;

      inline NOX::Nln::CorrectionType get_correction_type() const { return corr_type_; }

      /// @name access the eigenvalue data
      /// @{

      const Core::LinAlg::SerialDenseVector& get_jacobian_real_eigenvalues() const;
      const Core::LinAlg::SerialDenseVector& get_jacobian_imaginary_eigenvalues() const;
      double get_jacobian_max_real_eigenvalue() const;
      double get_jacobian_min_real_eigenvalue() const;

      /// @}

      //! returns the nox_nln_interface_required pointer
      Teuchos::RCP<const NOX::Nln::Interface::Required> get_nln_req_interface_ptr() const;

      //! returns the primary rhs norms
      virtual Teuchos::RCP<const std::vector<double>> get_rhs_norms(
          const std::vector<::NOX::Abstract::Vector::NormType>& type,
          const std::vector<NOX::Nln::StatusTest::QuantityType>& chQ,
          Teuchos::RCP<const std::vector<::NOX::StatusTest::NormF::ScaleType>> scale =
              Teuchos::null) const;

      //! returns the Root Mean Squares (abbr.: RMS) of the primary solution updates
      virtual Teuchos::RCP<std::vector<double>> get_solution_update_rms(
          const ::NOX::Abstract::Vector& xOld, const std::vector<double>& aTol,
          const std::vector<double>& rTol,
          const std::vector<NOX::Nln::StatusTest::QuantityType>& chQ,
          const std::vector<bool>& disable_implicit_weighting) const;

      double get_trial_update_norm(const ::NOX::Abstract::Vector& dir,
          const ::NOX::Abstract::Vector::NormType normtype, const StatusTest::QuantityType quantity,
          const StatusTest::NormUpdate::ScaleType scale = StatusTest::NormUpdate::Unscaled) const;

      //! returns the desired norm of the primary solution updates
      virtual Teuchos::RCP<std::vector<double>> get_solution_update_norms(
          const ::NOX::Abstract::Vector& xOld,
          const std::vector<::NOX::Abstract::Vector::NormType>& type,
          const std::vector<StatusTest::QuantityType>& chQ,
          Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale =
              Teuchos::null) const;

      //! returns the desired norm of the previous solution
      virtual Teuchos::RCP<std::vector<double>> get_previous_solution_norms(
          const ::NOX::Abstract::Vector& xOld,
          const std::vector<::NOX::Abstract::Vector::NormType>& type,
          const std::vector<StatusTest::QuantityType>& chQ,
          Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale) const;

      //! create a backup state
      void create_backup_state(const ::NOX::Abstract::Vector& dir) const;

      //! recover from stored backup state
      void recover_from_backup_state();

      //! @name reset the pre/post operator wrapper objects
      //! @{
      /*! \brief Resets the pre/post operator for the nln group
       *  Default call to the two parameter version, without resetting the isValid flags.
       *  @param[in] grpOptionParams   ParameterList which holds the new pre/post operator. */
      void reset_pre_post_operator(Teuchos::ParameterList& grpOptionParams)
      {
        reset_pre_post_operator(grpOptionParams, false);
      };

      /*! \brief Resets the pre/post operator wrapper for the nln group
       *  @param[in] grpOptionsParams   ParameterList which holds the new pre/post operator
       *  @param[in] resetIsValidFlag   If true, this forces the computeJacobian(), computeF() etc.
       * routines to reevaluate the linear system after setting a new pre/post operator. */
      void reset_pre_post_operator(
          Teuchos::ParameterList& grpOptionParams, const bool& resetIsValidFlags);

      /*! \brief Resets the pre/post operator wrapper for the nln linear system
       *  Default call to the two parameter version, without resetting the isValid flags.
       *  @param[in] linearSolverParams ParameterList which holds the new pre/post operator. */
      void reset_lin_sys_pre_post_operator(Teuchos::ParameterList& linearSolverParams)
      {
        reset_lin_sys_pre_post_operator(linearSolverParams, false);
      };

      /*! \brief Resets the pre/post operator wrapper for the nln linear system
       *  @param[in] linearSolverParams   ParameterList which holds the new pre/post operator.
       *  @param[in] resetIsValidFlag     If true, this forces the computeJacobian(), computeF()
       * etc. routines to reevaluate the linear system after setting a new pre/post operator. */
      void reset_lin_sys_pre_post_operator(
          Teuchos::ParameterList& linearSolverParams, const bool& resetIsValidFlags);
      //! @}

      //! @name PTC related methods
      //! @{
      //! adjust the pseudo time step length for the ptc nln solver
      void adjust_pseudo_time_step(double& delta, const double& stepSize,
          const ::NOX::Abstract::Vector& dir, const NOX::Nln::Solver::PseudoTransient& ptcsolver);
      void adjust_pseudo_time_step(double& delta, const double& stepSize,
          const ::NOX::Epetra::Vector& dir, const NOX::Nln::Solver::PseudoTransient& ptcsolver);

      //! get the lumped mass matrix
      Teuchos::RCP<const Epetra_Vector> get_lumped_mass_matrix_ptr() const;

      // Get element based scaling operator
      Teuchos::RCP<Core::LinAlg::SparseMatrix> get_contributions_from_element_level();
      //! @}

      //! @name XFEM related methods
      //! @{

      //! destroy the jacobian ptr in the linear system
      bool destroy_jacobian();

      //! @}

      //! compute and return some energy representative
      virtual double get_model_value(const enum MeritFunction::MeritFctName merit_func_type) const;

      //! compute contributions to a linear model
      virtual double get_linearized_model_terms(const ::NOX::Abstract::Vector& dir,
          enum NOX::Nln::MeritFunction::MeritFctName mf_type,
          enum NOX::Nln::MeritFunction::LinOrder linorder,
          enum NOX::Nln::MeritFunction::LinType lintype) const;

      /// return the DOF map of the solution vector
      const Epetra_BlockMap& get_dof_map() const;

      /// get jacobian range map from row/column block rbid/cbid
      const Epetra_Map& get_jacobian_range_map(unsigned rbid, unsigned cbid) const;

      /// return a copy of the Jacobian diagonal block \c diag_bid
      Teuchos::RCP<Epetra_Vector> get_diagonal_of_jacobian(unsigned diag_bid) const;

      /// replace the Jacobian diagonal block \c diag_bid
      void replace_diagonal_of_jacobian(const Epetra_Vector& new_diag, unsigned diag_bid) const;

      /// compute the condition number of the Jacobian matrix (serial mode)
      ::NOX::Abstract::Group::ReturnType compute_serial_jacobian_condition_number(
          const NOX::Nln::LinSystem::ConditionNumber condnum_type, bool printOutput);

      /// compute the eigenvalues of the Jacobian matrix (serial mode)
      ::NOX::Abstract::Group::ReturnType compute_serial_jacobian_eigenvalues(bool printOutput);

      /// reset is valid Newton member
      inline void reset_is_valid_newton() { isValidNewton = false; };

      /// allow to set isValidNewton manually
      inline void set_is_valid_newton(const bool value) { isValidNewton = value; };

      /// allow to set isValidRHS manually
      inline void set_is_valid_rhs(const bool value) { isValidRHS = value; };

     protected:
      //! resets the isValid flags to false
      void resetIsValid() override;

     private:
      //! Throw an NOX_error
      void throw_error(const std::string& functionName, const std::string& errorMsg) const;

     protected:
      /*! flag whether update of x vector should be skipped
       *  (e.g. because it has already be done in preComputeX as might be the case if we
       *  need a multiplicative update of some beam elements' rotation (pseudo-)vector DOFs) */
      bool skipUpdateX_;

      /// correction system type
      NOX::Nln::CorrectionType corr_type_;

      //! pointer to an user defined wrapped NOX::Nln::Abstract::PrePostOperator object.
      Teuchos::RCP<NOX::Nln::GROUP::PrePostOperator> prePostOperatorPtr_;

     private:
      /// container for eigenvalue info
      struct Eigenvalues
      {
        /// assign operator
        Eigenvalues& operator=(const Eigenvalues& src);

        /// real part of the eigenvalues
        Core::LinAlg::SerialDenseVector realpart_;

        /// imaginary part of the eigenvalues
        Core::LinAlg::SerialDenseVector imaginarypart_;

        /// maximal real part
        double real_max_ = 0.0;

        /// minimal real part
        double real_min_ = 0.0;

        /// Are the eigenvalues valid?
        bool isvalid_ = false;
      };

      /// instance of the Eigenvalue container
      Eigenvalues ev_;
    };  // class Group
  }     // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
