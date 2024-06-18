/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all constraint submodel evaluators.


\level 3
*/
#ifndef FOUR_C_CONSTRAINT_FRAMEWORK_SUBMODELEVALUATOR_BASE_HPP
#define FOUR_C_CONSTRAINT_FRAMEWORK_SUBMODELEVALUATOR_BASE_HPP

#include "4C_config.hpp"

#include "4C_constraint_framework_equation_mpc.hpp"
#include "4C_structure_new_model_evaluator_generic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONSTRAINTS::SUBMODELEVALUATOR
{
  /*! \brief Interface class of all submodel evaluators managing
   *  constraint terms
   *
   * This is the abstract base class of all submodel evaluators for constraint
   * problems
   *
   *  This class summarizes the functionality which all submodel evaluators share
   *  and/or have to implement. Look in the derived classes for examples.
   *
   * Definitions:
   * \f$ \begin{bmatrix} K_{sys} + Q_{dd} & Q_{dL} \\ Q_{Ld} & 0 \end{bmatrix} +
   * \begin{bmatrix}  \Delta D \\ \lambda \end{bmatrix} =
   * \begin{bmatrix}  r_d \\  r_L \end{bmatrix}\f$
   */
  class ConstraintBase
  {
   public:
    //! constructor
    ConstraintBase() = default;

    //! destructor
    virtual ~ConstraintBase() = default;

    /*! Evaluate the current right-hand-side vector and tangential stiffness matrix at \f$t_{n+1}\f$
     */
    bool evaluate_force_stiff(Teuchos::RCP<Core::LinAlg::SparseMatrix> me_stiff_ptr,
        Teuchos::RCP<Epetra_Vector> me_force_ptr);

    //! Evaluate the matrices of the saddle-point system
    virtual void evaluate_coupling_terms(STR::TimeInt::BaseDataGlobalState& gstate);

    //! Reset the sub model evaluator
    virtual void reset() = 0;

    //! Return the Penalty-Parameter
    double& get_penalty_parameter_ptr() { return penalty_parameter_; }

   private:
    //! Colum Map
    Teuchos::RCP<Epetra_Map> n_condition_map_;

    //! Penalty parameter
    double penalty_parameter_;

   protected:
    //! Vector containing all multipoint constraint and related constraint equation objects
    std::vector<Teuchos::RCP<MultiPointConstraintEquationBase>> listMPCs_;

    //! Pointer to the structural stiffness matrix \f$ K_{sys} \f$
    Core::LinAlg::SparseMatrix* stiff_ptr_;

    //! Enforcement Strategy
    enum Inpar::RveMpc::EnforcementStrategy strategy_ = Inpar::RveMpc::EnforcementStrategy::penalty;

    //! Pointer to the discretization
    Teuchos::RCP<const Core::FE::Discretization> discret_ptr_ = Teuchos::null;

    //! System Matrix \f$ Q_{dd} \f$
    Teuchos::RCP<Core::LinAlg::SparseMatrix> Q_dd_ = Teuchos::null;

    //! System Matrix \f$ Q_{dL} \f$
    Teuchos::RCP<Core::LinAlg::SparseMatrix> Q_dL_ = Teuchos::null;

    //! System Matrix \f$ Q_{Ld} \f$
    Teuchos::RCP<Core::LinAlg::SparseMatrix> Q_Ld_ = Teuchos::null;

    //! coupling conditions evaluated at current displacements
    Teuchos::RCP<Epetra_Vector> constraint_vector_;
  };

}  // namespace CONSTRAINTS::SUBMODELEVALUATOR

FOUR_C_NAMESPACE_CLOSE
#endif
