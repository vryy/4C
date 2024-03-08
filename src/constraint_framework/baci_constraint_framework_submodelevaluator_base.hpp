/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all constraint submodel evaluators.


\level 3
*/
#ifndef BACI_CONSTRAINT_FRAMEWORK_SUBMODELEVALUATOR_BASE_HPP
#define BACI_CONSTRAINT_FRAMEWORK_SUBMODELEVALUATOR_BASE_HPP

#include "baci_config.hpp"

#include "baci_constraint_framework_equation_mpc.hpp"
#include "baci_structure_new_model_evaluator_generic.hpp"

BACI_NAMESPACE_OPEN

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
    bool EvaluateForceStiff(Teuchos::RCP<CORE::LINALG::SparseMatrix> me_stiff_ptr,
        Teuchos::RCP<Epetra_Vector> me_force_ptr);

    //! Evaluate the matrices of the saddle-point system
    virtual void EvaluateCouplingTerms(STR::TIMINT::BaseDataGlobalState& gstate);

    //! Reset the sub model evaluator
    virtual void Reset() = 0;

    //! Return the Penalty-Parameter
    double& GetPenaltyParameterPtr() { return penaltyParameter_; }

   private:
    //! Colum Map
    Teuchos::RCP<Epetra_Map> n_condition_map_;

    //! Penalty parameter
    double penaltyParameter_;

   protected:
    //! Vector containing all multipoint constraint and related constraint equation objects
    std::vector<Teuchos::RCP<MultiPointConstraintEquationBase>> listMPCs_;

    //! Pointer to the structural stiffness matrix \f$ K_{sys} \f$
    CORE::LINALG::SparseMatrix* stiff_ptr_;

    //! Enforcement Strategy
    enum INPAR::RVE_MPC::enforcementStrategy strategy_ =
        INPAR::RVE_MPC::enforcementStrategy::penalty;

    //! Pointer to the discretization
    Teuchos::RCP<const DRT::Discretization> discret_ptr_ = Teuchos::null;

    //! System Matrix \f$ Q_{dd} \f$
    Teuchos::RCP<CORE::LINALG::SparseMatrix> Q_dd_ = Teuchos::null;

    //! System Matrix \f$ Q_{dL} \f$
    Teuchos::RCP<CORE::LINALG::SparseMatrix> Q_dL_ = Teuchos::null;

    //! System Matrix \f$ Q_{Ld} \f$
    Teuchos::RCP<CORE::LINALG::SparseMatrix> Q_Ld_ = Teuchos::null;

    //! coupling conditions evaluated at current displacements
    Teuchos::RCP<Epetra_Vector> constraint_vector_;
  };

}  // namespace CONSTRAINTS::SUBMODELEVALUATOR

BACI_NAMESPACE_CLOSE
#endif
