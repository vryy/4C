/*---------------------------------------------------------------------*/
/*!
\file meshtying_penalty_strategy.H

\brief Penalty mesh-tying solving strategy.

\level 2

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#ifndef MESHTYING_PENALTY_STRATEGY_H
#define MESHTYING_PENALTY_STRATEGY_H

#include "meshtying_abstract_strategy.H"


namespace CONTACT
{
  /*!
   \brief Meshtying solving strategy with regularization of Lagrangian multipliers,
   also known as Penalty Method or regularization. An Augmented Lagrangian version
   based on the Uzawa algorithm is included, too.

   This is a specialization of the abstract meshtying algorithm as defined in MtAbstractStrategy.
   For a more general documentation of the involved functions refer to MtAbstractStrategy.

   \author popp (popp@lnm.mw.tum.de)
   */
  class MtPenaltyStrategy : public MtAbstractStrategy
  {
   public:
    /*!
    \brief Standard Constructor

    */
    MtPenaltyStrategy(const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<MORTAR::MortarInterface>> interface,
        int dim, Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~MtPenaltyStrategy(){};

    //! @name Access methods

    /*!
    \brief Return L2-norm of active constraints

    */
    double ConstraintNorm() const { return constrnorm_; }

    /*!
    \brief Return initial penalty parameter

    */
    double InitialPenalty() { return initialpenalty_; }

    //@}

    //! @name Evaluation methods

    /*!
    \brief Do mortar coupling in reference configuration

    Only do this ONCE for meshtying upon initialization!

    */
    void MortarCoupling(const Teuchos::RCP<Epetra_Vector> dis);

    /*!
    \brief Mesh intialization for rotational invariance

    Only do this ONCE for meshtying upon initialization!

    */
    Teuchos::RCP<Epetra_Vector> MeshInitialization();

    /*!
    \brief Evaluate meshtying

    This is the main routine of our meshtying algorithms on a global level.
    It contains the setup of the global linear system including meshtying.

    For a penalty strategy this includes the evaluation of regularized forces
    and results in a simple addition of extra stiffness contributions to kteff
    and extra meshtying forces to feff.

    \param kteff (in/out): effective stiffness matrix (without -> with contact)
    \param feff (in/out): effective residual / force vector (without -> with contact)
    \param dis (in): current displacement state

    */
    void EvaluateMeshtying(Teuchos::RCP<LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis);

    /*!
    \brief Initialize Uzawa step

    LM is updated to z = zuzawa - pp * gap. This mehtod is called at the
    beginning of the second, third, ... Uzawa iterarion in order to
    create an out-of-balance force again.

    */
    void InitializeUzawa(
        Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff);

    /*!
    \brief Reset penalty parameter to intial value

    When applying an Augmented Lagrangian version of the penalty approach,
    the penalty parameter is sometimes updated during the Uzawa steps in
    order to accelerate convergence of the constraint norm. This increase
    in penalty stiffness can be dealt with, because at the time it is applied
    the constraint norm is already quite low. Yet, for a new time step, we have
    to come back to the initial penalty parameter. Thus, this method is called
    at the beginning of each time step and resets the penalty parameter to its initial value.

    */
    void ResetPenalty();

    void ModifyPenalty();

    /*!
    \brief Compute L2-norm of active constraints

    In a classical penalty approach, the constraint norm is only monitored.
    When applying an Augmented Lagrangian version, the constraint norm is the
    relevant stopping criterion of the Uzawa iteration. In order to accelerate
    convergence, a heuristic update formula for the penalty parameter is applied
    in this method, too.

    */
    void UpdateConstraintNorm(int uzawaiter = 0);

    /*!
    \brief Store Lagrange multipliers for next Uzawa step

    A method ONLY called for the Uzawa Augmented Lagrangian version of the penalty method.
    At the end of an Uzawa step, the converged Lagrange multiplier value is stored
    in the variable zuzawa_, which is then used in the next Uzawa step.

    */
    void UpdateUzawaAugmentedLagrange();

    //@}

    //! @name Empty functions (Lagrange meshtying)

    // All these functions only have functionality in Lagrange meshtying simulations,
    // thus they are defined empty here in the case of Penalty meshtying.

    void Recover(Teuchos::RCP<Epetra_Vector> disi) {}
    void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps, int numiter,
        Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol,
        Teuchos::RCP<Epetra_Vector>& blockrhs)
    {
      std::cout << "WARNING: no saddlepoint solver in penalty meshtying strategy" << std::endl;
    };
    void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Epetra_Vector> blocksol)
    {
      std::cout << "WARNING: no saddlepoint solver in penalty meshtying strategy" << std::endl;
    };
    void EvalConstrRHS()
    {
      std::cout << "Warning: No constraint RHS in contact penalty strategy" << std::endl;
    }
    //@}

   protected:
    // don't want = operator and cctor
    MtPenaltyStrategy operator=(const MtPenaltyStrategy& old);
    MtPenaltyStrategy(const MtPenaltyStrategy& old);

    double constrnorm_;      // L2-norm of normal contact constraints
    double initialpenalty_;  // initial penalty parameter

    Teuchos::RCP<LINALG::SparseMatrix> mtm_;  // mortar matrix product: M(^T)*M
    Teuchos::RCP<LINALG::SparseMatrix> mtd_;  // mortar matrix product: M(^T)*D
    Teuchos::RCP<LINALG::SparseMatrix> dtm_;  // mortar matrix product: D(^T)*M
    Teuchos::RCP<LINALG::SparseMatrix> dtd_;  // mortar matrix product: D(^T)*D


  };  // class MtPenaltyStrategy
}  // namespace CONTACT

#endif  // #ifndef MESHTYING_PENALTY_STRATEGY_H
