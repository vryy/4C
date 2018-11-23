/*---------------------------------------------------------------------*/
/*!
\file meshtying_penalty_strategy.H

\brief Penalty mesh-tying solving strategy.

\level 2

\maintainer Matthias Mayr

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
        int dim, const Teuchos::RCP<const Epetra_Comm>& comm, double alphaf, int maxdof);

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
    void MortarCoupling(const Teuchos::RCP<const Epetra_Vector>& dis);

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

    /*!
    \brief Tell that this is a penalty strategy
    */
    bool IsPenalty() const override { return true; };
    //@}

    //! @name Empty functions (Lagrange meshtying)

    // All these functions only have functionality in Lagrange meshtying simulations,
    // thus they are defined empty here in the case of Penalty meshtying.

    void Recover(Teuchos::RCP<Epetra_Vector> disi) {}
    void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) override
    {
      dserror(
          "A penalty approach does not have Lagrange multiplier DOFs. So, saddle point system "
          "makes no sense here.");
    };
    void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) override
    {
      dserror(
          "A penalty approach does not have Lagrange multiplier DOFs. So, saddle point system "
          "makes no sense here.");
    };
    void EvalConstrRHS()
    {
      std::cout << "Warning: No constraint RHS in contact penalty strategy" << std::endl;
    }
    //@}

    //! @name New time integration

    /*! \brief Evaluate residual
     *
     * @param[in] dis Current displacement field
     * @return Boolean flag indicating successfull evaluation
     */
    bool EvaluateForce(const Teuchos::RCP<const Epetra_Vector> dis) override;

    /*! \brief Evaluate stiffness term
     *
     * \note We assume that the stiffness matrix has already been evaluated,
     * so we do nothing in here.
     *
     * @param[in] dis Current displacement field
     * @return Boolean flag indicating successfull evaluation
     */
    bool EvaluateStiff(const Teuchos::RCP<const Epetra_Vector> dis) override;

    /*! \brief Evaluate residual and stiffness matrix
     *
     * @param[in] dis Current displacement field
     * @return Boolean flag indicating successfull evaluation
     */
    bool EvaluateForceStiff(const Teuchos::RCP<const Epetra_Vector> dis) override;

    //! Return the desired right-hand-side block pointer (read-only) [derived]
    Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum DRT::UTILS::VecBlockType& bt) const override;

    //! Return the desired matrix block pointer (read-only) [derived]
    Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt) const override;

    //@}

   protected:
    // don't want = operator and cctor
    MtPenaltyStrategy operator=(const MtPenaltyStrategy& old);
    MtPenaltyStrategy(const MtPenaltyStrategy& old);

    //! L2-norm of normal contact constraints
    double constrnorm_;

    //! Initial penalty parameter
    double initialpenalty_;

    //! Mortar matrix product \f$M^T M\f$
    Teuchos::RCP<LINALG::SparseMatrix> mtm_;

    //! Mortar matrix product \f$M^T D\f$
    Teuchos::RCP<LINALG::SparseMatrix> mtd_;

    //! Mortar matrix product \f$D^T M\f$
    Teuchos::RCP<LINALG::SparseMatrix> dtm_;

    //! Mortar matrix product \f$D^T D\f$
    Teuchos::RCP<LINALG::SparseMatrix> dtd_;


    //! Global stiffness matrix
    Teuchos::RCP<LINALG::SparseMatrix> stiff_;

    /*! \brief Global residual vector
     *
     * \todo Is this the residual or the right-hand side vector?
     */
    Teuchos::RCP<Epetra_Vector> force_;


  };  // class MtPenaltyStrategy
}  // namespace CONTACT

#endif  // #ifndef MESHTYING_PENALTY_STRATEGY_H
