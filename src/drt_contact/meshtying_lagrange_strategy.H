/*---------------------------------------------------------------------*/
/*!
\file meshtying_lagrange_strategy.H

\brief Mesh-tying solving strategy with (standard/dual) Lagrangian multipliers

\level 2

\maintainer Matthias Mayr

*/
/*---------------------------------------------------------------------*/
#ifndef MESHTYING_LAGRANGE_STRATEGY_H
#define MESHTYING_LAGRANGE_STRATEGY_H

#include "meshtying_abstract_strategy.H"


namespace CONTACT
{
  /*!
   \brief Meshtying solving strategy with (standard/dual) Lagrangian multipliers.

   This is a specialization of the abstract meshtying algorithm as defined in MtAbstractStrategy.
   For a more general documentation of the involved functions refer to MtAbstractStrategy.

   \sa MtAbstractStrategy

   \author popp (popp@lnm.mw.tum.de)
   */
  class MtLagrangeStrategy : public MtAbstractStrategy
  {
   public:
    /*!
    \brief Standard Constructor

    \param[in] DofRowMap Dof row map of underlying problem
    \param[in] NodeRowMap Node row map of underlying problem
    \param[in] params List of contact/parameters
    \param[in] interface All contact interface objects
    \param[in] spatialDim Spatial dimension of the problem
    \param[in] comm Communicator
    \param[in] alphaf Mid-point for Generalized-alpha time integration
    \param[in] maxdof Highest DOF number in global problem
    */
    MtLagrangeStrategy(const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<MORTAR::MortarInterface>> interface,
        const int spatialDim, const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf,
        const int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~MtLagrangeStrategy() = default;

    /*!
    \brief Do mortar coupling in reference configuration

    Only do this ONCE for meshtying upon initialization!

    */
    void MortarCoupling(const Teuchos::RCP<const Epetra_Vector>& dis) override;

    /*!
    \brief Mesh initialization for rotational invariance

    Compute necessary modifications to the reference configuration of interface nodes, such that the
    weighted gap in the modified reference configuration is zero.

    \note Only do this \em once for meshtying upon initialization!

    \warning This is only implemented for mortar coupling. No implementation for node-to-segment
    approach.

    \return Vector with modified nodal positions
    */
    Teuchos::RCP<const Epetra_Vector> MeshInitialization() override;

    /*!
    \brief Evaluate meshtying

    This is the main routine of our meshtying algorithms on a global level.
    It contains the setup of the global linear system including meshtying.

    For a Lagrangian strategy this involves heavy modification to the initial kteff and feff.
    Hence, they are in fact build from scratch here. The application of modifications to
    groups of dofs (slave, master, etc.) results in some matrix and vector splitting and a
    lot of matrix-vector calculation in here!

    \param kteff (in/out): effective stiffness matrix (without -> with contact)
    \param feff (in/out): effective residual / force vector (without -> with contact)
    \param dis (in): current displacement state

    */
    void EvaluateMeshtying(Teuchos::RCP<LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis);

    /*!
    \brief Build 2x2 saddle point system

    \param kdd (in): the displacement dof stiffness (upper left block)
    \param fd (in): the displacement dof r.h.s. (upper block)
    \param sold (in): the displacement dof solution increment
    \param dirichtoggle (in): toggle vector for dirichlet conditions
    \param blockMat (out): Epetra_Operator containing the 2x2 block sparse matrix object
    \param mergedsol (out): Epetra_Vector for merged solution vector
    \param mergedrhs (out): Epetra_Vector for merged right hand side vector
    */
    void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) override;

    /*!
    \brief Update internal member variables after solving the 2x2 saddle point contact system

    \param sold (out): the displacement dof solution increment (associated with displacement dofs)
    \param mergedsol (in): Epetra_Vector for merged solution vector (containing the new solution
    vector of the full merged linear system)
    */
    void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) override;


    void EvalConstrRHS()
    {
      std::cout << "Warning: The EvalConstrRHS() function is not yet implemented for meshtying."
                << std::endl;
    }


    /*!
    \brief Recovery method

    We only recover the Lagrange multipliers here, which had been statically condensed during
    the setup of the global problem!

    */
    void Recover(Teuchos::RCP<Epetra_Vector> disi) override;

    //@}

    /*! @name Empty functions (Penalty meshtying)
     *
     * All these functions only have functionality in Penalty meshtying simulations,
     * thus they are defined as dserror here in the case of Lagrange meshtying.
     */

    double ConstraintNorm() const { return 0.0; }
    void InitializeUzawa(
        Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
    {
    }
    double InitialPenalty() { return 0.0; }
    void ResetPenalty() {}
    void ModifyPenalty() {}
    void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) {}
    void UpdateUzawaAugmentedLagrange() {}
    void UpdateConstraintNorm(int uzawaiter = 0) {}
    bool IsPenalty() const override { return false; };

    //@}

    //! @name New time integration
    //!@{

    /*! \brief Evaluate residual
     *
     * @param[in] dis Current displacement field
     * @return Boolean flag indicating successfull evaluation
     */
    bool EvaluateForce(const Teuchos::RCP<const Epetra_Vector> dis) override;

    /*! \brief Evaluate stiffness term
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

    /*! \brief Modify system before linear solve
     *
     * Perform the static condensation of mortar terms if a condensed formulation is used.
     *
     * This exploit the biorthogonality condition of the dual shape functions.
     *
     * @param[in/out] kteff Stiffness matrix
     * @param[in/out] rhs right-hand side vector
     *
     * \warning This only works for dual shape functions. Standard shape functions are prohibited
     * as they are too expensive.
     *
     * \todo Is this really the right-hand side vector or the residual?
     */
    virtual void RunPreApplyJacobianInverse(
        Teuchos::RCP<LINALG::SparseMatrix> kteff, Epetra_Vector& rhs);

    virtual void RunPostApplyJacobianInverse(Epetra_Vector& result);

    void RunPostComputeX(
        const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

    void RemoveCondensedContributionsFromRhs(Epetra_Vector& rhs) const override;

    //!@}

   protected:
    //! @name Accessors
    //!@{

    //! Access to #mhatmatrix_
    virtual Teuchos::RCP<const LINALG::SparseMatrix> GetMHat() { return mhatmatrix_; };

    //! Access to #invd_
    virtual Teuchos::RCP<const LINALG::SparseMatrix> GetDInverse() { return invd_; };

    //!@}

   private:
    //! don't want = operator
    MtLagrangeStrategy operator=(const MtLagrangeStrategy& old);

    //! don't want copy constructor
    MtLagrangeStrategy(const MtLagrangeStrategy& old);

    //! Constraint matrix for saddle point system
    Teuchos::RCP<LINALG::SparseMatrix> conmatrix_;

    //! Mortar projection matrix \f$P = D^{-1} M\f$
    Teuchos::RCP<LINALG::SparseMatrix> mhatmatrix_;

    //! Slave side effective forces (needed for Lagrange multipliers)
    Teuchos::RCP<Epetra_Vector> fs_;

    //! Inverse \f$D^{-1}\f$ of Mortar matrix \f$D\f$ (needed for Lagrange multipliers)
    Teuchos::RCP<LINALG::SparseMatrix> invd_;

    /*! @name Blocks for Jacobian matrix
     *
     * Subscripts are defined as follows
     * - s: slave
     * - m: master
     * - n: non-mortar nodes (i.e. all non-interface nodes in all subdomains)
     */
    //!@{

    //! Stiffness block \f$K_{sn}\f$ (needed for Lagrange multipliers)
    Teuchos::RCP<LINALG::SparseMatrix> ksn_;

    //! Stiffness block \f$K_{sm}\f$ (needed for Lagrange multipliers)
    Teuchos::RCP<LINALG::SparseMatrix> ksm_;

    //! Stiffness block \f$K_{ss}\f$ (needed for Lagrange multipliers)
    Teuchos::RCP<LINALG::SparseMatrix> kss_;

    //!@}

  };  // class MtLagrangeStrategy
}  // namespace CONTACT


#endif  // #ifndef MESHTYING_LAGRANGE_STRATEGY_H
