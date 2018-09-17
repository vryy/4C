/*---------------------------------------------------------------------*/
/*!
\file meshtying_lagrange_strategy.H

\brief Mesh-tying solving strategy with (standard/dual) Lagrangian
       multipliers.

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

   \author popp (popp@lnm.mw.tum.de)
   */
  class MtLagrangeStrategy : public MtAbstractStrategy
  {
   public:
    /*!
    \brief Standard Constructor

    */
    MtLagrangeStrategy(const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<MORTAR::MortarInterface>> interface,
        int dim, Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~MtLagrangeStrategy(){};

    /*!
    \brief Do mortar coupling in reference configuration

    Only do this ONCE for meshtying upon initialization!

    */
    void MortarCoupling(const Teuchos::RCP<Epetra_Vector> dis);

    /*!
    \brief Mesh intialization for rotational invariance

    Only do this ONCE for meshtying upon initialization!
    (only implemented for DUAL shape functions)

    */
    Teuchos::RCP<Epetra_Vector> MeshInitialization();

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
    \param numiter (in): current nonlinear iteration step
    \param blockMat (out): Epetra_Operator containing the 2x2 block sparse matrix object
    \param mergedsol (out): Epetra_Vector for merged solution vector
    \param mergedrhs (out): Epetra_Vector for merged right hand side vector
    */
    void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps, int numiter,
        Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol,
        Teuchos::RCP<Epetra_Vector>& blockrhs);

    /*!
    \brief Update internal member variables after solving the 2x2 saddle point contact system

    \param sold (in): the displacement dof solution increment (associated with displacement dofs)
    \param mergedsol (in): Epetra_Vector for merged solution vector (containing the new solution
    vector of the full merged linear system)
    */
    void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Epetra_Vector> blocksol);


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
    void Recover(Teuchos::RCP<Epetra_Vector> disi);

    //@}

    //! @name Empty functions (Penalty meshtying)

    // All these functions only have functionality in Penalty meshtying simulations,
    // thus they are defined as dserror here in the case of Lagrange meshtying..

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

    //@}

   protected:
    // don't want = operator and cctor
    MtLagrangeStrategy operator=(const MtLagrangeStrategy& old);
    MtLagrangeStrategy(const MtLagrangeStrategy& old);

    Teuchos::RCP<LINALG::SparseMatrix> conmatrix_;   // constraint matrix for saddle point system
    Teuchos::RCP<LINALG::SparseMatrix> mhatmatrix_;  // product of global Mortar matrices inv(D)*M

    Teuchos::RCP<Epetra_Vector> fs_;           // slave side effective forces (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> invd_;  // inverse of Mortar matrix D (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> ksn_;   // stiffness block K_sn (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> ksm_;   // stiffness block K_sm (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> kss_;   // stiffness block K_ss (needed for LM)

  };  // class MtLagrangeStrategy
}  // namespace CONTACT


#endif  // #ifndef MESHTYING_LAGRANGE_STRATEGY_H
