/*!----------------------------------------------------------------------
\file contact_augmented_strategy.H

<pre>
Created on: Apr 7, 2014

Maintainer: Michael Hiermeier
            hiermeier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15268
</pre>

*----------------------------------------------------------------------*/

#ifndef CONTACT_AUGMENTED_STRATEGY_H_
#define CONTACT_AUGMENTED_STRATEGY_H_

#include "../drt_contact/contact_lagrange_strategy.H"

namespace CONTACT
{
// forward declarartion
class AugmentedInterface;

class AugmentedLagrangeStrategy : public CoLagrangeStrategy
{
  public:
  /*!
  \brief Standard Constructor

  */
  AugmentedLagrangeStrategy(
      const Epetra_Map* DofRowMap,
      const Epetra_Map* NodeRowMap,
      Teuchos::ParameterList params,
      std::vector<Teuchos::RCP<CONTACT::CoInterface> > interfaces,
      int dim,
      Teuchos::RCP<Epetra_Comm> comm,
      double alphaf,
      int maxdof);

  /*!
  \brief Destructor

  */
  virtual ~AugmentedLagrangeStrategy() {};

  /*!
  \brief Return number of augmented active nodes

  */
  virtual int NumberOfActiveNodes()
  {
    if (gAugActiveSlaveNodes_!=Teuchos::null)
      return gAugActiveSlaveNodes_->NumGlobalElements();
    return 0;
  }

  /*!
  \brief Return global active node row map

  */
  virtual const Teuchos::RCP<Epetra_Map> ActiveRowNodes() { return gAugActiveSlaveNodes_; }

  /*!
  \brief Return global active dof row map

  */
  virtual const Teuchos::RCP<Epetra_Map> ActiveRowDofs() { return gAugActiveSlaveDofs_; }

  /*!
  \brief Return augmented constraint rhs vector

  */
  virtual Teuchos::RCP<Epetra_Vector> ConstrRhs() { return augConstrRhs_; }

  /*!
  \brief Do the active set Update. Build the local and global active set
  Epetra_maps and check for convergence of the active set.

  it : Number of the current interface

  */
  void BuildGlobalAugActiveSet(const int it);

  /*!
  \brief Initialize all matrices needed for the linearization of the
  augmented Lagrange formulation.

  */
  void Initialize();

  /*!
  \brief Evaluate all frictionless contact contributions for
  the augmented Lagrange formulation

  */
  void EvaluateContact(Teuchos::RCP<LINALG::SparseOperator>& kteff,
      Teuchos::RCP<Epetra_Vector>& feff);

  /*!
  \brief Outsourcing of the structural right hand side update.

  */
  void UpdateStructuralRHS(Teuchos::RCP<Epetra_Vector>& feff);

  /*!
  \brief Outsourcing of the structural stiffness matrix update.

  */
  void UpdateStructuralStiff(Teuchos::RCP<LINALG::SparseOperator>& kteff);

  /*!
  \brief The entries of the constraint right-hand side are
  calculated. This is implemented in a extra function cause
  the vector is needed for the line search algorithm for example.
  */
  void EvalConstrRHS();

  /*!
  \brief Evaluate the augmented contact forces on slave and master side

  */
  void AugForces(Epetra_Vector& augfs_lm, Epetra_Vector& augfs_g,
                 Epetra_Vector& augfm_lm, Epetra_Vector& augfm_g);

  /*!
  \brief Projection of the nodal LM values in the nodal normal and
  tangential direction. Pay attention, that the Lagrange multipliers
  don't have any directions, they are scalar values. This function is just
  used to visualize the LM values.
  A better and more realistic visualization is the corresponding contact force.

  Function overwrites the basic function in the abstract lagrange strategy.

  */
  void OutputStresses();

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
  void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd, Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<LINALG::MapExtractor> dbcmaps, int numiter,Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs);

   /*!
    \brief Update internal member variables after solving the 2x2 saddle point contact system

    \param sold (in): the displacement dof solution increment (associated with displacement dofs)
    \param mergedsol (in): Epetra_Vector for merged solution vector (containing the new solution vector of the full merged linear system)
  */
  void UpdateDisplacementsAndLMincrements(Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Epetra_Vector> blocksol);

  //! @name Debugging and visualization methods

  /*!
  \brief Check linear and angular momentum conservation

  */
  void CheckConservationLaws(const Epetra_Vector& augfs, const Epetra_Vector& augfm);

  /*!
  \brief Finite Difference check of the augmented Lagrange terms at GP level

  */
  void AugFDCheckGP();

  /*!
  \brief Finite Difference check of the augmented Lagrange terms at global level

  */
  void AugFDCheckGlobal();

  /*!
  \brief Evaluate the global FD check w.r.t. the displacements
  This function is only suitable for the augmented Lagrange formulation due to
  the special structure of the EvaluateContact implementation!

  */
  void EvalFDCheckGlobalDispl(LINALG::SparseMatrix derivMatrix,
                          Teuchos::RCP<Epetra_Vector>& rhsVector);
  //@}

  protected:
  // Don't want = operator and cctor
  AugmentedLagrangeStrategy operator = (const AugmentedLagrangeStrategy& old);
  AugmentedLagrangeStrategy(const AugmentedLagrangeStrategy& old);

  // booleans
  bool gFdCheck_;                                       // indicates if the global finite difference check is active

  // >>> interfaces <<<
  std::vector<Teuchos::RCP<CONTACT::AugmentedInterface> > interface_;

  // >>> matrices <<<
  Teuchos::RCP<LINALG::SparseMatrix> dGLmSlLinMatrix_;  // global matrix dGLmSlLin
  Teuchos::RCP<LINALG::SparseMatrix> dGLmMaLinMatrix_;  // global matrix dGLmMaLin
  Teuchos::RCP<LINALG::SparseMatrix> dGGSlLinMatrix_;   // global matrix dGGSlLin
  Teuchos::RCP<LINALG::SparseMatrix> dGGMaLinMatrix_;   // global matrix dGGMaLin

  Teuchos::RCP<LINALG::SparseMatrix> augDnMatrix_;
  Teuchos::RCP<LINALG::SparseMatrix> augMnMatrix_;

  Teuchos::RCP<LINALG::SparseMatrix> augDnLinMatrix_;
  Teuchos::RCP<LINALG::SparseMatrix> augMnLinMatrix_;

  Teuchos::RCP<LINALG::SparseMatrix> aWGapLinMatrix_;
  Teuchos::RCP<LINALG::SparseMatrix> dLmNWGapLinMatrix_;
  Teuchos::RCP<LINALG::SparseMatrix> dLmTLmTMatrix_;
  Teuchos::RCP<LINALG::SparseMatrix> dLmTLmTLinMatrix_;

  Teuchos::RCP<LINALG::SparseMatrix> augInactiveMatrix_;
  Teuchos::RCP<LINALG::SparseMatrix> augInactiveLinMatrix_;

  // >>> vectors <<<
  Teuchos::RCP<Epetra_Vector> augLm_;
  Teuchos::RCP<Epetra_Vector> aWGapRhs_;
  Teuchos::RCP<Epetra_Vector> dLmNWGapRhs_;
  Teuchos::RCP<Epetra_Vector> dLmNGapRhs_;
  Teuchos::RCP<Epetra_Vector> dLmNLmNRhs_;
  Teuchos::RCP<Epetra_Vector> dLmTLmTRhs_;
  Teuchos::RCP<Epetra_Vector> augInactiveRhs_;
  Teuchos::RCP<Epetra_Vector> augConstrRhs_;

  // vectors for debugging
  // slave side
  Teuchos::RCP<Epetra_Vector> augfs_lm_;
  Teuchos::RCP<Epetra_Vector> augfs_g_;
  // master side
  Teuchos::RCP<Epetra_Vector> augfm_lm_;
  Teuchos::RCP<Epetra_Vector> augfm_g_;

  // global maps
  Teuchos::RCP<Epetra_Map> gAugActiveSlaveNodes_;  // global map of all active slave nodes
  Teuchos::RCP<Epetra_Map> gAugActiveSlaveDofs_;   // global map of all active slave dofs
  Teuchos::RCP<Epetra_Map> gAugActiveSlaveNDofs_;  // global map of all active slave dofs in normal direction
  Teuchos::RCP<Epetra_Map> gAugActiveSlaveTDofs_;  // global map of all active slave dofs in tangential direction

};
}


#endif /* #ifndef CONTACT_AUGMENTED_STRATEGY_H_ */
