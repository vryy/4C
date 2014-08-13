/*!----------------------------------------------------------------------
\file contact_augmented_interface.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Created on: Apr 16, 2014

Maintainer: Michael Hiermeier
            hiermeier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15268
</pre>

*----------------------------------------------------------------------*/

#ifndef CONTACT_AUGMENTED_INTERFACE_H_
#define CONTACT_AUGMENTED_INTERFACE_H_

#include "../drt_contact/contact_interface.H"

namespace CONTACT
{
class AugmentedInterface : public CoInterface
{
public:

  /*!
  \brief Constructor

  */
  AugmentedInterface(const int id, const Epetra_Comm& comm, const int dim,
      const Teuchos::ParameterList& icontact,
      bool selfcontact, INPAR::MORTAR::RedundantStorage redundant);

  /*!
  \brief Destructor

  */
  virtual ~AugmentedInterface() {}

  //! @name Access methods

  /*!
  \brief Get row map of augmented active slave nodes
  */
  const Teuchos::RCP<Epetra_Map> AugActiveSlaveNodes() const
  {
    if (Filled()) return augActiveSlaveNodes_;
    else dserror("CONTACT::AugmentedInterface::FillComplete was not called");
    exit(EXIT_FAILURE);
  }

  /*!
  \brief Get row map of augmented active slave dofs
  */
  const Teuchos::RCP<Epetra_Map> AugActiveSlaveDofs() const
  {
    if (Filled()) return augActiveSlaveDofs_;
    else dserror("CONTACT::AugmentedInterface::FillComplete was not called");
    exit(EXIT_FAILURE);
  }

  /*!
  \brief Get row map of augmented active normal slave dofs
  */
  const Teuchos::RCP<Epetra_Map> AugActiveSlaveNDofs() const
  {
    if (Filled()) return augActiveSlaveNDofs_;
    else dserror("CONTACT::AugmentedInterface::FillComplete was not called");
    exit(EXIT_FAILURE);
  }

  /*!
  \brief Get row map of augmented active tangential slave dofs
  */
  const Teuchos::RCP<Epetra_Map> AugActiveSlaveTDofs() const
  {
    if (Filled()) return augActiveSlaveTDofs_;
    else dserror("CONTACT::AugmentedInterface::FillComplete was not called");
    exit(EXIT_FAILURE);
  }

  /*!
  \brief Get row map of augmented active master nodes
  */
  const Teuchos::RCP<Epetra_Map> AugActiveMasterNodes() const
  {
    if (Filled()) return augActiveMasterNodes_;
    else dserror("CONTACT::AugmentedInterface::FillComplete was not called");
    exit(EXIT_FAILURE);
  }

  /*!
  \brief Get row map of augmented active master dofs
  */
  const Teuchos::RCP<Epetra_Map> AugActiveMasterDofs() const
  {
    if (Filled()) return augActiveMasterDofs_;
    else dserror("CONTACT::AugmentedInterface::FillComplete was not called");
    exit(EXIT_FAILURE);
  }

  //@}

  //! @name Initialize and evaluate interface element contributions

  /*!
  \brief Initialization of all augmented contact related quantities

  */
  void Initialize();

  /*!
  \brief Reduced evaluate of the contact interface. We don't need any
  search algorithm and we don't have to build the nodal normals again.
  */

  void RedEvaluate();

  //@}

  //! @name Assemble partial global matrices

  /*----------------------------------------------------------------*
   |  Assemble right hand side vectors, Dn and Mn                   |
   *----------------------------------------------------------------*/

  /*!
  \brief Assemble Dn and Mn matrices (equivalent to D and M in the Lagrange formulation)

  */
  void AssembleAugDnMnMatrix(LINALG::SparseMatrix& augDnMatrix,
                             LINALG::SparseMatrix& augMnMatrix);

  /*!
  \brief Assemble the augmented Lagrange multiplier vector

  */
  void AssembleAugLmVector(Epetra_Vector& augLmVec);

  /*!
  \brief Assemble averaged weighted gap vector
  (necessary for the conservation check and the augmented force output)

  */
  void AssembleAWGapRhs(Epetra_Vector& aWGapRhs);

  /*!
  \brief Assemble the active normal constraint right hand side

  */
  void AssembleDLmNWGapRhs(Epetra_Vector& dLmNWGapRhs);

  /*!
  \brief Assemble tangential constraint equation for active nodes (dLmTLmTrhs, ACTIVE)

  */
  void AssembleDLmTLmTRhs(Epetra_Vector& dLmTLmTRhs);

  /*!
  \brief Assemble inactive constraint rhs
  Combination of tangential and normal parts
  */
  void AssembleAugInactiveRhs(Epetra_Vector& augInactiveRhs);

  /*----------------------------------------------------------------*
   |  Assemble linearization matrices                               |
   *----------------------------------------------------------------*/

  /*!
  \brief Assemble DGLmLinMatrix
  Linearization w.r.t. the displacements

  */
  void AssembleDGLmLinMatrix(LINALG::SparseMatrix& dGLmSlLinMatrix,
                             LINALG::SparseMatrix& dGLmMaLinMatrix);

  /*!
  \brief Assemble DGGLinMatrix
  Linearization w.r.t. the displacements

  */
  void AssembleDGGLinMatrix(LINALG::SparseMatrix& dGGSlLinMatrix,
                            LINALG::SparseMatrix& dGGMaLinMatrix);

  /*!
  \brief Assemble DLmNWGapLinMatrix
  Linearization w.r.t. the displ.

  */
  void AssembleDLmNWGapLinMatrix(LINALG::SparseMatrix& dLmNWGapLinMatrix);

  /*!
  \brief Assemble DLmTLmTMatrix
  Linearization w.r.t. the LM

  */
  void AssembleDLmTLmTMatrix(LINALG::SparseMatrix& dLmTLmTMatrix);

  /*!
  \brief Assemble DLmTLmTLinMatrix
  Linearization w.r.t. the displ.

  */
  void AssembleDLmTLmTLinMatrix(LINALG::SparseMatrix& dLmTLmTLinMatrix);

  /*!
  \brief Assemble AugInactiveMatrix
  Linearization w.r.t. the LM

  */
  void AssembleAugInactiveMatrix(LINALG::SparseMatrix& augInactiveMatrix);

  /*!
  \brief Assemble AugInactiveLinMatrix
  Linearization w.r.t. the displ.

  */
  void AssembleAugInactiveLinMatrix(LINALG::SparseMatrix& augInactiveLinMatrix);

  //@}

  //! @name Calculate nodal quantities
  /*!
  \brief Calculate the nodal weighted gap values

  */
  void WGap();

  /*!
  \brief Calculate the linearization of the averaged weighted gap

  */
  void AWGapLin();

  //@}

  //! @name Augmented active set decisions

  /*!
  \brief Augmented active set decision
  NOTE: Differs from the standard active set decision, because we make use
  of the averaged weighted gap and we do the active set decision earlier.

  */
  bool UpdateAugActiveSetSemiSmooth();

  /*!
  \brief Update augmented active set for each interface
  */
  void BuildAugActiveSet(bool init = false);

  /*!
  \brief Split augmented active dof set for each interface
  */
  void SplitAugActiveDofs();

  //@}

  //! @name Debugging methods of the augmented Lagrange formulation

  /*!
  \brief Check linearization of kappa with finite differences

  */
  void FDCheckKappaLin();

  /*!
  \brief Check linearization of augA with finite differences

  */
  void FDCheckAugALin();

  /*!
  \brief Check linearization of the averaged weighted gap with finite differences

  */
  void FDCheckAWGapLin();

  /*!
  \brief Check linearization of varWGapSl with finite differences

  */
  void FDCheckVarWGapLinSl();

  /*!
  \brief Check linearization of varWGapMa with finite differences

  */
  void FDCheckVarWGapLinMa();

  /*!
  \brief Update Interfaces for global finite difference-check
  */
  bool UpdateInterfaces(int gid,int dof,double delta,bool forward);
  //@}

protected:
  // don't want = operator and cctor
  AugmentedInterface operator = (const AugmentedInterface& old);
  AugmentedInterface(const AugmentedInterface& old);

  Teuchos::RCP<Epetra_Map>    augActiveSlaveNodes_;
  Teuchos::RCP<Epetra_Map>    augActiveSlaveDofs_;
  Teuchos::RCP<Epetra_Map>    augActiveSlaveNDofs_;
  Teuchos::RCP<Epetra_Map>    augActiveSlaveTDofs_;
  Teuchos::RCP<Epetra_Map>    augActiveMasterNodes_;
  Teuchos::RCP<Epetra_Map>    augActiveMasterDofs_;

  Teuchos::RCP<Epetra_Map>    activeMasterNodes_;
  Teuchos::RCP<Epetra_Map>    activeMasterDofs_;
};
}


#endif /* CONTACT_AUGMENTED_INTERFACE_H_ */
