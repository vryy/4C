/*!----------------------------------------------------------------------
\file drt_contact_lagrange_strategy.H

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
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef DRT_CONTACT_LAGRANGE_STRATEGY_H
#define DRT_CONTACT_LAGRANGE_STRATEGY_H

#include "Teuchos_RefCountPtr.hpp"
#include "drt_contact_abstract_strategy.H"
#include "../drt_lib/drt_discret.H"

using namespace std;
using namespace Teuchos;

namespace CONTACT
{

/*!
 \brief Contact solving strategy with (standard/dual) Lagrangian multipliers.

 This is a specialization of the abstract contact algorithm as defined in AbstractStrategy.
 For a more general documentation of the involved functions refer to Abstract Strategy.

 \author popp (popp@lnm.mw.tum.de)
 */
class LagrangeStrategy : public AbstractStrategy
{
  public:

    /*!
    \brief Standard Constructor

    */
    explicit LagrangeStrategy(RCP<Epetra_Map> problemrowmap, Teuchos::ParameterList params,
                              vector<RCP<CONTACT::Interface> > interface, int dim,
                              RCP<Epetra_Comm> comm, double alphaf);

    /*!
    \brief Destructor

    */
    virtual ~LagrangeStrategy() {};

    /*!
    \brief Initialize general contact variables for next Newton step
     
    For a lagrangian strategy this includes the global normal / tangent matrices N and T, 
    the global derivative matrices S and P and Tresca friction matrix L + vector r.
    
    */
    void Initialize();

    /*!
    \brief Evaluate contact
     
    For a lagrangian strategy this involves heavy modification to the initial kteff and feff.
    Hence, they are infact build from scratch here.
    The application of modifications to groups of dofs (slave, master, active etc.)
    results in some matrix and vector splitting and a lot of matrix-vector calculation in here!
    
    */
    void EvaluateContact(RCP<LINALG::SparseOperator>& kteff, RCP<Epetra_Vector>& feff);

    /*!
    \brief Evaluate frictional contact
     
    */
    void EvaluateFriction(RCP<LINALG::SparseOperator>& kteff, RCP<Epetra_Vector>& feff);

    /*!
    \brief Recovery method

    We only recover the Lagrange multipliers here, which had been
    statically condensated during the setup of the global problem!
    Optionally satisfaction or violation of the contact boundary
    conditions can be checked, too.
    
    */
    void Recover(RCP<Epetra_Vector> disi);

    /*!
    \brief Update active set and check for convergence
     
    In this function we loop over all interfaces and then over all
    slave nodes to check, whether the assumption of them being active
    or inactive respectively has been correct. If a single node changes
    state, the active set is adapted accordingly and the convergence
    flag is kept on false. Here we have the "standard" case of two
    nested iteration loops, and as a consequence this method is
    called AFTER convergence of the inner Newton iteration. If there
    is a change in the active set, another full Newton iteration has
    to be performed for the current time / load step.
     
    */
    void UpdateActiveSet();

    /*!
    \brief Update active set and check for convergence
     
    In this function we loop over all interfaces and then over all
    slave nodes to check, whether the assumption of them being active
    or inactive respectively has been correct. If a single node changes
    state, the active set is adapted accordingly and the convegence
    flag is kept on false. Here we have the semi-smooth Newton case
    with one combined iteration loop for active set search and large
    deformations. As a consequence this method is called AFTER each
    (not yet converged) Newton step. If there is a change in the active
    set or the residual and disp norm are still above their limits,
    another Newton step has to be performed.

    */
    void UpdateActiveSetSemiSmooth();

  protected:

    RCP<LINALG::SparseMatrix> mhatmatrix_; // product of global Mortar matrices inv(D)*M    
    RCP<LINALG::SparseMatrix> mhata_;      // active part of mhatmatrix_ (needed for HüeberCheck)

    RCP<LINALG::SparseMatrix> nmatrix_; // global Matrix N containing active node normals
    RCP<LINALG::SparseMatrix> tmatrix_; // global Matrix T containing active node tangents
    RCP<LINALG::SparseMatrix> smatrix_; // global Matrix S containing normal+D+M derivatives
    RCP<LINALG::SparseMatrix> pmatrix_; // global Matrix P containing tangent derivatives

    RCP<Epetra_Vector> fs_;          // slave side effective forces (needed for LM)
    RCP<LINALG::SparseMatrix> invd_; // inverse of Mortar matrix D (needed for LM)
    RCP<LINALG::SparseMatrix> ksn_;  // stiffness block K_sn (needed for LM)
    RCP<LINALG::SparseMatrix> ksm_;  // stiffness block K_sm (needed for LM)
    RCP<LINALG::SparseMatrix> kss_;  // stiffness block K_ss (needed for LM)

    RCP<LINALG::SparseMatrix>  lmatrix_;     // global matrix L for Tresca friction slip nodes
    RCP<Epetra_Vector>         r_;           // r.h.s vector for Tresca friction slip nodes
    RCP<Epetra_Vector>         rstick_;      // r.h.s vector for fricton stick nodes
    RCP<LINALG::SparseMatrix>  linslipLM_;   // global matrix containing derivatives (LM) of slip condition
    RCP<LINALG::SparseMatrix>  linslipDIS_;  // global matrix containing derivatives (DIS) of slip condition
    RCP<Epetra_Vector>         linslipRHS_;  // r.h.s vector for Tresca friction slip nodes
    RCP<LINALG::SparseMatrix>  linstickLM_;  // global matrix containing derivatives (LM) of slip condition
    RCP<LINALG::SparseMatrix>  linstickDIS_; // global matrix containing derivatives (DIS) of stick condition
    RCP<Epetra_Vector>         linstickRHS_; // r.h.s vector for stick condition
    
    RCP<Epetra_Map> zigzagone_;   // active node set of last active set try
    RCP<Epetra_Map> zigzagtwo_;   // active node set of second-last active set try
    RCP<Epetra_Map> zigzagthree_; // active node set of third-last active set try

    RCP<Epetra_Map> zigzagslipone_;   // slip node set of last slip set try
    RCP<Epetra_Map> zigzagsliptwo_;   // slip node set of second-last slip set try
    RCP<Epetra_Map> zigzagslipthree_; // slip node set of third-last slip set try

}; // class LagrangeStrategy
} // namespace CONTACT


#endif  // #ifndef DRT_CONTACT_LAGRANGE_STRATEGY_H
#endif  // #ifdef CCADISCRET
