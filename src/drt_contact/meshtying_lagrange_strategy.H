/*!----------------------------------------------------------------------
\file meshtying_lagrange_strategy.H

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

#ifndef MESHTYING_LAGRANGE_STRATEGY_H
#define MESHTYING_LAGRANGE_STRATEGY_H

#include "meshtying_abstract_strategy.H"

using namespace std;
using namespace Teuchos;

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
    MtLagrangeStrategy(DRT::Discretization& discret,RCP<Epetra_Map> problemrowmap,
                       Teuchos::ParameterList params, vector<RCP<MORTAR::MortarInterface> > interface,
                       int dim, RCP<Epetra_Comm> comm, double alphaf, int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~MtLagrangeStrategy() {};

    /*!
    \brief Do mortar coupling in reference configuration
     
    Only do this ONCE for meshtying upon initialization!
    
    */
    void MortarCoupling(const RCP<Epetra_Vector> dis);

    /*!
    \brief Mesh intialization for rotational invariance
    
    Only do this ONCE for meshtying upon initialization!
    (only implemented for DUAL shape functions)
    
    */
    void MeshInitialization();
      
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
    void EvaluateMeshtying(RCP<LINALG::SparseOperator>& kteff,
                           RCP<Epetra_Vector>& feff, RCP<Epetra_Vector> dis);

    /*!
    \brief Solve linear system of saddle point type
   
    \param solver (in): LINALG::Solver object
    \param kdd (in): the displacement dof stiffness (upper left block)
    \param fd (in): the displacement dof r.h.s. (upper block)
    \param sold (in/out): the displacement dof solution increment
    \param dirichtoggle (in): toggle vector for dirichlet conditions
    \param numiter (in): current nonlinear iteration step
    
    This routine is called by the <tt>LinearSolve</tt> function in the (contact) structure time integration algorithm.
    The <tt>SaddlePointSolve</tt> problem is only used for saddle point problems (contact or meshtying). The solver must be
    a LINALG::Solver object with a block preconditioner (e.g. CheapSIMPLE). The fallback solver is used whenever there is no 
    contact and a pure (symmetric) structural problem has to be solved. Therefore per default the STRUCT SOLVER block parameters
    are used.
    */
    void SaddlePointSolve(LINALG::Solver& solver,
                          LINALG::Solver& fallbacksolver,
                          RCP<LINALG::SparseOperator> kdd,
                          RCP<Epetra_Vector>          fd,
                          RCP<Epetra_Vector>          sold,
                          RCP<Epetra_Vector>          dirichtoggle,
                          int numiter);

    /*!
    \brief Solve linear system

    \param solver (in): LINALG::Solver object
    \param kdd (in): the displacement dof stiffness (upper left block)
    \param fd (in): the displacement dof r.h.s. (upper block)
    \param sold (in/out): the displacement dof solution increment
    \param dirichtoggle (in): toggle vector for dirichlet conditions
    \param numiter (in): current nonlinear iteration step

    This routine is called by the <tt>LinearSolve</tt> function in the (contact) structure time integration algorithm.
    The <tt>Solve</tt> method is only used for condensed systems (no saddle point structure). Per default, the solver
    parameter is a <tt>LINALG::SOLVER</tt> object from the CONTACT SOLVER block in the BACI dat file, that is used for 
    contact problems. If there is no contact, the fallback solver is applied (per default: STRUCT SOLVER from BACI dat file).
    */
    void Solve(LINALG::Solver& solver,
               LINALG::Solver& fallbacksolver,
               RCP<LINALG::SparseOperator> kdd,
               RCP<Epetra_Vector>          fd,
               RCP<Epetra_Vector>          sold,
               RCP<Epetra_Vector>          dirichtoggle,
               int numiter);
    
    /*!
    \brief Recovery method

    We only recover the Lagrange multipliers here, which had been statically condensed during
    the setup of the global problem!
    
    */
    void Recover(RCP<Epetra_Vector> disi);
    
    //@}
    
    //! @name Empty functions (Penalty meshtying)
        
    // All these functions only have functionality in Penalty meshtying simulations,
    // thus they are defined as dserror here in the case of Lagrange meshtying..
    
    double ConstraintNorm() {return 0.0;}
    void InitializeUzawa(RCP<LINALG::SparseOperator>& kteff, RCP<Epetra_Vector>& feff) {}
    double InitialPenalty() {return 0.0;}   
    void ResetPenalty() {}
    void SaveReferenceState(const RCP<Epetra_Vector> dis) {}
    void UpdateAugmentedLagrange() {}
    void UpdateConstraintNorm(int uzawaiter = 0) {}
    
    //@}
    
  protected:

    // don't want = operator and cctor
    MtLagrangeStrategy operator = (const MtLagrangeStrategy& old);
    MtLagrangeStrategy(const MtLagrangeStrategy& old);

    RCP<LINALG::SparseMatrix> conmatrix_;   // constraint matrix for saddle point system
    RCP<LINALG::SparseMatrix> mhatmatrix_;  // product of global Mortar matrices inv(D)*M    
    
    RCP<Epetra_Vector>        fs_;          // slave side effective forces (needed for LM)
    RCP<LINALG::SparseMatrix> invd_;        // inverse of Mortar matrix D (needed for LM)
    RCP<LINALG::SparseMatrix> ksn_;         // stiffness block K_sn (needed for LM)
    RCP<LINALG::SparseMatrix> ksm_;         // stiffness block K_sm (needed for LM)
    RCP<LINALG::SparseMatrix> kss_;         // stiffness block K_ss (needed for LM)

}; // class MtLagrangeStrategy
} // namespace CONTACT


#endif  // #ifndef MESHTYING_LAGRANGE_STRATEGY_H
#endif  // #ifdef CCADISCRET
