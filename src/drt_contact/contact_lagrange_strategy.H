/*!----------------------------------------------------------------------
\file contact_lagrange_strategy.H

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
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#ifndef CONTACT_LAGRANGE_STRATEGY_H
#define CONTACT_LAGRANGE_STRATEGY_H

#include "contact_abstract_strategy.H"


namespace CONTACT
{

/*!
 \brief Contact solving strategy with (standard/dual) Lagrangian multipliers.

 This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
 For a more general documentation of the involved functions refer to CoAbstract Strategy.

 \author popp (popp@lnm.mw.tum.de)
 */
class CoLagrangeStrategy : public CoAbstractStrategy
{
  public:

    /*!
    \brief Standard Constructor

    */
    CoLagrangeStrategy(DRT::Discretization& probdiscret,
                       Teuchos::ParameterList params,
                       std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface, int dim,
                       Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~CoLagrangeStrategy() {};
    
    //! @name Access methods
    
    /*!
    \brief Return convergence status of semi-smooth active set search

    If this Lagrange contact strategy is not based on a semi-smooth
    Newton approach, but on a fixed-point approach with two nested
    loops, then this method simply returns true, of course. Convergence
    of the active set is monitored with the flag activesetconv_ in
    this case and activesetssconv_ is meaningless.
    
    */
    bool ActiveSetSemiSmoothConverged()
    {
      bool semismooth = DRT::INPUT::IntegralValue<int>(Params(),"SEMI_SMOOTH_NEWTON");
      if (semismooth) return activesetssconv_;
      else            return true;
    }
    
    /*!
    \brief Return convergence status of fixed-point active set search

    If this Lagrange contact strategy is based on a semi-smooth
    Newton approach and not on a fixed-point approach with two nested
    loops, then this method simply returns true, of course. Convergence
    of the active set is monitored with the flag activesetssconv_ in
    this case and activesetconv_ is meaningless.
    
    */
    bool ActiveSetConverged()
    {
      bool semismooth = DRT::INPUT::IntegralValue<int>(Params(),"SEMI_SMOOTH_NEWTON");
      if (!semismooth) return activesetconv_;
      else             return true;
    }
    
    /*!
    \brief Return no. of fixed-point active sets in this time step

    */
    int ActiveSetSteps() { return activesetsteps_; }
    
    //@}

    //! @name Evaluation methods
    
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
    void EvaluateContact(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff);

    /*!
    \brief Evaluate frictional contact
     
    */
    void EvaluateFriction(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff);

    /*!
    \brief Solve linear system of saddle point type
   
    \param solver (in): LINALG::Solver object
    \param solver (in): LINALG::Solver object (fallback solver)
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
                          Teuchos::RCP<LINALG::SparseOperator> kdd,
                          Teuchos::RCP<Epetra_Vector>          fd,
                          Teuchos::RCP<Epetra_Vector>          sold,
                          Teuchos::RCP<Epetra_Vector>          dirichtoggle,
                          int numiter);
    

    /*!
    \brief Recovery method

    We only recover the Lagrange multipliers here, which had been
    statically condensated during the setup of the global problem!
    Optionally satisfaction or violation of the contact boundary
    conditions can be checked, too.
    
    */
    void Recover(Teuchos::RCP<Epetra_Vector> disi);

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
    
    /*!
    \brief Reset active set status for next time step

    */
    void ResetActiveSet()
    {
      activesetssconv_ = false;
      activesetconv_ = false;
      activesetsteps_ = 1;
      return;
    }
    
    /*!
    \brief Return matrix T

    */
    Teuchos::RCP<LINALG::SparseMatrix> TMatrix() { return tmatrix_; }
    
    //@}
    
    //! @name Empty methods
      
    /*!
    \brief Empty methods only relevant for other strategies
    
    For a Lagrange strategy these are functions without functionality.
    Call them whenever you like.
    
    */
    double ConstraintNorm() {return 0.0;}
    void EvaluateRelMovPredict() {}
    double InitialPenalty() {return 0.0;}
    void InitializeUzawa(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) {}
    void ResetPenalty() {}
    void SaveReferenceState(const Teuchos::RCP<Epetra_Vector> dis) {}
    void UpdateAugmentedLagrange() {}
    void UpdateConstraintNorm(int uzawaiter = 0) {}
    
    //@}

  protected:

    // don't want = operator and cctor
    CoLagrangeStrategy operator = (const CoLagrangeStrategy& old);
    CoLagrangeStrategy(const CoLagrangeStrategy& old);

    bool activesetssconv_;    // convergence flag for semi-smooth active set search
    bool activesetconv_;      // convergence flag for fixed-point active set search      
    int  activesetsteps_;     // number of fixed-point active set steps in this time step
    
    Teuchos::RCP<LINALG::SparseMatrix> mhatmatrix_; // product of global Mortar matrices inv(D)*M    

    Teuchos::RCP<LINALG::SparseMatrix> tmatrix_; // global Matrix T containing active node tangents
    Teuchos::RCP<LINALG::SparseMatrix> smatrix_; // global Matrix S containing normal+D+M derivatives
    Teuchos::RCP<LINALG::SparseMatrix> pmatrix_; // global Matrix P containing tangent derivatives

    Teuchos::RCP<LINALG::SparseMatrix> wlinmatrix_;   // global Matrix Wg containing wear-lm derivatives
    Teuchos::RCP<LINALG::SparseMatrix> wlinmatrixsl_; // global Matrix Wsl containing wear-lm slip derivatives
    Teuchos::RCP<LINALG::SparseMatrix> wlinmatrixst_; // global Matrix Wst containing wear-lm stick derivatives


    Teuchos::RCP<Epetra_Vector> fs_;          // slave side effective forces (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> invd_; // inverse of Mortar matrix D (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> ksn_;  // stiffness block K_sn (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> ksm_;  // stiffness block K_sm (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> kss_;  // stiffness block K_ss (needed for LM)

    Teuchos::RCP<LINALG::SparseMatrix>  linslipLM_;   // global matrix containing derivatives (LM) of slip condition
    Teuchos::RCP<LINALG::SparseMatrix>  linslipDIS_;  // global matrix containing derivatives (DIS) of slip condition
    Teuchos::RCP<Epetra_Vector>         linslipRHS_;  // r.h.s vector friction slip nodes
    Teuchos::RCP<LINALG::SparseMatrix>  linstickLM_;  // global matrix containing derivatives (LM) of slip condition
    Teuchos::RCP<LINALG::SparseMatrix>  linstickDIS_; // global matrix containing derivatives (DIS) of stick condition
    Teuchos::RCP<Epetra_Vector>         linstickRHS_; // r.h.s vector for friction stick condition
    
    Teuchos::RCP<Epetra_Map> zigzagone_;   // active node set of last active set try
    Teuchos::RCP<Epetra_Map> zigzagtwo_;   // active node set of second-last active set try
    Teuchos::RCP<Epetra_Map> zigzagthree_; // active node set of third-last active set try

}; // class CoLagrangeStrategy
} // namespace CONTACT


#endif  // #ifndef CONTACT_LAGRANGE_STRATEGY_H
