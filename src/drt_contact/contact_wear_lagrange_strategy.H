/*!----------------------------------------------------------------------
\file contact_wear_lagrange_srategy.H

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              farah 09/13 |
 *----------------------------------------------------------------------*/
#ifndef CONTACT_WEAR_LAGRANGE_STRATEGY_H_
#define CONTACT_WEAR_LAGRANGE_STRATEGY_H_

/*----------------------------------------------------------------------*
 | header                                                   farah 09/13 |
 *----------------------------------------------------------------------*/
#include "contact_lagrange_strategy.H"
#include <Epetra_FEVector.h>



namespace CONTACT
{
// forward declarations
class WearInterface;


class WearLagrangeStrategy : public CoLagrangeStrategy
{

public:

  /*!
  \brief Standard Constructor

  */
  WearLagrangeStrategy(DRT::Discretization& probdiscret,
                     Teuchos::ParameterList params,
                     std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface, int dim,
                     Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof);

  /*!
  \brief Destructor

  */
  virtual ~WearLagrangeStrategy() {};


  /*!
  \brief Evaluate frictional contact

  */
  void EvaluateFriction(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff);


  /*!
  \brief Condense discr. wear and lm. for frictional contact

  */
  void CondenseWearDiscr(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff,
                         Teuchos::RCP<Epetra_Vector>& gact);

  /*!
  \brief Condense lm. for frictional contact with explicit/implicit wear algorithm

  */
  void CondenseWearImplExpl(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff,
                         Teuchos::RCP<Epetra_Vector>& gact);

  /*!
  \brief Prepare SaddlePointSystem

  */
  void PrepareSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff);
  /*!
  \brief Initialize and evaluate Mortar stuff for the next Newton step

  This method first checks if we are dealing with self contact and updates
  the interface slave and master sets if so. Then it resets the global
  Mortar matrices D and M and the global gap vector g accordingly.

  The nodal quantites computed in InitEvalInterface() are then assembled
  to global matrices and vectors respectively. No setup of the global system
  is to be done here yet, so there is no need to pass in the effective
  stiffness K or the effective load vector f.

  */
  void InitMortar();
  void AssembleMortar();

  /*!
  \brief Initialize general contact variables for next Newton step

  For a lagrangian strategy this includes the global normal / tangent matrices N and T,
  the global derivative matrices S and P and Tresca friction matrix L + vector r.

  */
  void Initialize();

  /*!
  \brief Recovery method

  We only recover the Lagrange multipliers here, which had been
  statically condensated during the setup of the global problem!
  Optionally satisfaction or violation of the contact boundary
  conditions can be checked, too.

  */
  void Recover(Teuchos::RCP<Epetra_Vector> disi);

  /*!
  \brief Redistribute all contact interfaces in parallel

  We hand in the current global displacement state so that
  a contact search can be performed and set state called.

  */
  void RedistributeContact(Teuchos::RCP<Epetra_Vector> dis);

  /*!
  \brief Setup this strategy object (maps, vectors, etc.)

  All global maps and vectors are initialized by collecting
  the necessary information from all interfaces. In the case
  of a parallel redistribution, this method is called again
  to re-setup the above mentioned quantities. In this case
  we set the input parameter redistributed=TRUE. Moreover,
  when called for the first time (in the constructor) this
  method is given the input parameter init=TRUE to account
  for initialization of the active set.

  */
  void SetupWear(bool redistributed, bool init);

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
                        Teuchos::RCP<LINALG::MapExtractor>   dbcmaps,
                        int numiter);

  /*!
  \brief Evaluate wear vector

  Evaluates the unweighted wear vector.
  Refer also to the Semesterarbeit of Karl Wichmann 2010

  */
  void OutputWear();


  /*!
  \brief Perform a write restart

  A write restart is initiated by the contact manager. However, the manager has no
  direct access to the nodal quantities. Hence, a portion of the restart has to be
  performed on the level of the contact algorithm, for short: here's the right place.

  */
  void DoWriteRestart(Teuchos::RCP<Epetra_Vector>& activetoggle,
                      Teuchos::RCP<Epetra_Vector>& sliptoggle,
                      Teuchos::RCP<Epetra_Vector>& weightedwear,
                      Teuchos::RCP<Epetra_Vector>& realwear,
                      bool forcedrestart = false);

  /*!
  \brief Perform a write restart

  A write restart is initiated by the contact manager. However, the manager has no
  direct access to the nodal quantities. Hence, all the restart action has to be
  performed on the level of the contact algorithm, for short: here's the right place.

  */
  void DoReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis);


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


  void UpdateWearDiscretIterate(bool store);
  void UpdateWearDiscretAccumulation(bool wearaccumulation);
  /*!
  \brief Return vector of wear (t_n+1)

  */
  Teuchos::RCP<Epetra_Vector> ContactWear() { return wearoutput_; }   // for slave side
  Teuchos::RCP<Epetra_Vector> ContactWear2() { return wearoutput2_; } // for master side

  /*!
  \brief Return wear interfaces

  */
  std::vector<Teuchos::RCP<CONTACT::WearInterface> > WearInterfaces() { return interface_; }

  const Teuchos::RCP<Epetra_Map> MasterSlipNodes() {return gmslipnodes_;};
  const Teuchos::RCP<Epetra_Map> MasterActiveNodes() {return gmactivenodes_;};

protected:
  // don't want = operator and cctor
  WearLagrangeStrategy operator = (const WearLagrangeStrategy& old);
  WearLagrangeStrategy(const WearLagrangeStrategy& old);

  std::vector<Teuchos::RCP<CONTACT::WearInterface> > interface_;

  // implicit wear algorithm
  Teuchos::RCP<LINALG::SparseMatrix> wlinmatrix_;   // global Matrix Wg containing wear-lm derivatives
  Teuchos::RCP<LINALG::SparseMatrix> wlinmatrixsl_; // global Matrix Wsl containing wear-lm slip derivatives
  Teuchos::RCP<LINALG::SparseMatrix> wlinmatrixst_; // global Matrix Wst containing wear-lm stick derivatives

  // both-sided wear weak dirich cond
  Teuchos::RCP<LINALG::SparseMatrix> d2matrix_;     // global Mortar matrix D2

  Teuchos::RCP<Epetra_Map> gminvolvednodes_;        // global involved master node row map (of all interfaces)
  Teuchos::RCP<Epetra_Map> gminvolveddofs_;         // global involved master dof row map (of all interfaces)
  Teuchos::RCP<Epetra_Map> gslipn_;                 // global row map of matrix N for slip dofs (of all interfaces)
  Teuchos::RCP<Epetra_Map> gwinact_;                // global row map of matrix N for slip dofs (of all interfaces)
  Teuchos::RCP<Epetra_Map> gmslipn_;                // global row map of matrix N for slip dofs (of all interfaces)
  Teuchos::RCP<Epetra_Map> gwminact_;               // global row map of matrix N for slip dofs (of all interfaces)

  Teuchos::RCP<Epetra_Map> gwmdofrowmap_;           // global master wear dof row map (of all interfaces) -active
  Teuchos::RCP<Epetra_Map> gwdofrowmap_;            // global slave wear dof row map (of all interfaces) -active
  Teuchos::RCP<Epetra_Map> gsdofnrowmap_;           // global slave wear dof row map (of all interfaces)
  Teuchos::RCP<Epetra_Map> gmdofnrowmap_;           // global master wear dof row map (of all interfaces)
  Teuchos::RCP<Epetra_Map> galldofnrowmap_;           // global master wear dof row map (of all interfaces)
  Teuchos::RCP<Epetra_Map> gwalldofrowmap_;         // all
  Teuchos::RCP<Epetra_Map> gmslipnodes_;            // global master slip nodes
  Teuchos::RCP<Epetra_Map> gmactivenodes_;          // global master active nodes

  Teuchos::RCP<Epetra_Vector> wearoutput_;          // vector of unweighted wear at t_n+1  -- slave
  Teuchos::RCP<Epetra_Vector> wearoutput2_;         // vector of unweighted wear at t_n+1  -- master
  Teuchos::RCP<Epetra_Vector> wearvector_;          // global weighted wear vector w

  int  maxdofwear_;                                 // highest dof number in problem discretization

  bool wearimpl_;                                   // weartype: implicit
  bool weardiscr_;                                  // bool for wear with own discretization
  bool wearbothdiscr_;                              // bool for both-sided discrete wear
  bool weartimescales_;                             // bool for different time scales
  bool sswear_;                                     // bool steady state wear

  // discrete wear algorithm (SLAVE)
  Teuchos::RCP<LINALG::SparseMatrix> twmatrix_;     // global Mortar wear matrix T
  Teuchos::RCP<LINALG::SparseMatrix> ematrix_;      // global Mortar wear matrix E
  Teuchos::RCP<LINALG::SparseMatrix> eref_;         // global Mortar wear matrix E
  Teuchos::RCP<LINALG::SparseMatrix> lintdis_;      // Lin T w.r.t. displ: Lin(T*n*lm)
  Teuchos::RCP<LINALG::SparseMatrix> lintlm_;       // Lin T w.r.t. lm: (T*n)
  Teuchos::RCP<LINALG::SparseMatrix> linedis_;      // Lin E w.r.t. displ: Lin(E*w)
  Teuchos::RCP<LINALG::SparseMatrix> linslipW_;     // global matrix containing derivatives (LM) of slip condition
  Teuchos::RCP<Epetra_Vector> inactiveWearRhs_;     // inactive wear rhs: -w_i
  Teuchos::RCP<Epetra_Vector> WearCondRhs_;         // rhs wear condition: -E*w_i + k*T*n*lm_i

  // discrete wear algorithm (MASTER)
  Teuchos::RCP<LINALG::SparseMatrix> twmatrixM_;     // global Mortar wear matrix T
  Teuchos::RCP<LINALG::SparseMatrix> ematrixM_;      // global Mortar wear matrix E
  Teuchos::RCP<LINALG::SparseMatrix> lintdisM_;      // Lin T w.r.t. displ: Lin(T*n*lm)
  Teuchos::RCP<LINALG::SparseMatrix> lintlmM_;       // Lin T w.r.t. lm: (T*n)
  Teuchos::RCP<LINALG::SparseMatrix> linedisM_;      // Lin E w.r.t. displ: Lin(E*w)
  Teuchos::RCP<LINALG::SparseMatrix> linslipWM_;     // global matrix containing derivatives (LM) of slip condition
  Teuchos::RCP<Epetra_FEVector> inactiveWearRhsM_;   // inactive wear rhs: -w_i
  Teuchos::RCP<Epetra_FEVector> WearCondRhsM_;       // rhs wear condition: -E*w_i + k*T*n*lm_i

  // matrix blocks for recovering
  Teuchos::RCP<LINALG::SparseMatrix> dnblock_;
  Teuchos::RCP<LINALG::SparseMatrix> dmblock_;
  Teuchos::RCP<LINALG::SparseMatrix> diblock_;
  Teuchos::RCP<LINALG::SparseMatrix> dablock_;
  Teuchos::RCP<Epetra_Vector> fw_;

  Teuchos::RCP<Epetra_Map> gidofs_;


}; //class


}


#endif /* CONTACT_WEAR_H_ */
