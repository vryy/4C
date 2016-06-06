/*---------------------------------------------------------------------*/
/*!
\file contact_penalty_strategy.H

\brief Penalty contact solving strategy: The contact constrains are enforced
       by a penalty formulation.

\level 2

\maintainer Philipp Farah, Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_PENALTY_STRATEGY_H
#define CONTACT_PENALTY_STRATEGY_H

#include "contact_abstract_strategy.H"


namespace CONTACT
{

// forward declaration
//class WearInterface;
/*!
 \brief Contact solving strategy with regularization of Lagrangian multipliers,
 also known as Penalty Method or regularization. An Augmented Lagrangian version
 based on the Uzawa algorithm is included, too.

 This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
 For a more general documentation of the involved functions refer to CoAbstract Strategy.

 Refer also to the Semesterarbeit of Bernd Budich, 2009

 \author popp (popp@lnm.mw.tum.de)
 */
class CoPenaltyStrategy : public CoAbstractStrategy
{
  public:

  //! Standard constructor
  CoPenaltyStrategy(
      const Epetra_Map* DofRowMap,
      const Epetra_Map* NodeRowMap,
      Teuchos::ParameterList params,
      std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
      int dim,
      Teuchos::RCP<Epetra_Comm> comm,
      double alphaf,
      int maxdof);

  //! Shared data constructor
  CoPenaltyStrategy(
      const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
      const Epetra_Map* DofRowMap,
      const Epetra_Map* NodeRowMap,
      Teuchos::ParameterList params,
      std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
      int dim,
      Teuchos::RCP<Epetra_Comm> comm,
      double alphaf,
      int maxdof);

  /*!
  \brief Destructor

  */
  virtual ~CoPenaltyStrategy() {};

  //! @name Access methods

  /*!
  \brief Return L2-norm of active constraints

  */
  double ConstraintNorm() const { return constrnorm_;}

  /*!
  \brief Return L2-norm of slip constraints

  */
  double ConstraintNormTan() { return constrnormtan_;}


  /*!
  \brief Return initial penalty parameter for non-penetration

  */
  double InitialPenalty() { return initialpenalty_;}

  /*!
  \brief Return initial penalty parameter for tangential direction

  */
  double InitialPenaltyTan() { return initialpenaltytan_;}


  //@}

  //! @name Evaluation methods

  /*!
  \brief Save nodal kappa-coefficients

  Before starting with the time integration, we have to calculate a nodal scaling factor,
  which will compensate the different integration area for computing the nodal weighted
  gap. Omitting this scaling, nodes on edges or boundaries would have a smaller weighted
  gap, even in case of a uniform physical gap. Hence, this scaling is of crucial importance
  for a penalty strategy since the weighted gap determines the lagrangian multipliers.

  */
  void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis);

  /*!
  \brief Evaluate relative movement of contact bodies in predictor

  This is a tiny control routine for evaluating the relative movement of
  contact bodies in the predictor of an implicit time integration scheme.
  This evaluation (resetting) is ONLY necessary for penalty strategy and
  Uzawa augmented lagrange strategy, thus this tiny routine here.

  */

  void EvaluateRelMovPredict();

  /*!
  \brief Initialize general contact variables for next Newton step

  For a penalty strategy this involves the derivative matrix for the regularized lagrange multipliers.

  */
  void Initialize();

  /*!
  \brief Evaluate contact

  For a penalty strategy this includes the evaluation of regularized forces
  in normal and tangential direction and results in a simple addition of extra
  stiffness contributions to kteff and extra contact forces to feff.

  */
  void EvaluateContact(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff);

  /*!
  \brief Evaluate frictional contact

  This includes the evaluation of of the frictional contact forces.

  */
  void EvaluateFriction(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff);

  /*!
  \brief Reset penalty parameter to intial value

  When applying an Uzawa Augmented Lagrangian version of the penalty approach,
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
  \brief Initialize Uzawa step


  This method is called at the beginning of the second, third, ... Uzawa
  iterarion in order to create an of an out-of-balance force again. First,
  the contact force and stiffness terms are removed from feff and kteff.
  Then the LM and derivatives are updated (Uzawa AugmentedLagrange) and the new
  contact forces and stiffness terms are created by calling Initialize()
  and finally Evaluate().

  */
  void InitializeUzawa(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff);

  /*!
  \brief Compute L2-norm of active constraints

  In a classical penalty approach, the constraint norm is only monitored.
  When applying an Uzawa Augmented Lagrangian version, the constraint norm is the
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
  \brief Apply the GPTS contact force and stiffness to the global system

  */
  void ApplyGPTSforces(Teuchos::RCP<LINALG::SparseOperator> kteff, Teuchos::RCP<Epetra_Vector> feff);

  virtual void SetParentState(const std::string& statename,
      const Teuchos::RCP<Epetra_Vector> vec,
      const Teuchos::RCP<DRT::Discretization> dis);
  //@}

  //! @name Empty functions (Lagrange contact)

  // All these functions only have functionality in Lagrange contact simulations,
  // thus they are defined empty here in the case of Penalty contact.

  //! Get the active node row map of the previous Newton step
  virtual Teuchos::RCP<const Epetra_Map> GetOldActiveRowNodes() const
  { dserror("Currently unsupported for the penalty strategy!"); return Teuchos::null; };
  virtual Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const
  { dserror("Currently unsupported for the penalty strategy!"); return Teuchos::null; };
  bool ActiveSetSemiSmoothConverged() const {return true;}
  bool ActiveSetConverged() {return true;}
  int ActiveSetSteps() {return 0;}
  void ResetActiveSet() {}
  void Recover(Teuchos::RCP<Epetra_Vector> disi) {}
  void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd, Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<LINALG::MapExtractor> dbcmaps, int numiter,Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) { std::cout << __FILE__ << ":" << __LINE__ << "WARNING: saddle point system makes no sense here" << std::endl;}
  void UpdateDisplacementsAndLMincrements(Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Epetra_Vector> blocksol) { std::cout << __FILE__ << ":" << __LINE__ << "WARNING: saddle point system makes no sense here" << std::endl;}
  void EvalConstrRHS(){}
  void UpdateActiveSet() {}
  void UpdateActiveSetSemiSmooth() {}

protected:

  // don't want = operator and cctor
  CoPenaltyStrategy operator = (const CoPenaltyStrategy& old);
  CoPenaltyStrategy(const CoPenaltyStrategy& old);

  Teuchos::RCP<LINALG::SparseMatrix> linzmatrix_;        // global matrix LinZ with derivatives of LM
  double                             constrnorm_;        // L2-norm of normal contact constraints
  double                             constrnormtan_;     // L2-norm of tangential contact constraints
  double                             initialpenalty_;    // initial penalty parameter
  double                             initialpenaltytan_; // initial tangential penalty parameter

}; // class CoPenaltyStrategy
}  // namespace CONTACT

#endif  // #ifndef CONTACT_PENALTY_STRATEGY_H
