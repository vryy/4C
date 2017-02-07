/*---------------------------------------------------------------------*/
/*!
\file contact_nitsche_strategy.H

\brief Nitsche contact solving strategy

\level 3

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_NITSCHE_STRATEGY_H
#define CONTACT_NITSCHE_STRATEGY_H

#include "contact_abstract_strategy.H"
namespace CONTACT
{

/*!
 \brief Contact solving strategy with Nitsche's method.

 This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
 For a more general documentation of the involved functions refer to CoAbstract Strategy.

 */
class CoNitscheStrategy : public CoAbstractStrategy
{
public:

  //! Standard constructor
  CoNitscheStrategy(
      const Epetra_Map* DofRowMap,
      const Epetra_Map* NodeRowMap,
      Teuchos::ParameterList params,
      std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
      int dim,
      Teuchos::RCP<Epetra_Comm> comm,
      double alphaf,
      int maxdof)
: CoAbstractStrategy(Teuchos::rcp(new CONTACT::AbstractStratDataContainer()),
    DofRowMap,NodeRowMap,params,interface,dim,comm,alphaf,maxdof){}

  //! Shared data constructor
  CoNitscheStrategy(
      const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
      const Epetra_Map* DofRowMap,
      const Epetra_Map* NodeRowMap,
      Teuchos::ParameterList params,
      std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
      int dim,
      Teuchos::RCP<const Epetra_Comm> comm,
      double alphaf,
      int maxdof)
  : CoAbstractStrategy(data_ptr,DofRowMap,NodeRowMap,params,interface,dim,comm,alphaf,maxdof){}

  virtual void ApplyForceStiffCmt(
      Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<LINALG::SparseOperator>& kt,
      Teuchos::RCP<Epetra_Vector>& f, const int step, const int iter,
      bool predictor);

  virtual void DoReadRestart(IO::DiscretizationReader& reader,
      Teuchos::RCP<const Epetra_Vector> dis,
      Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr);

  virtual bool IsSaddlePointSystem() const {return false;}

  virtual bool IsCondensedSystem() const {return false;}

  void Evaluate(
      CONTACT::ParamsInterface& cparams,
      const std::vector<Teuchos::RCP<const Epetra_Vector> >& eval_vec);

  void Integrate(CONTACT::ParamsInterface& cparams);

  virtual Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
     const enum STR::VecBlockType& bt) const;

  virtual Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
      const enum STR::MatBlockType& bt) const ;

  /*! \brief Setup this strategy object (maps, vectors, etc.)

   derived from contact abstract strategy.
   The Nitsche strategy does not have
    */
  virtual void Setup(bool redistributed, bool init);

  virtual void UpdateTraceIneqEtimates();

  /*! \brief Get dirichlet B.C. status and store into CoNodes

   This is called once at the beginning of the simulation
   to set the D.B.C. status in each CNode.

   \param dbcmaps (in): MapExtractor carrying global dbc map */
  virtual void StoreDirichletStatus(Teuchos::RCP<const LINALG::MapExtractor> dbcmaps)
  {/* we don't care about dirichlet for now */};
  virtual void Update(Teuchos::RCP<const Epetra_Vector> dis);
  virtual void EvaluateReferenceState(Teuchos::RCP<const Epetra_Vector> vec);
  virtual void DoWriteRestart(std::map<std::string,Teuchos::RCP<Epetra_Vector> >& restart_vectors,
        bool forcedrestart = false) const
  {/* nothing stored in nitsche strategy that would need to be written */};
  virtual void OutputStresses()
  {/* nothing stress output in nitsche strategy yet */};
  virtual void ReconnectParentElements();

  /*!
\brief Destructor

   */
  virtual ~CoNitscheStrategy() {};

  virtual void SetState(const enum MORTAR::StateType& statename,
        const Epetra_Vector& vec);

  virtual void SetParentState(const std::string& statename,
      const Epetra_Vector& vec);

  virtual Teuchos::RCP<const Epetra_Vector> GetLagrMultN(const bool& redist) const
      {;return Teuchos::null;}
  virtual Teuchos::RCP<const Epetra_Vector> GetLagrMultNp(const bool& redist) const
      {return Teuchos::null;}
  virtual Teuchos::RCP<Epetra_Vector> LagrMultOld()
      {return Teuchos::null;}
  virtual Teuchos::RCP<const Epetra_Map> LMDoFRowMapPtr(const bool& redist) const
      {return Teuchos::null;}
  // All these functions only have functionality in Lagrange contact simulations,
  // thus they are defined empty here in the case of Penalty contact.

  //! Get the active node row map of the previous Newton step
  virtual Teuchos::RCP<const Epetra_Map> GetOldActiveRowNodes() const
    { return Teuchos::null; };
  virtual Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const
    { dserror("Currently unsupported for the penalty strategy!"); return Teuchos::null; };
  virtual bool IsNitsche() const {return true;}
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
  void EvaluateRelMovPredict() {}
  void ModifyPenalty() {}
  void UpdateUzawaAugmentedLagrange() {}
  void UpdateConstraintNorm(int uzawaiter = 0) {}
  void Initialize(){};
  void EvaluateContact(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
  {dserror("not supported in this strategy");}
  void EvaluateFriction(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
  {dserror("not supported in this strategy");}
  void InitializeUzawa(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) {}
  void ResetPenalty() {}
  void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) {}
  double InitialPenalty() {return 0.0;}
  double ConstraintNorm() const {return 0.0;}

protected:

  // don't want = operator and cctor
  CoNitscheStrategy operator = (const CoNitscheStrategy& old);
  CoNitscheStrategy(const CoNitscheStrategy& old);

  Teuchos::RCP<Epetra_Vector> curr_state_;
  bool                        curr_state_eval_;
};
}
#endif //  CONTACT_NITSCHE_STRATEGY_H
