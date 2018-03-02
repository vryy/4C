/*---------------------------------------------------------------------*/
/*!
\file contact_nitsche_strategy_tsi.H

\brief Nitsche contact solving strategy

\level 3

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_NITSCHE_STRATEGY_TSI_H
#define CONTACT_NITSCHE_STRATEGY_TSI_H

#include "contact_nitsche_strategy.H"

namespace ADAPTER {class Coupling;}

namespace CONTACT
{

/*!
 \brief Contact solving strategy with Nitsche's method.

 This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
 For a more general documentation of the involved functions refer to CoAbstract Strategy.

 */
class CoNitscheStrategyTsi : public CoNitscheStrategy
{
public:

  //! Standard constructor
  CoNitscheStrategyTsi(
      const Epetra_Map* DofRowMap,
      const Epetra_Map* NodeRowMap,
      Teuchos::ParameterList params,
      std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
      int dim,
      Teuchos::RCP<Epetra_Comm> comm,
      double alphaf,
      int maxdof)
: CoNitscheStrategy(DofRowMap,NodeRowMap,params,interface,dim,comm,alphaf,maxdof)
  ,fix_redistribution_(true)
{}

  //! Shared data constructor
  CoNitscheStrategyTsi(
      const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
      const Epetra_Map* DofRowMap,
      const Epetra_Map* NodeRowMap,
      Teuchos::ParameterList params,
      std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
      int dim,
      Teuchos::RCP<const Epetra_Comm> comm,
      double alphaf,
      int maxdof)
  : CoNitscheStrategy(data_ptr,DofRowMap,NodeRowMap,params,interface,dim,comm,alphaf,maxdof)
  ,fix_redistribution_(true)
  {}

  virtual void ApplyForceStiffCmt(
      Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<LINALG::SparseOperator>& kt,
      Teuchos::RCP<Epetra_Vector>& f, const int step, const int iter,
      bool predictor){dserror("not implemented");}

  void Integrate(CONTACT::ParamsInterface& cparams);

  /*! \brief Setup this strategy object (maps, vectors, etc.)

   derived from contact abstract strategy.
   The Nitsche strategy does not have
    */
  virtual void Setup(bool redistributed, bool init);

  virtual void UpdateTraceIneqEtimates();

  /*!
\brief Destructor

   */
  virtual ~CoNitscheStrategyTsi() {};

  virtual void SetState(const enum MORTAR::StateType& statename,
        const Epetra_Vector& vec);

  virtual void SetParentState(const enum MORTAR::StateType& statename,
      const Epetra_Vector& vec);

  void SetCoupling(Teuchos::RCP<ADAPTER::Coupling> coupST){coupST_=coupST;};

  /*!
   \brief Set time integration parameter from Thermo time integration

   */
  void SetAlphafThermo(const Teuchos::ParameterList& tdyn);

  virtual Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
     const enum DRT::UTILS::VecBlockType& bt) const;

  virtual Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
      const enum DRT::UTILS::MatBlockType& bt,
      const CONTACT::ParamsInterface* cparams=NULL) const ;

  virtual void RedistributeContact(Teuchos::RCP<const Epetra_Vector> dis)
  {
    if (fix_redistribution_) return;
    CONTACT::CoAbstractStrategy::RedistributeContact(dis);
  }

  void EnableRedistribution() {fix_redistribution_=false;}
protected:

  // don't want = operator and cctor
  CoNitscheStrategyTsi operator = (const CoNitscheStrategyTsi& old);
  CoNitscheStrategyTsi(const CoNitscheStrategyTsi& old);

  //time integration
  double thermo_alpha_;

  // do not reditribute (during constructor phase)
  bool fix_redistribution_;

  Teuchos::RCP<Epetra_Vector> curr_state_temp_;

  // pointer to TSI coupling object
  Teuchos::RCP<ADAPTER::Coupling> coupST_;


  Teuchos::RCP<Epetra_FEVector> ft_;
  Teuchos::RCP<LINALG::SparseMatrix> ktt_;
  Teuchos::RCP<LINALG::SparseMatrix> ktd_;
  Teuchos::RCP<LINALG::SparseMatrix> kdt_;
};
}
#endif //  CONTACT_NITSCHE_STRATEGY_TSI_H
