/*---------------------------------------------------------------------*/
/*!
\file contact_nitsche_strategy_poro.H

\brief Nitsche poro contact solving strategy

\level 3

\maintainer Christoph Ager

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_NITSCHE_STRATEGY_PORO_H
#define CONTACT_NITSCHE_STRATEGY_PORO_H

#include "contact_nitsche_strategy.H"

namespace ADAPTER {class Coupling;}

namespace CONTACT
{

/*!
 \brief Contact solving strategy with Nitsche's method.

 This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
 For a more general documentation of the involved functions refer to CoAbstract Strategy.

 */
class CoNitscheStrategyPoro : public CoNitscheStrategy
{
public:

  //! Standard constructor
  CoNitscheStrategyPoro(
      const Epetra_Map* DofRowMap,
      const Epetra_Map* NodeRowMap,
      Teuchos::ParameterList params,
      std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
      int dim,
      Teuchos::RCP<Epetra_Comm> comm,
      double alphaf,
      int maxdof)
: CoNitscheStrategy(DofRowMap,NodeRowMap,params,interface,dim,comm,alphaf,maxdof)
,no_penetration_(params.get<bool>("CONTACTNOPEN"))
,dd_dv_(params.get<double>("porotimefac"))
{}

  //! Shared data constructor
  CoNitscheStrategyPoro(
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
  ,no_penetration_(params.get<bool>("CONTACTNOPEN"))
  ,dd_dv_(params.get<double>("porotimefac"))
  {}

  virtual void ApplyForceStiffCmt(
      Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<LINALG::SparseOperator>& kt,
      Teuchos::RCP<Epetra_Vector>& f, const int step, const int iter,
      bool predictor);

//  void Integrate(CONTACT::ParamsInterface& cparams);

  /*!
\brief Destructor

   */
  virtual ~CoNitscheStrategyPoro() {};

  virtual void SetState(const enum MORTAR::StateType& statename,
        const Epetra_Vector& vec);

  virtual void SetParentState(const enum MORTAR::StateType& statename,
      const Epetra_Vector& vec);

  virtual Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
     const enum DRT::UTILS::VecBlockType& bt) const;

  virtual Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
      const enum DRT::UTILS::MatBlockType& bt) const ;


  //Flag for Poro No Penetration Condition
  bool HasPoroNoPenetration() const {return no_penetration_;}

protected:

  // don't want = operator and cctor
  CoNitscheStrategyPoro operator = (const CoNitscheStrategyPoro& old);
  CoNitscheStrategyPoro(const CoNitscheStrategyPoro& old);

  // create an appropriate vector for the RHS
  virtual Teuchos::RCP<Epetra_FEVector> SetupRhsBlockVec(
       const enum DRT::UTILS::VecBlockType& bt)const;

  // create an appropriate matrix block
  virtual Teuchos::RCP<LINALG::SparseMatrix> SetupMatrixBlockPtr(
      const enum DRT::UTILS::MatBlockType& bt);

  // complete matrix block with correct maps
  virtual void CompleteMatrixBlockPtr(
      const enum DRT::UTILS::MatBlockType& bt,
      Teuchos::RCP<LINALG::SparseMatrix> kc);

  bool no_penetration_;
  double dd_dv_;


  Teuchos::RCP<Epetra_FEVector> fp_;
  Teuchos::RCP<LINALG::SparseMatrix> kpp_;
  Teuchos::RCP<LINALG::SparseMatrix> kpd_;
  Teuchos::RCP<LINALG::SparseMatrix> kdp_;
};
}
#endif //  CONTACT_NITSCHE_STRATEGY_PORO_H
