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
    Teuchos::RCP<Epetra_Comm> comm,
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


/*!
\brief Destructor

*/
virtual ~CoNitscheStrategy() {};


virtual void SetParentState(const std::string& statename,
    const Teuchos::RCP<Epetra_Vector> vec,
    const Teuchos::RCP<DRT::Discretization> dis);

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

};
}
#endif //  CONTACT_NITSCHE_STRATEGY_H
