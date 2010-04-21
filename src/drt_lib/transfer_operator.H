/*
 * transfer_operator.H
 *
 *  Created on: Apr 20, 2010
 *      Author: wiesner
 */

#ifdef CCADISCRET

#ifndef TRANSFER_OPERATOR_H_
#define TRANSFER_OPERATOR_H_

// Trilinos includes
#include <Epetra_CrsMatrix.h>
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Time.h"

#include "linalg_sparsematrix.H"

#include "drt_dserror.H"
#include "linalg_mapextractor.H"
#include "linalg_sparsematrix.H"
#include "linalg_blocksparsematrix.H"
#include "linalg_utils.H"
#include "linalg_ana.H"

using namespace std;
using namespace Teuchos;

namespace LINALG
{

class TransferOperator : public virtual Epetra_Operator
{

  public:
    explicit TransferOperator(const RCP<SparseMatrix>& A, FILE* outfile = NULL);
    virtual ~TransferOperator();

    virtual const char* Label() const { return "TransferOperator"; }
    virtual const Epetra_Comm& Comm() const { return A_->Comm(); }

    virtual const Epetra_Map& OperatorDomainMap() const { return A_->OperatorDomainMap(); } // TODO fix me
    virtual const Epetra_Map& OperatorRangeMap()  const { return A_->OperatorRangeMap();  } // TODO fix me

    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      dserror("ApplyInverse not implemented"); return -1;
    }

    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      dserror("Apply not implemented"); return -1;
    }

    virtual int SetUseTranspose(bool UseTranspose)
    {
      dserror("SetUseTranspose not implemented");
      return -1;
    }

    virtual bool UseTranspose() const
    {
      dserror("UseTranspose not implemented");
      return false;
    }

    virtual bool HasNormInf() const
    {
      dserror("HasNormInf not implemented");
      return false;
    }

    virtual double NormInf() const
    {
      dserror("NormInf not implemented");
      return -1.0;
    }

    virtual ostream& Print(std::ostream& os) const;

    /*!
    \brief build transfer operators (prolongator and restrictor) for current level

    This function calculates the prolongator and restrictor from aggregation info.
    We only need the complete aggregation info and a multi vector with a prediction of the nullspace
    for the current level. The domainoffset parameter is optional.
    The function returns a multivector with a prediction for the nullspace of the next coarser level.

    \param aggs (in): aggregation information (rowmap of aggs should be the same as rowmap)
    \param naggs_local (in): number of aggregates for current proc
    \param params (in): parameter list
    \param ThisNS (in): nullspace for this level
    \param domainoffset (in): offset for domain maps gids of next level (only interesting for Braess-Sarazin velocity-pressure splitting)

    \return nullspace for next level

    \author Tobias Wiesner (wiesner@lnm.mw.tum.de)
    */
    virtual RCP<Epetra_MultiVector> buildTransferOperators(const RCP<Epetra_IntVector> aggs, int naggs_local, Teuchos::ParameterList& params, const RCP<Epetra_MultiVector>& ThisNS, const int domainoffset = 0) = 0;

    RCP<SparseMatrix> Prolongator() {return prolongator_; };
    RCP<SparseMatrix> Restrictor() {return restrictor_; };
    const SparseMatrix& P() {return *prolongator_; };
    const SparseMatrix& R() {return *restrictor_; };

  protected:
    const Teuchos::RCP<SparseMatrix>       A_;              ///< Input matrix
    Teuchos::RCP<SparseMatrix>       prolongator_;    ///< prolongator matrix
    Teuchos::RCP<SparseMatrix>       restrictor_;     ///< restrictor matrix
    FILE*                                 outfile_;

};


class TransferOperatorFactory
{
  public:
    TransferOperatorFactory() {};
    ~TransferOperatorFactory() {};

    static RCP<TransferOperator> Create(const string TransferOperatorType, const Teuchos::RCP<SparseMatrix>& A, FILE* outfile = NULL);
};

} // namespace LINALG

#endif /* TRANSFER_OPERATOR_H_ */
#endif /* CCADISCRET */
