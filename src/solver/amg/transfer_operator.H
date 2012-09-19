/*
 * transfer_operator.H
 *
 *  Created on: Apr 20, 2010
 *      Author: wiesner
 */

/*!
  \file transfer_operator.H
  \brief basis class for multigrid transfer operators
 */


#ifndef TRANSFER_OPERATOR_H_
#define TRANSFER_OPERATOR_H_

// Trilinos includes
#include <Epetra_CrsMatrix.h>
//#include "Teuchos_RefCountPtr.hpp"
//#include "Teuchos_ParameterList.hpp"
#include "Epetra_Time.h"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_ana.H"

using namespace std;
using namespace Teuchos;

namespace LINALG
{

class MapExtractor;
class BlockSparseMatrixBase;
class SparseMatrix;


/*!
 *  \brief basis class for multigrid transfer operators
 *
 *  class provides a virtual interface for concrete transfer operator implementations.
 *  a new transfer operator has to implement the buildTransferOperators function
 */
class TransferOperator : public virtual Epetra_Operator
{

  public:
    explicit TransferOperator(const RCP<SparseMatrix>& A, FILE* outfile = NULL);
    virtual ~TransferOperator();

    virtual const char* Label() const { return "TransferOperator"; }
    virtual const Epetra_Comm& Comm() const { return A_->Comm(); }

    virtual const Epetra_Map& OperatorDomainMap() const { return prolongator_->OperatorDomainMap(); } // TODO fix me
    virtual const Epetra_Map& OperatorRangeMap()  const { return prolongator_->OperatorRangeMap();  } // TODO fix me

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

    RCP<SparseMatrix> Prolongator() {return prolongator_; };///< returns the prolongation operator
    RCP<SparseMatrix> Restrictor() {return restrictor_; }; ///< returns the restriction operator
    const SparseMatrix& P() {return *prolongator_; }; ///< returns the prolongation operator
    const SparseMatrix& R() {return *restrictor_; };  ///< returns the restriction operator

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
