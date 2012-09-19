/*
 * transfer_operator_tentative.H
 *
 *  Created on: Apr 20, 2010
 *      Author: wiesner
 */

/*!
  \file transfer_operator_tentative.H
  \brief implementation for tentative transfer operators
 */

#ifndef TRANSFER_OPERATOR_TENTATIVE_H_
#define TRANSFER_OPERATOR_TENTATIVE_H_

#include "transfer_operator.H"



using namespace std;
using namespace Teuchos;

namespace LINALG
{

  /*! \brief concrete implementation of the tentative transfer opeartors
   * uses the aggregation information for building the tentative prolongator and restrictor
   */
  class TentativeTransferOperator : public TransferOperator
  {
    public:
      explicit TentativeTransferOperator(const RCP<SparseMatrix>& A, FILE* outfile);
      virtual ~TentativeTransferOperator() {}

      virtual const char* Label() const { return "TentativeTransferOperator"; }

      /*!
      \brief build transfer operators (prolongator and restrictor) for current level

      This function calculates the tentative prolongator and restrictor from aggregation info.
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
      RCP<Epetra_MultiVector> buildTransferOperators(const RCP<Epetra_IntVector> aggs, int naggs_local, Teuchos::ParameterList& params, const RCP<Epetra_MultiVector>& ThisNS, const int domainoffset = 0);
    protected:

      void GetPtent(const Epetra_Map& rowmap, const Epetra_IntVector aggvec, int naggs, ParameterList& List, const Epetra_MultiVector& ThisNS, RCP<SparseMatrix>& Ptent, RCP<Epetra_MultiVector>& NextNS, const int domainoffset);


  };

} /* namespace LINALG */

#endif /* TRANSFER_OPERATOR_TENTATIVE_H_ */
