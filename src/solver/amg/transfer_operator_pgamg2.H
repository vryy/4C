/*
 * transfer_operator_pgamg2.H
 *
 *  Created on: Apr 23, 2010
 *      Author: wiesner
 */


#ifndef TRANSFER_OPERATOR_PGAMG2_H_
#define TRANSFER_OPERATOR_PGAMG2_H_


#include "transfer_operator_saamg.H"

using namespace Teuchos;

namespace LINALG
{
  class PGAMG2TransferOperator : public SAAMGTransferOperator
  {
    public:
      explicit PGAMG2TransferOperator(const RCP<SparseMatrix>& A, FILE* outfile);
      virtual ~PGAMG2TransferOperator() {}

      virtual const char* Label() const { return "PGAMG2TransferOperator"; }

      /*!
      \brief build transfer operators (prolongator and restrictor) for current level

      This function calculates the PG-AMG prolongator and restrictor from aggregation info.
      We only need the complete aggregation info and a multi vector with a prediction of the nullspace
      for the current level. The domainoffset parameter is optional.
      The function returns a multivector with a prediction for the nullspace of the next coarser level.

      \param aggs (in): aggregation information
      \param naggs_local (in): number of aggregates for current proc
      \param params (in): parameter list
      \param ThisNS (in): nullspace for this level
      \param domainoffset (in): offset for domain maps gids of next level (only interesting for Braess-Sarazin velocity-pressure splitting)

      \return nullspace for next level

      \author Tobias Wiesner (wiesner@lnm.mw.tum.de)
      */
      RCP<Epetra_MultiVector> buildTransferOperators(const RCP<Epetra_IntVector> aggs, int naggs_local, Teuchos::ParameterList& params, const RCP<Epetra_MultiVector>& ThisNS, const int domainoffset = 0);
    protected:

    private:
      RCP<Epetra_Vector> CalculateColBasedOmegaP(const RCP<SparseMatrix>& A,const RCP<Epetra_Vector>& Dinv,const RCP<SparseMatrix>& Ptent,const RCP<SparseMatrix>& AP0);
      RCP<Epetra_Vector> CalculateRowBasedOmegaR(const RCP<SparseMatrix>& A,const RCP<Epetra_Vector>& Dinv,const RCP<SparseMatrix>& Ptent,const RCP<SparseMatrix>& AtP0);
      RCP<Epetra_Vector> TransformColBased2RowBasedOmegas(const RCP<SparseMatrix>& AP, const RCP<Epetra_Vector> colbasedomegas, double fallbackvalue = 0.0);

      void MultiplyAll(const RCP<SparseMatrix>& left, const RCP<SparseMatrix>& right, RCP<Epetra_Vector>& InnerProd);
      void MultiplySelfAll(const RCP<SparseMatrix>& Op,RCP<Epetra_Vector>& Column2Norm);


      int nVerbose_;    ///< verbosity level
  };

} /* namespace LINALG */


#endif /* TRANSFER_OPERATOR_PGAMG2_H_ */

