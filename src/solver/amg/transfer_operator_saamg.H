/*
 * transfer_operator_saamg.H
 *
 *  Created on: Apr 20, 2010
 *      Author: wiesner
 */

/*!
  \file transfer_operator_saamg.H
  \brief implementation for SA-AMG transfer operators
 */


#ifndef TRANSFER_OPERATOR_SAAMG_H_
#define TRANSFER_OPERATOR_SAAMG_H_

#include "transfer_operator_tentative.H"

using namespace Teuchos;

namespace LINALG
{
  /*! concrete implementation for the smoothed aggregation transfer operators
   * uses the tentative transfer operators and smoothes them using the SA-AMG method. This smoothing method is dedicated to symmetric positive definite problems only.
   */
  class SAAMGTransferOperator : public TentativeTransferOperator
  {
    public:
      SAAMGTransferOperator(const RCP<SparseMatrix>& A, FILE* outfile);
      virtual ~SAAMGTransferOperator() {}

      virtual const char* Label() const { return "SAAMGTransferOperator"; }

      /*!
      \brief build transfer operators (prolongator and restrictor) for current level

      This function calculates the smoothed aggregation (SAAMG) prolongator and restrictor from aggregation info.
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

     /*!
      \brief calculate max eigenvalue of A using a MLAPI call

      function adapted from MLAPI_Eig.cpp.
      calculates maximum eigenvalue of A using a Krylov CG method.
      if DiagonalScaling == false -> calc max EW of A
      if DiagonalScaling == true  -> calc max EW of D^{-1}A where D is the diagonal of A

      \param A (in): sparse matrix input
      \param DiagonalScaling (in): boolean parameter (see above)

      \return max eigenvalue

      \author Tobias Wiesner (wiesner@lnm.mw.tum.de)
      */
      double MaxEigCG(const SparseMatrix& A, const bool DiagonalScaling);

    protected:


  };

} /* namespace LINALG */


#endif /* TRANSFER_OPERATOR_SAAMG_H_ */ 
