/*
 * transfer_operator_pgamg.H
 *
 *  Created on: Apr 20, 2010
 *      Author: wiesner
 */

/*!
  \file transfer_operator_pgamg.H
  \brief implementation for a PG-AMG variant of transfer operators
 */
#ifndef TRANSFER_OPERATOR_PGAMG_H_
#define TRANSFER_OPERATOR_PGAMG_H_

#include "transfer_operator_saamg.H"

using namespace Teuchos;

namespace LINALG
{
  /*! \brief variant of PG-AMG
   *  this transfer operators is dedicated to nonsymmetric (and indefinite problems).
   *  \note not optimized for parallel use
   */
  class PGAMGTransferOperator : public SAAMGTransferOperator
  {
    public:
      explicit PGAMGTransferOperator(const Teuchos::RCP<SparseMatrix>& A, FILE* outfile);
      virtual ~PGAMGTransferOperator() {}

      virtual const char* Label() const { return "PGAMGTransferOperator"; }

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
      Teuchos::RCP<Epetra_MultiVector> buildTransferOperators(const Teuchos::RCP<Epetra_IntVector> aggs, int naggs_local, Teuchos::ParameterList& params, const Teuchos::RCP<Epetra_MultiVector>& ThisNS, const int domainoffset = 0);
    protected:

    private:
      /*!
      \brief internal routine for PG-AMG smoothing

      this function calculates prolongator and restrictor by minimizing the "energy"
      with respect to the (D^{-1} A)' D^{-1} A norm.

                   diag( R0 (A D^{-1}' D^{-1} A' D^{-1} A' R0' )
      omega = ---------------------------------------------------------
                diag( R0 A D^{-1} A D^{-1} D^{-1} A' D^{-1} A' R0' )

      there's no compression for these column-based omegas (immense communication in parallel case)
      if there's a problem in finding some row-based omegas we fall back to SA-AMG omega for these entries.

      \param A (in): matrix for smoothing
      \param P_tent (in): tentative prolongator
      \param R_tent (in): tentative restrictor
      \param P_smoothed (out): smoothed prolongator
      \param R_smoothed (out): smoothed restrictor

      \return void

      \author Tobias Wiesner (wiesner@lnm.mw.tum.de)

       */
      void PG_AMG(const Teuchos::RCP<SparseMatrix>& A, const Teuchos::RCP<SparseMatrix>& P_tent, const Teuchos::RCP<SparseMatrix>& R_tent, Teuchos::RCP<SparseMatrix>& P_smoothed, Teuchos::RCP<SparseMatrix>& R_smoothed);

      void MultiplyAll(const Teuchos::RCP<SparseMatrix>& left, const Teuchos::RCP<SparseMatrix>& right, Teuchos::RCP<Epetra_Vector>& InnerProd);
      void MultiplySelfAll(const Teuchos::RCP<SparseMatrix>& Op,Teuchos::RCP<Epetra_Vector>& Column2Norm);


      int nVerbose_;    ///< verbosity level
  };

} /* namespace LINALG */


#endif /* TRANSFER_OPERATOR_PGAMG_H_ */
