/*
 * saddlepointpreconditioner.H
 *
 *  Created on: Feb 16, 2010
 *      Author: wiesner
 */

/*!
  \file saddlepointpreconditoiner.H
  \brief implementation of multigrid saddlepoint preconditioner for fluid problems
 */

#ifndef SADDLEPOINTPRECONDITIONER_H_
#define SADDLEPOINTPRECONDITIONER_H_

#undef USE_MLAPI

// Trilinos includes
#include <Epetra_CrsMatrix.h>
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Time.h"

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_ana.H"

#include "amg/braesssarazin_smoother.H"
#include "amg/transfer_operator.H"

using namespace std;
using namespace Teuchos;

namespace LINALG
{

/*!
 * \brief class SaddlePointPreconditioner: implements saddle point preconditioner for fluid problems
 *
 * this class implements a multigrid saddle point preconditioner for fluid problems (Navier-Stokes equations)
 * We assume, that the number and ratio of velocity and pressure degrees of freedom per node is the same for all nodes.
 *
 * The "AMGBS Parameters" list handles following parameters
 * - "amgbs: smoother: pre or post": select if pre- and/or postsmoothing is done within V-cycle (default: "both", choose from: "both", "pre", "post")
 * - "amgbs: prolongator smoother (vel)": prolongation smoothing/improving strategy for velocity part (default: "PA-AMG", can also be "SA-AMG" and "PG-AMG")
 * - "amgbs: prolongator smoother (pre)": prolongation smoothing/improving strategy for pressure part (default: "PA-AMG", can also be "SA-AMG" and "PG-AMG")
 * - "max coarse dimension": the maximal dimension of coarse grid system before aggregation process stops (default: 20)
 * - "ML output": verbosity level between 0 and 10 (default: 0 = no debug output)
 * - "PDE equations": number of underlaying PDE equations. corresponds to number of degrees of freedom per node (should be 3 for 2D and 4 for 3D problems)
 * - "aggregation: type": algorithm for aggregation process (default: "Uncoupled")
 * - "aggregation: threshold": threshold parameter for aggregation process (default: 0.0)
 * - furthermore all ML parameters are transferred to ML aggregation method
 */
class SaddlePointPreconditioner : public virtual Epetra_Operator
{
  public:
    /*! constructor */
    SaddlePointPreconditioner(RCP<Epetra_Operator> A, const ParameterList& params, FILE* outfile = NULL);
    virtual ~SaddlePointPreconditioner();

    virtual const char* Label() const { return &label_[0]; }
    const Epetra_Comm& Comm() const { return (Ainput_->Comm()); }

    const Epetra_Map& OperatorDomainMap() const { return Ainput_->FullDomainMap(); }
    const Epetra_Map& OperatorRangeMap() const { return Ainput_->FullRangeMap(); }

    /*! entry function for V-Cycle */
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    int SetUseTranspose(bool UseTranspose)
    {
         dserror("SetUseTranspose not impl.");
         return -1;
    }

    double NormInf() const
    {
         dserror("NormInf not impl.");
         return(-1.0);
    }

    bool UseTranspose() const
    {
         dserror("UseTranspose not impl.");
         return false;
    }

    bool HasNormInf() const
    {
         dserror("HasNormInf not impl.");
         return false;
    }

    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      dserror("Apply does not make sense for LINALG::SaddlePointPreconditioner");
      return(-1);
    }

  private:
    // don't want copy-ctor and = operator
    SaddlePointPreconditioner(SaddlePointPreconditioner& old);
    SaddlePointPreconditioner operator = (const SaddlePointPreconditioner& old);

    /*! \brief setup phase of multigrid algorithm
     * performs the multigrid setup
     *
     * -# for each level first the aggregates for the velocity part are generated using ML aggregation
     * -# the aggregation information is transferred to the pressure part. here we assume, that we can find corresponding pressure degrees of freedom for all aggregated velocity degrees of freedom
     * -# next the transfer operators (separately for velocity and pressure) are built.
     * -# with the transfer operators the system matrix for the next level is calculated.
     * -# the BraessSarazin_Smoother objects for pre- and postsmoothing in the V-cycle are generated.
     * -# finally it's checked if the setup process can be stopped
     *
     * \param A: Epetra_Operator: should be a BlockSparseMatrixBase object or a sparse matrix or at least a SparseMatrix.
     * \param origlist: parameter list for saddle-point preconditioner
     *
     * \note if A is a SparseMatrix we have to split the matrix in 2x2 blocks every time this function is called. That's too expensive.
     *
     * \author Tobias Wiesner (wiesner@lnm.mw.tum.de)
     */
    void Setup(RCP<Epetra_Operator>& A, const ParameterList& origlist);

    /*! \brief implementation of the V-Cycle
     * this function is called recursively
     * \author Tobias Wiesner (wiesner@lnm.mw.tum.de)
     */
    int VCycle(const Epetra_MultiVector& Xvel, const Epetra_MultiVector& Xpre, Epetra_MultiVector& Yvel, Epetra_MultiVector& Ypre, const int level) const;

    /*!
    \brief fast Matrix-Matrix-Matrix multiplication

           uses Michaels MLMultiply with ML routines for matrix-matrix multiplication
           returnes ret = A * B * C
    \param A (in): Matrix A
    \param B (in): Matrix B
    \param C (in): Matrix C
    \param bComplete (in): call FillComplete in the end

    \return this = A*B*C

    \author Tobias Wiesner (wiesner@lnm.mw.tum.de)
    */
    RCP<SparseMatrix> Multiply(const SparseMatrix& A, const SparseMatrix& B, const SparseMatrix& C, bool bComplete = true);




    string                                              label_;      ///< name of this operator

    ParameterList                                       params_;    ///< parameter list for SaddlePoint Preconditioner
    FILE*                                               outfile_;   ///< error file to write output to or NULL

    int                                                 nmaxlevels_;  // /<number of maximal allowed multigrid levels
    int                                                 nlevels_;     ///< number of used multigrid levels

    RCP<BlockSparseMatrixBase>                          Ainput_; ///< input matrix

    std::vector<RCP<SparseMatrix> >                    A11_;    ///< vector with splitted blocks of matrix for all multigrid levels
    std::vector<RCP<SparseMatrix> >                    A12_;    ///< vector with splitted blocks of matrix for all multigrid levels
    std::vector<RCP<SparseMatrix> >                    A21_;    ///< vector with splitted blocks of matrix for all multigrid levels
    std::vector<RCP<SparseMatrix> >                    A22_;    ///< vector with splitted blocks of matrix for all multigrid levels

    std::vector<RCP<TransferOperator> >                Tvel_;
    std::vector<RCP<TransferOperator> >                Tpre_;

    std::vector<RCP<BraessSarazin_Smoother> >          preS_; // presmoother (bad: we only support BraessSarazin, as we cannot use the standard ApplyInverse)
    std::vector<RCP<BraessSarazin_Smoother> >          postS_; // postsmoother

    RCP<BraessSarazin_Smoother>                        coarsestSmoother_; // we also use BraessSarazin on coarsest level


    MultiMapExtractor                                  mmex_;    // a  multimapetxractor to handle extracts

    int                                                nVerbose_; // verbosity level (set in Setup)

    bool                                               bPresmoothing_;        ///< presmoothing in V-Cycle
    bool                                               bPostsmoothing_;       ///< postsmoothing in V-Cycle

};

}

#endif /* SADDLEPOINTPRECONDITIONER_H_ */
