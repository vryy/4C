/*
 * saddlepointpreconditioner.H
 *
 *  Created on: Feb 16, 2010
 *      Author: wiesner
 */



#ifdef CCADISCRET
#ifndef SADDLEPOINTPRECONDITIONER_H_
#define SADDLEPOINTPRECONDITIONER_H_

#undef USE_MLAPI

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

#include "MLAPI_Error.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI_Eig.h"

#include "braesssarazin_smoother.H"

using namespace std;
using namespace Teuchos;

namespace LINALG
{

class SaddlePointPreconditioner : public virtual Epetra_Operator
{
  public:
    SaddlePointPreconditioner(RCP<Epetra_Operator> A, const ParameterList& params, const ParameterList& pressurelist, FILE* outfile = NULL);
    virtual ~SaddlePointPreconditioner();

    virtual const char* Label() const { return "SaddlePoint Preconditioner"; }
    const Epetra_Comm& Comm() const { return (Ainput_->Comm()); }

    const Epetra_Map& OperatorDomainMap() const { return Ainput_->FullDomainMap(); }
    const Epetra_Map& OperatorRangeMap() const { return Ainput_->FullRangeMap(); }

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

    void Setup(RCP<Epetra_Operator> A, const ParameterList& origlist, const ParameterList& origplist);
    int VCycle(const Epetra_MultiVector& Xvel, const Epetra_MultiVector& Xpre, Epetra_MultiVector& Yvel, Epetra_MultiVector& Ypre, const int level) const;

    /*!
    \brief Call ML aggregation on A according to parameters supplied in List. Return
           aggregates in aggrinfo

           On input, map of aggrinfo has to map row map of A. On output, aggrinfo[i]
           contains number of aggregate the row belongs to, where aggregates are
           numbered starting from 0.
           Return value is the processor-local number of aggregates build.
           If aggrinfo[i] >= return-value, then i is a processor local row
           of a row that ML has detected to be on a Dirichlet BC.

    \param A (in): Matrix to be aggregated on
    \param List (in): ParameterList containing ML options
    \param ThisNS (in): MultiVector with null space of current level
    \param aggrinfo(out): vector containing aggregation information

    \note Map of aggrinfo has to match rowmap of A on input.

    \return returns processor-local number of aggregates

    \author Michael Gee (gee@lnm.mw.tum.de)
    */
    int GetAggregates(SparseMatrix& A, ParameterList& List, const Epetra_MultiVector& ThisNS, Epetra_IntVector& aggrinfo);

    /*!
    \brief Call ML aggregation on A according to parameters supplied in List. Return
           aggregates in aggrinfo

           On input, map of aggrinfo has to map row map of A. On output, aggrinfo[i]
           contains number of global aggregate the row belongs to, where aggregates are
           numbered starting from 0 globally.
           Return value is the processor-local number of aggregates build.
           If aggrinfo[i] < 0, then i is a processor local row
           that ML has detected to be on a Dirichlet BC.
           if aggrinfo[i] >= 0, then i is a processor local row and aggrinfo[i] is
           a global aggregate id.

    \param A (in): Matrix to be aggregated on
    \param List (in): ParameterList containing ML options
    \param ThisNS (in): MultiVector with null space of current level
    \param aggrinfo(out): vector containing aggregation information in global numbering Epetra_IntVector with aggregate info (maps dof GID -> aggregate id)
    \param naggregates_local(out): number of local aggregates (for current proc)

    \note Map of aggrinfo has to match rowmap of A on input.

    \return global number of aggregates (over all procs)

    \author Michael Gee (gee@lnm.mw.tum.de)
    \author Tobias Wiesner (wiesner@lnm.mw.tum.de)
    */
    int GetGlobalAggregates(SparseMatrix& A, ParameterList& List, const Epetra_MultiVector& ThisNS, Epetra_IntVector& aggrinfo, int& naggregates_local);

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
    void GetPtent(const Epetra_Map& rowmap, const Epetra_IntVector aggvec, int naggs, ParameterList& List, const Epetra_MultiVector& ThisNS, RCP<SparseMatrix>& Ptent, RCP<Epetra_MultiVector>& NextNS, const int domainoffset);

    double MaxEigCG(const SparseMatrix& A, const bool DiagonalScaling);
    void SA_AMG(const RCP<SparseMatrix>& A, const RCP<SparseMatrix>& P_tent, const RCP<SparseMatrix>& R_tent, RCP<SparseMatrix>& P_smoothed, RCP<SparseMatrix>& R_smoothed);

    void GetPtent(const Epetra_Map& rowmap, const Epetra_IntVector aggvec, int naggs, ParameterList& List, const Epetra_MultiVector& ThisNS, RCP<Epetra_CrsMatrix>& Ptent, RCP<Epetra_MultiVector>& NextNS, const int domainoffset);


    ParameterList                                       params_;   // list for velocity solver
    ParameterList                                       pressureparams_;   // list for pressure solver
    FILE*                                               outfile_; // error file to write output to or NULL

    int                                                 nmaxlevels_;  // number of maximal allowed multigrid levels
    int                                                 nlevels_;     // number of used multigrid levels

    RCP<BlockSparseMatrixBase>                          Ainput_; // input matrix

    std::vector<RCP<SparseMatrix> >                    A11_;    // vector with splitted blocks of matrix
    std::vector<RCP<SparseMatrix> >                    A12_;
    std::vector<RCP<SparseMatrix> >                    A21_;
    std::vector<RCP<SparseMatrix> >                    A22_;
    std::vector<RCP<SparseMatrix> >                    Pvel_;
    std::vector<RCP<SparseMatrix> >                    Ppre_;
    std::vector<RCP<SparseMatrix> >                    Rvel_;
    std::vector<RCP<SparseMatrix> >                    Rpre_;

    std::vector<RCP<BraessSarazin_Smoother> >          preS_; // presmoother (bad: we only support BraessSarazin, as we cannot use the standard ApplyInverse)
    std::vector<RCP<BraessSarazin_Smoother> >          postS_; // postsmoother

    RCP<BraessSarazin_Smoother>                        coarsestSmoother_; // we also use BraessSarazin on coarsest level


    MultiMapExtractor                                  mmex_;    // a  multimapetxractor to handle extracts

    bool                                               bPresmoothing_;        // presmoothing in V-Cycle
    bool                                               bPostsmoothing_;       // postsmoothing in V-Cycle

};

}

#endif /* SADDLEPOINTPRECONDITIONER_H_ */
#endif /* CCADISCRET */
