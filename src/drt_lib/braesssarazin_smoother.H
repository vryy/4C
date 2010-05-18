/*
 * braesssarazin_smoother.H
 *
 *  Created on: Feb 22, 2010
 *      Author: wiesner
 */

#ifdef CCADISCRET

#ifndef BRAESSSARAZIN_SMOOTHER_H_
#define BRAESSSARAZIN_SMOOTHER_H_

#include "drt_dserror.H"
#include "linalg_mapextractor.H"
#include "linalg_blocksparsematrix.H"
#include "linalg_solver.H"
#include "linalg_utils.H"
#include "linalg_ana.H"


#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace std;
using namespace Teuchos;

namespace LINALG
{

class BraessSarazin_Smoother : public virtual Epetra_Operator
{
  public:
    explicit BraessSarazin_Smoother(RCP<const SparseMatrix> A11, RCP<const SparseMatrix> A12, RCP<const SparseMatrix> A21, RCP<const SparseMatrix> A22, const ParameterList& params);
    virtual ~BraessSarazin_Smoother() {};

    virtual const char* Label() const { return "Braess-Sarazin Smoother"; }

    const Epetra_Map& OperatorDomainMap() const {return *(LINALG::MergeMap(velmap_,premap_));}
    const Epetra_Map& OperatorRangeMap() const {return *(LINALG::MergeMap(velmap_,premap_));}

    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    int ApplyInverse(const Epetra_MultiVector& velrhs, const Epetra_MultiVector& prerhs, Epetra_MultiVector& velsol, Epetra_MultiVector& presol) const;

    const Epetra_Comm& Comm() const { return F_->Comm(); }

    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
    dserror("Apply does not make sense for LINALG::BraessSarazin_Smoother");
    return(-1);
    }

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

  private:
    bool Setup();

    RCP<SparseMatrix> Multiply(const SparseMatrix& A, bool transA, const SparseMatrix& B, bool transB, const SparseMatrix& C, bool transC, bool bComplete = true); ///< matrix-matrix-matrix multiplication
#if 0
    RCP<SparseMatrix>                  diagFinv_;   // inverse diagonal of fluid operator
#else
    RCP<Epetra_Vector>                 diagFinv_;   ///< inverse diagonal of fluid operator (vector)
#endif
    RCP<const SparseMatrix>                  F_;    ///< fluid operator
    RCP<const SparseMatrix>                  G_;    ///< pressure gradient operator
    RCP<const SparseMatrix>                  D_;    ///< divergence operator
    RCP<const SparseMatrix>                  Z_;    ///< pressure stabilization terms or null block
    RCP<SparseMatrix>                        S_;    ///< Schur complement matrix of F_

    ParameterList             params_;              ///< parameter list for Braess-Sarazin smoother

    double                    omega_;               ///< damping factor for Braess Sarazin

    RCP<Epetra_Operator>      Pp_;                  ///< "preconditioner" for "solving" pressure correction equation (ILU, Jacobi, Gauss-Seidel etc)

    const Epetra_Map          velmap_;              ///< map for velocity dofs (primary variables)
    const Epetra_Map          premap_;              ///< map for pressure dofs (Lagrange multipliers)

    int                       nSweeps_;             ///< number of Braess-Sarazin sweeps

};

} // namespace LINALG

#endif /* BRAESSSARAZIN_SMOOTHER_H_ */
#endif /* CCADISCRET */
