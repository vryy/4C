/*
 * braesssarazin_smoother.H
 *
 *  Created on: Feb 22, 2010
 *      Author: wiesner
 */

/*!
  \file braesssarazin_smoother.H
  \brief implementation of Braess-Sarazin smoother
 */


#ifndef BRAESSSARAZIN_SMOOTHER_H_
#define BRAESSSARAZIN_SMOOTHER_H_

#include "../../linalg/linalg_solver.H"
#include "../../linalg/linalg_utils.H"

using Teuchos::ParameterList;

namespace LINALG
{
/*!
 * \brief class BraessSarazin_Smoother: implements Braess-Sarazin smoothing algorithm
 *
 * BraessSarazin_Smoother implements the Braess-Sarazin smoothing algorithm for a multigrid saddlepoint preconditioner
 *
 * the Braess-Sarazn smoother supports the following parameters
 * - "braess-sarazin: sweeps": number of Braess-Sarazin sweeps (default: 1)
 * - "braess-sarazin: damping factor": damping factor for Braess-Sarazin (default: 1.3)
 * - "pressure correction approx: type": type of pressure correction approximation (default: "ILU")
 *  possible parameters for this are "ILU", "Jacobi", "Gauss-Seidel", "symmetric Gauss-Seidel", "Jacobi stand-alone", "Gauss-Seidel stand-alone", "symmetric Gauss-Seidel stand-alone".
 *  These methods use the parameters in the "IFPACK parameters" sublist.
 *  Furthermore you can use "KLU" and "Umfpack" as direct solvers.
 *
 *  For BACI users: The "IFPACK parameters" sublist is filled using the IFPACK specific parameters from the dat file.
 *
 */
class BraessSarazin_Smoother : public virtual Epetra_Operator
{
  public:
    /*!
     *  constructor
     *  \param A11: matrix block (1,1)
     *  \param A12: matrix block (1,2)
     *  \param A21: matrix block (2,1)
     *  \param A22: matrix block (2,2)
     *  \param params: parameter list
     */
    explicit BraessSarazin_Smoother(Teuchos::RCP<const SparseMatrix> A11, Teuchos::RCP<const SparseMatrix> A12, Teuchos::RCP<const SparseMatrix> A21, Teuchos::RCP<const SparseMatrix> A22, const Teuchos::ParameterList& params);
    virtual ~BraessSarazin_Smoother() {};

    virtual const char* Label() const { return "Braess-Sarazin Smoother"; }

    /*! returns the domain map of the overall operator
     *  \return Epetra_Map for the overall operator domain map
     */
    const Epetra_Map& OperatorDomainMap() const {return *(LINALG::MergeMap(velmap_,premap_));}

    /*! returns the range map of the overall operator
     *  \return Epetra_Map for the overall operator range map
     */
    const Epetra_Map& OperatorRangeMap() const {return *(LINALG::MergeMap(velmap_,premap_));}

    /*! not supported */
    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /*! application of Braess-Sarazin smoother to (velsol,presol) vector with given rhs (velrhs, prerhs)
     *  performs nSweeps_ sweeps with Braess-Sarazin smoother.
     *  After exit, velsol and presol contain the smoothed solution vector (velocity and pressure part).
     *  The pressure-correction equation is solved approximately with the Pp_ operator, that is generated in Setup.
     *
     *  \param velrhs: rhs vector (velocity part)
     *  \param prerhs: rhs vector (pressure part)
     *  \param velsol: solution vector (velocity part)
     *  \param presol: solution vector (pressure part)
     *  \return 0
     */
    int ApplyInverse(const Epetra_MultiVector& velrhs, const Epetra_MultiVector& prerhs, Epetra_MultiVector& velsol, Epetra_MultiVector& presol, int level = 0) const;

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

    Teuchos::RCP<SparseMatrix> Multiply(const SparseMatrix& A, bool transA, const SparseMatrix& B, bool transB, const SparseMatrix& C, bool transC, bool bComplete = true); ///< matrix-matrix-matrix multiplication
#if 0
    Teuchos::RCP<SparseMatrix>                  diagFinv_;   // inverse diagonal of fluid operator
#else
    Teuchos::RCP<Epetra_Vector>                 diagFinv_;   ///< inverse diagonal of fluid operator (vector)
#endif
    Teuchos::RCP<const SparseMatrix>                  F_;    ///< fluid operator
    Teuchos::RCP<const SparseMatrix>                  G_;    ///< pressure gradient operator
    Teuchos::RCP<const SparseMatrix>                  D_;    ///< divergence operator
    Teuchos::RCP<const SparseMatrix>                  Z_;    ///< pressure stabilization terms or null block
    Teuchos::RCP<SparseMatrix>                        S_;    ///< Schur complement matrix of F_

    Teuchos::ParameterList             params_;              ///< parameter list for Braess-Sarazin smoother

    double                             omega_;               ///< damping factor for Braess Sarazin

    Teuchos::RCP<Epetra_Operator>      Pp_;                  ///< "preconditioner" for "solving" pressure correction equation (ILU, Jacobi, Gauss-Seidel etc)

    const Epetra_Map                   velmap_;              ///< map for velocity dofs (primary variables)
    const Epetra_Map                   premap_;              ///< map for pressure dofs (Lagrange multipliers)

    int                                nSweeps_;             ///< number of Braess-Sarazin sweeps

};

} // namespace LINALG

#endif /* BRAESSSARAZIN_SMOOTHER_H_ */
