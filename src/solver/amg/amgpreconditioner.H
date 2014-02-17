/*
 * amgpreconditioner.H
 *
 *  Created on: Jun 7, 2010
 *      Author: wiesner
 */

#ifndef AMGPRECONDITIONER_H_
#define AMGPRECONDITIONER_H_

// Trilinos includes
#include "Teuchos_ParameterList.hpp"

// BACI includes
#include "../../linalg/linalg_sparsematrix.H"

using namespace Teuchos;

namespace LINALG
{
  class TransferOperator;
  class Smoother;

  /*!
   * AMG Preconditioner
   *
   * \note why new AMG preconditioner class and not reuse AMG_Operator: we don't want MLAPI
   */
  class AMGPreconditioner : public virtual Epetra_Operator
  {
    public:

    explicit AMGPreconditioner(Teuchos::RCP<Epetra_Operator> A, const Teuchos::ParameterList& params, FILE* outfile = NULL);
    virtual ~AMGPreconditioner() {}

    const char* Label() const { return &label_[0]; }

    const Epetra_Comm & Comm() const {return(Ainput_->Comm());}

    ParameterList Params() const {return params_;}

    Teuchos::RCP<SparseMatrix> A() {return Ainput_;}

    const Epetra_Map & OperatorDomainMap() const {return Ainput_->OperatorDomainMap();}

    const Epetra_Map & OperatorRangeMap() const {return Ainput_->OperatorRangeMap();}

    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;

    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      dserror("Apply does not make sense for LINALG::AMGPreconditioner");
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

    //void Vcycle(const Epetra_MultiVector& b_f, Epetra_MultiVector& x_f, const int level) const;
    void Vcycle(const Epetra_MultiVector& rhs, Epetra_MultiVector& sol, const int level) const;

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
    Teuchos::RCP<SparseMatrix> Multiply(const SparseMatrix& A, const SparseMatrix& B, const SparseMatrix& C, bool bComplete = true);

    // don't want copy-ctor and = operator
    AMGPreconditioner(AMGPreconditioner& old);
    AMGPreconditioner operator = (const AMGPreconditioner& old);


    void Setup();

    string                  label_;      // name of this operator
    ParameterList           params_;     // ml parameter list (input)
    Teuchos::RCP<SparseMatrix>       Ainput_;      ///< input matrix
    int                     nlevel_;      ///< no. of grids in hierarchy
    int                     nmaxlevels_;  ///< no. of maximal allowed number of levels
    std::vector<Teuchos::RCP<SparseMatrix> >  A_;          ///< the matrix hierarchy

    std::vector<Teuchos::RCP<Smoother> >          preS_; // presmoother
    std::vector<Teuchos::RCP<Smoother> >          postS_; // postsmoother
    Teuchos::RCP<Epetra_Operator>                   coarsestSmoother_; // we also use BraessSarazin on coarsest level

    std::vector<Teuchos::RCP<TransferOperator> >  T_;          ///< transfer operators

    int nVerbose_;                               ///< verbosity level
  };
}


#endif /* AMGPRECONDITIONER_H_ */
