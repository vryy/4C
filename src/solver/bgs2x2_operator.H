/*!----------------------------------------------------------------------
\file bgs2x2_operator.H

<pre>
\brief Declaration
\level 1
\maintainer Lena Wiechert
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef BGS2X2_OPERATOR_H
#define BGS2X2_OPERATOR_H

// Trilinos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// BACI headers
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_mapextractor.H"       // for mmex_ ( mmex is a reference by intention)
#include "../linalg/linalg_blocksparsematrix.H"  // TODO remove me

namespace LINALG
{
  // forward declaration
  class Preconditioner;

  /// Block Gauss-Seidel preconditioner for a 2x2 system
  class BGS2x2_Operator : public Epetra_Operator
  {
   public:
    /*!
    \brief Standard Constructor
    */
    explicit BGS2x2_Operator(Teuchos::RCP<Epetra_Operator> A, const Teuchos::ParameterList& list1,
        const Teuchos::ParameterList& list2, int global_iter, double global_omega, int block1_iter,
        double block1_omega, int block2_iter, double block2_omega, bool fliporder,
        FILE* outfile = NULL);

    /*!
    \brief Destructor
    */
    virtual ~BGS2x2_Operator()
    {
      solver1_ = Teuchos::null;
      solver2_ = Teuchos::null;
    }

    /// Llabel of this class.
    const char* Label() const { return "LINALG::BGS2x2_Operator"; }

    /// Comm of this class
    const Epetra_Comm& Comm() const { return (A_->Comm()); }


    /// Operator domain map
    const Epetra_Map& OperatorDomainMap() const { return A_->FullDomainMap(); }

    /// Operator range map
    const Epetra_Map& OperatorRangeMap() const { return A_->FullRangeMap(); }

    /// Setup of preconditioners for individual blocks
    void SetupBlockPreconditioners();

    /// Apply inverse of the preconditioner
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /// not implemented
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      dserror("Apply does not make sense for LINALG::BGS2x2_Operator");
      return (-1);
    }

    /// not implemented
    int SetUseTranspose(bool UseTranspose)
    {
      dserror("SetUseTranspose not impl.");
      return -1;
    }

    /// not implemented
    double NormInf() const
    {
      dserror("NormInf not impl.");
      return (-1.0);
    }


    /// not implemented
    bool UseTranspose() const
    {
      dserror("UseTranspose not impl.");
      return false;
    }

    /// not implemented
    bool HasNormInf() const
    {
      dserror("HasNormInf not impl.");
      return false;
    }

   private:
    // don't want copy-ctor and = operator
    BGS2x2_Operator(BGS2x2_Operator& old);
    BGS2x2_Operator operator=(const BGS2x2_Operator& old);

    /// Richardson iteration on one block using the given flags
    void LocalBlockRichardson(Teuchos::RCP<Preconditioner> solver, const SparseMatrix& Op,
        Teuchos::RCP<Epetra_MultiVector> x, Teuchos::RCP<Epetra_MultiVector> y,
        Teuchos::RCP<Epetra_MultiVector> tmpx, int iter, double omega) const;


    FILE* outfile_;                 // error file to write output to or NULL
    Teuchos::ParameterList list1_;  // list for solver of first diagonal block
    Teuchos::ParameterList list2_;  // list for solver of second diagonal block

    MultiMapExtractor mmex_;                 // a  multimapetxractor to handle extracts
    Teuchos::RCP<BlockSparseMatrixBase> A_;  // 2x2 block matrix

    Teuchos::RCP<Preconditioner> solver1_;  // solver of block 1
    Teuchos::RCP<Preconditioner> solver2_;  // solver of block 2

    int global_iter_;
    double global_omega_;
    int block1_iter_;
    double block1_omega_;
    int block2_iter_;
    double block2_omega_;

    int firstind_;  /// index of block "0" in Gauss-Seidel procedure
    int secind_;    /// index of block "1" in Gauss-Seidel procedure

  };  // class  BGS2x2_Operator : public virtual Epetra_Operator
}  // namespace LINALG


#endif  // BGS2x2_OPERATOR_H
