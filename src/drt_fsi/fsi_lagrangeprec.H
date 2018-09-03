

#ifndef FSI_LAGRANGEPREC_H
#define FSI_LAGRANGEPREC_H

#include "fsi_overlapprec.H"

// forward declarations

namespace ADAPTER
{
  class FSIStructureWrapper;
  class Fluid;
  class AleFsiWrapper;
}  // namespace ADAPTER

namespace FSI
{
  /// special version of block matrix that includes the FSI block preconditioner
  class LagrangianBlockMatrix : public BlockPreconditioningMatrix
  {
   public:
    /// construction
    LagrangianBlockMatrix(const LINALG::MultiMapExtractor& maps,
        ::ADAPTER::FSIStructureWrapper& structure, ::ADAPTER::Fluid& fluid,
        ::ADAPTER::AleFsiWrapper& ale, int symmetric, double omega = 1.0, int iterations = 1,
        double somega = 1.0, int siterations = 0, double fomega = 1.0, int fiterations = 0,
        double aomega = 1.0, int aiterations = 0, FILE* err = NULL);

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    //@}

    /// setup of block preconditioners
    void SetupPreconditioner();

   private:
    void SGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    ::ADAPTER::FSIStructureWrapper& structure_;
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocks_;
  };

}  // namespace FSI

#endif
