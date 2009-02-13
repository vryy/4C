#ifdef CCADISCRET

#ifndef FSI_LAGRANGEPREC_H
#define FSI_LAGRANGEPREC_H

#include "../drt_lib/linalg_blocksparsematrix.H"
#include "../drt_lib/linalg_precond.H"
#include "../drt_lib/linalg_solver.H"

#include "fsi_overlapprec.H"

namespace FSI
{

  /// special version of block matrix that includes the FSI block preconditioner
  class LagrangianBlockMatrix : public BlockPreconditioningMatrix
  {
  public:

    /// construction
    LagrangianBlockMatrix(const LINALG::MultiMapExtractor& maps,
                          ADAPTER::Structure& structure,
                          ADAPTER::Fluid& fluid,
                          ADAPTER::Ale& ale,
                          int symmetric,
                          double omega=1.0,
                          int iterations=1,
                          double somega=1.0,
                          int siterations=0,
                          double fomega=1.0,
                          int fiterations=0,
                          FILE* err=NULL);

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    //@}

    /// setup of block preconditioners
    void SetupPreconditioner();

  private:

    void SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;
  };

}

#endif
#endif
