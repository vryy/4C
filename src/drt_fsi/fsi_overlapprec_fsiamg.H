#ifdef CCADISCRET

#ifndef FSI_OVERLAPPREC_FSIAMG_H
#define FSI_OVERLAPPREC_FSIAMG_H

#include "fsi_overlapprec.H"
#include "MLAPI_Operator.h"
#include "MLAPI_InverseOperator.h"


namespace FSI
{


  /// special version of block matrix that includes the FSI block preconditioner
  /*!
      This one does coupling of coarse grids for fluid and structure
   */
  class OverlappingBlockMatrixFSIAMG : public OverlappingBlockMatrix
  {
  public:

    /// construction
    OverlappingBlockMatrixFSIAMG(const LINALG::MultiMapExtractor& maps,
                                 ADAPTER::Structure& structure,
                                 ADAPTER::Fluid& fluid,
                                 ADAPTER::Ale& ale,
                                 bool structuresplit,
                                 int symmetric,
                                 vector<double>& omega,
                                 vector<int>& iterations,
                                 double somega,
                                 int siterations,
                                 double fomega,
                                 int fiterations,
                                 FILE* err=NULL);

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    //@}

    /// setup of block preconditioners
    virtual void SetupPreconditioner();

  protected:

    /// symmetric Gauss-Seidel block preconditioner
    virtual void SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// do list for MLAPI smoother
    void SelectMLAPISmoother(std::string& type,
                             const int level,
                             Teuchos::ParameterList& subp,
                             Teuchos::ParameterList& p,
                             Teuchos::ParameterList& pushlist);
    
    /// wrap ILU smoother from ML
    void WrapILUSmoother(ML* ml,
                         MLAPI::Operator& A,
                         MLAPI::InverseOperator& S,const int level);
    
    
    /// generic Vcycle that works on all fields
    virtual void Vcycle(const int level,
                     const int nlevel,
                     MLAPI::MultiVector& z,
                     const MLAPI::MultiVector& b,
                     const vector<MLAPI::Operator>& A,
                     const vector<MLAPI::InverseOperator>& S,
                     const vector<MLAPI::Operator>& P,
                     const vector<MLAPI::Operator>& R,
                     const bool trigger = false) const;
    
    /// FSIAMG multigrid
    virtual void BlockVcycle(const int level,
                     const int nlevel,
                     MLAPI::MultiVector& mlsy,
                     MLAPI::MultiVector& mlfy,
                     MLAPI::MultiVector& mlay,
                     const MLAPI::MultiVector& mlsx,
                     const MLAPI::MultiVector& mlfx,
                     const MLAPI::MultiVector& mlax) const;
    
    /// FSIAMG multigrid with explicit coarse off-diagonals
    virtual void ExplicitBlockVcycle(const int level,
                     const int nlevel,
                     MLAPI::MultiVector& mlsy,
                     MLAPI::MultiVector& mlfy,
                     MLAPI::MultiVector& mlay,
                     const MLAPI::MultiVector& mlsx,
                     const MLAPI::MultiVector& mlfx,
                     const MLAPI::MultiVector& mlax) const;

    /// block Gauss-Seidel smoother within one level
    virtual void BlockGaussSeidelSmoother(
                     const int level,
                     MLAPI::MultiVector& mlsy,
                     MLAPI::MultiVector& mlfy,
                     MLAPI::MultiVector& mlay,
                     const MLAPI::MultiVector& mlsx,
                     const MLAPI::MultiVector& mlfx,
                     const MLAPI::MultiVector& mlax,
                     const bool amgsolve) const;
    
    /// block Gauss-Seidel smoother within one level (explicit off-diagonals)
    virtual void ExplicitBlockGaussSeidelSmoother(
                     const int level,
                     MLAPI::MultiVector& mlsy,
                     MLAPI::MultiVector& mlfy,
                     MLAPI::MultiVector& mlay,
                     const MLAPI::MultiVector& mlsx,
                     const MLAPI::MultiVector& mlfx,
                     const MLAPI::MultiVector& mlax,
                     const bool amgsolve) const;

    /// prolongate to fine grid, multiply and restrict solution back
    virtual MLAPI::MultiVector ProlongateMultiplyRestrict(
                     const int                      level,
                     const MLAPI::MultiVector&      coarse,
                     const vector<MLAPI::Operator>& R,
                     const MLAPI::Operator&         A,
                     const vector<MLAPI::Operator>& P) const;
    
    /// triple matrix product with fine level operator only
    void RAPfine(MLAPI::Operator& RAP,
                 const MLAPI::Operator& R, 
                 Teuchos::RCP<Epetra_CrsMatrix> A,
                 const MLAPI::Operator& P);

    /// triple matrix product with coarse level operators
    void RAPcoarse(MLAPI::Operator& RAP,
                   const MLAPI::Operator& R, 
                   const MLAPI::Operator& A,
                   const MLAPI::Operator& P);

    void RAPoffdiagonals();


    int                             minnlevel_; /// min of the below nlevel_
    vector<int>                     pciter_;
    vector<double>                  pcomega_;

    int                             snlevel_;   /// num level in structure AMG
    Teuchos::ParameterList          sparams_;
    mutable vector<MLAPI::Operator> Ass_;
    vector<MLAPI::InverseOperator>  Sss_;
    vector<MLAPI::Operator>         Pss_;
    vector<MLAPI::Operator>         Rss_;

    int                             fnlevel_;   /// num level in fluid AMG
    Teuchos::ParameterList          fparams_;
    mutable vector<MLAPI::Operator> Aff_;
    vector<MLAPI::InverseOperator>  Sff_;
    vector<MLAPI::Operator>         Pff_;
    vector<MLAPI::Operator>         Rff_;

    int                             anlevel_;   /// num level in ale AMG
    Teuchos::ParameterList          aparams_;
    mutable vector<MLAPI::Operator> Aaa_;
    vector<MLAPI::InverseOperator>  Saa_;
    vector<MLAPI::Operator>         Paa_;
    vector<MLAPI::Operator>         Raa_;

    // these are for implicit off-diagonal in block smoother
    mutable MLAPI::Operator         Asf_;
    mutable MLAPI::Operator         Afs_;
    mutable MLAPI::Operator         Afa_;
    mutable MLAPI::Operator         Aaf_;

    // these are for explicit coarse off-diagonals in block smoother
    mutable vector<MLAPI::Operator> ASF_;
    mutable vector<MLAPI::Operator> AFS_;
    mutable vector<MLAPI::Operator> AFA_;
    mutable vector<MLAPI::Operator> AAF_;

  };


}

#endif
#endif
