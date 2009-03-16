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
    virtual void SetupPreconditioner();

  protected:

    /// symmetric Gauss-Seidel block preconditioner
    virtual void SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// do list for MLAPI smoother
    void SelectMLAPISmoother(std::string& type,
                             Teuchos::ParameterList& subp,
                             Teuchos::ParameterList& p);
    
    /// generic Vcycle that works on all fields
    void Vcycle(const int level,
                const int nlevel,
                MLAPI::MultiVector& z,
                const MLAPI::MultiVector& b,
                const vector<MLAPI::Operator>& A,
                const vector<MLAPI::InverseOperator>& S,
                const vector<MLAPI::Operator>& P,
                const vector<MLAPI::Operator>& R) const;
    
    int                            minnlevel_; /// min of the below nlevel_

    int                            snlevel_;   /// num level in structure AMG
    Teuchos::ParameterList         sparams_;
    vector<MLAPI::Operator>        Ass_;
    vector<MLAPI::InverseOperator> Sss_;
    vector<MLAPI::Operator>        Pss_;
    vector<MLAPI::Operator>        Rss_;

    int                            fnlevel_;   /// num level in fluid AMG
    Teuchos::ParameterList         fparams_;
    vector<MLAPI::Operator>        Aff_;
    vector<MLAPI::InverseOperator> Sff_;
    vector<MLAPI::Operator>        Pff_;
    vector<MLAPI::Operator>        Rff_;

    int                            anlevel_;   /// num level in ale AMG
    Teuchos::ParameterList         aparams_;
    vector<MLAPI::Operator>        Aaa_;
    vector<MLAPI::InverseOperator> Saa_;
    vector<MLAPI::Operator>        Paa_;
    vector<MLAPI::Operator>        Raa_;

    MLAPI::Operator                Asf_;
    MLAPI::Operator                Afs_;
    MLAPI::Operator                Afa_;
    MLAPI::Operator                Aaf_;
  };


}

#endif
#endif
