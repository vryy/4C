#ifdef CCADISCRET

#ifndef FSI_OVERLAPPREC_FSIAMG_H
#define FSI_OVERLAPPREC_FSIAMG_H

#include "fsi_overlapprec.H"
#include "MLAPI_Operator.h"
#include "MLAPI_LoadBalanceInverseOperator.h"

#define FSIAMG_STRENGTH  0.85   ///< emphasis on strength instead of speed on a 0...1 scale in analysis
#define FSIAMG_ANALYSIS  0      ///< level 0(off),1,2,3,4 switch on live analysis during actual run

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
                                 ALE::Ale& ale,
                                 bool structuresplit,
                                 int symmetric,
                                 vector<double>& omega,
                                 vector<int>& iterations,
                                 vector<double>& somega,
                                 vector<int>& siterations,
                                 vector<double>& fomega,
                                 vector<int>& fiterations,
                                 vector<double>& aomega,
                                 vector<int>& aiterations,
                                 int analyze,
                                 INPAR::FSI::LinearBlockSolver strategy,
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
                         MLAPI::LoadBalanceInverseOperator& S,
                         const int level);


#if 1 // no longer in use
    /// generic Vcycle that works on all fields
    virtual void Vcycle(const int level,
                     const int nlevel,
                     MLAPI::MultiVector& z,
                     const MLAPI::MultiVector& b,
                     const vector<MLAPI::Operator>& A,
                     const vector<RCP<MLAPI::InverseOperator> >& S,
                     const vector<MLAPI::Operator>& P,
                     const vector<MLAPI::Operator>& R,
                     const bool trigger = false) const;
#endif


#if 1 // no longer in use
    /// FSIAMG multigrid with explicit coarse off-diagonals
    virtual void ExplicitBlockVcycle(const int level,
                     const int nlevel,
                     MLAPI::MultiVector& mlsy,
                     MLAPI::MultiVector& mlfy,
                     MLAPI::MultiVector& mlay,
                     const MLAPI::MultiVector& mlsx,
                     const MLAPI::MultiVector& mlfx,
                     const MLAPI::MultiVector& mlax) const;
#endif

#if 1 // no longer in use
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
#endif

#if 1 // no longer in use
    /// iterate on the field individual blocks within the block Gauss Seidel smoother
    virtual void LocalBlockRichardson(
                     const int iterations,
                     const double omega,
                     const int level,
                     const bool amgsolve,
                     const int nlevel,
                     MLAPI::MultiVector& z,
                     const MLAPI::MultiVector& b,
                     const vector<MLAPI::Operator>& A,
                     const vector<RCP<MLAPI::InverseOperator> >& S,
                     const vector<MLAPI::Operator>& P,
                     const vector<MLAPI::Operator>& R
                     ) const;
#endif

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

    // methods for performance analysis of FSIAMG preconditioner
    inline bool Analyze() { return (bool)analyze_; }

    class AnalyzeBest
    {
     public:
       AnalyzeBest(const int nlevel) : 
       nlevel_(nlevel),
       besttype_(6,""),
       bestdamp_(6,1.0),
       bestpoly_(6,1),
       bestsweeps_(6,1),
       S_(7,Teuchos::null)
       {
         if (nlevel-1>6) dserror("Can only analyze V cycles upto 7 levels");
       }
       inline int             Nlevel() { return nlevel_; }
       inline vector<string>& Type()   { return besttype_; }
       inline vector<double>& Damp()   { return bestdamp_; }
       inline vector<int>&    Poly()   { return bestpoly_; }
       inline vector<int>&    Sweeps() { return bestsweeps_; }
       inline vector<RCP<MLAPI::InverseOperator> >& S() { return S_; }
     private:
       int                                  nlevel_;
       vector<string>                       besttype_;
       vector<double>                       bestdamp_;
       vector<int>                          bestpoly_;
       vector<int>                          bestsweeps_;
       vector<RCP<MLAPI::InverseOperator> > S_;
    };
    
    // main routine analysis of FSIAMG
    void AnalyzeFSIAMG(const int myrank,
                       const int snlevel,
                       const int fnlevel,
                       const int anlevel,
                       Teuchos::ParameterList& sparams,
                       Teuchos::ParameterList& fparams,
                       Teuchos::ParameterList& aparams,
                       vector<MLAPI::Operator>& Ass,
                       vector<MLAPI::Operator>& Aff,
                       vector<MLAPI::Operator>& Aaa,
                       vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                       vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                       vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                       vector<MLAPI::Operator>& Asf,
                       vector<MLAPI::Operator>& Afs,
                       vector<MLAPI::Operator>& Afa,
                       vector<MLAPI::Operator>& Aaf,
                       ML* sml,
                       ML* fml,
                       ML* aml);
                       
    // analyze a single field
    void Analyse_SingleField(
                             const string fieldname,
                             const int myrank,
                             string field,
                             const int nlevel,
                             Teuchos::ParameterList& params,
                             vector<MLAPI::Operator>& A,
                             vector<MLAPI::Operator>& P, 
                             vector<MLAPI::Operator>& R,
                             AnalyzeBest& best);
                            
    // analyse BGS(AMG)
    void Analyse_BGSAMG(
                       const int myrank,
                       AnalyzeBest& sbest,
                       AnalyzeBest& fbest,
                       AnalyzeBest& abest,  
                       vector<MLAPI::Operator>& Ass,
                       vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                       vector<MLAPI::Operator>& Aff,
                       vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                       vector<MLAPI::Operator>& Aaa,
                       vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                       vector<MLAPI::Operator>& Asf,
                       vector<MLAPI::Operator>& Afs,
                       vector<MLAPI::Operator>& Afa,
                       vector<MLAPI::Operator>& Aaf
                       );
                       
    // analyse AMG(BGS)
    void Analyse_AMGBGS(
                       const int myrank,
                       AnalyzeBest& sbest,
                       AnalyzeBest& fbest,
                       AnalyzeBest& abest,  
                       vector<MLAPI::Operator>& Ass,
                       vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                       vector<MLAPI::Operator>& Aff,
                       vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                       vector<MLAPI::Operator>& Aaa,
                       vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                       vector<MLAPI::Operator>& Asf,
                       vector<MLAPI::Operator>& Afs,
                       vector<MLAPI::Operator>& Afa,
                       vector<MLAPI::Operator>& Aaf
                       );

    
    double Rate(const int myrank, 
                double t, double r, double initr, double l) const
    {
      //r /= sqrt(l); initr /= sqrt(l);
      double linrate = t*(r/initr);
      double lograte = 0.0;
      if (r/initr>1.0) lograte = 2.0;
      else             lograte = t/(-log(r/initr));
      double rate = FSIAMG_STRENGTH*linrate + (1.-FSIAMG_STRENGTH)*lograte;
#if (FSIAMG_ANALYSIS>=4)
      if (r>initr) 
        if (!myrank) printf("**residual increase in individual field** r0/r:     %8.4e/%8.4e\n",initr,r);
#endif
      return rate;
    }

    
    // a single field single level Richardson iteration with smoother
    double RichardsonS( 
                    const string field,
                    const int myrank,
                    const int level,
                    const int sweeps,
                    const double damp,
                    const MLAPI::Operator& A,
                    const MLAPI::InverseOperator& S,
                    MLAPI::MultiVector& x,
                    const MLAPI::MultiVector& f,
                    bool initiguesszero=false,
                    bool analysis=false,
                    bool silent=true) const;

    // a non-Trilinos preconditioner Richardson
    double RichardsonMixed( 
                    const string field,
                    const int myrank,
                    const int level,
                    const int sweeps,
                    const double damp,
                    const MLAPI::Operator& A,
                    const LINALG::SparseMatrix& matrix,
                    const Teuchos::RCP<LINALG::Preconditioner>& solver,
                    MLAPI::MultiVector& x,
                    const MLAPI::MultiVector& f,
                    int& run,
                    bool initiguesszero=false,
                    bool analysis=false,
                    bool silent=true) const;


    // a single field V cycle analysis Richardson iteration
    double RichardsonV(
                    const string field,
                    const int myrank,
                    int sweeps,
                    const double damp,
                    vector<int>& levelsweeps,
                    vector<double>& leveldamps,
                    vector<MLAPI::Operator>& A,
                    vector<RCP<MLAPI::InverseOperator> >& S,
                    vector<MLAPI::Operator>& P, 
                    vector<MLAPI::Operator>& R,
                    const int level,
                    const int nlevel,
                    MLAPI::MultiVector& x,
                    const MLAPI::MultiVector& f,
                    bool initiguesszero=false,
                    bool analysis=false,
                    bool silent=false) const;
                    
    // a single field V cycle
    void Vcycle(
                const string field,
                const int myrank,
                vector<int>& sweeps,
                vector<double>& damps,
                const int level,
                const int nlevel,
                MLAPI::MultiVector& z,
                const MLAPI::MultiVector& b,
                const vector<MLAPI::Operator>& A,
                const vector<RCP<MLAPI::InverseOperator> >& S,
                const vector<MLAPI::Operator>& P,
                const vector<MLAPI::Operator>& R) const;



    // BGS(AMG) Richardson iteration with Block Gauss Seidel (BGS)
    // and V cycle (V) for individual fields
    double RichardsonBGS_V(
                          const int myrank,
                          const int sweeps,
                          const double damp,
                          vector<int>& blocksweeps,
                          vector<double>& blockdamps,
                          AnalyzeBest& sbest,
                          AnalyzeBest& fbest,
                          AnalyzeBest& abest,
                          MLAPI::MultiVector& sy,
                          MLAPI::MultiVector& fy,
                          MLAPI::MultiVector& ay,
                          const MLAPI::MultiVector& sf,
                          const MLAPI::MultiVector& ff,
                          const MLAPI::MultiVector& af,  
                          vector<MLAPI::Operator>& Ass,
                          vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                          vector<MLAPI::Operator>& Aff,
                          vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                          vector<MLAPI::Operator>& Aaa,
                          vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                          vector<MLAPI::Operator>& Asf,
                          vector<MLAPI::Operator>& Afs,
                          vector<MLAPI::Operator>& Afa,
                          vector<MLAPI::Operator>& Aaf,
                          bool initiguesszero=false,
                          bool analysis=false,
                          bool silent=false
                          ) const;
    
    // BGS(AMG) Richardson iteration with Block Gauss Seidel (BGS)
    // and V cycle (V) or other for individual fields
    double RichardsonBGS_Mixed(
                          const int myrank,
                          const int sweeps,
                          const double damp,
                          vector<int>& blocksweeps,
                          vector<double>& blockdamps,
                          const bool sisamg,
                          const bool fisamg,
                          const bool aisamg,
                          AnalyzeBest& sbest,
                          AnalyzeBest& fbest,
                          AnalyzeBest& abest,
                          MLAPI::MultiVector& sy,
                          MLAPI::MultiVector& fy,
                          MLAPI::MultiVector& ay,
                          const MLAPI::MultiVector& sf,
                          const MLAPI::MultiVector& ff,
                          const MLAPI::MultiVector& af,  
                          vector<MLAPI::Operator>& Ass,
                          vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                          vector<MLAPI::Operator>& Aff,
                          vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                          vector<MLAPI::Operator>& Aaa,
                          vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                          vector<MLAPI::Operator>& Asf,
                          vector<MLAPI::Operator>& Afs,
                          vector<MLAPI::Operator>& Afa,
                          vector<MLAPI::Operator>& Aaf,
                          bool initiguesszero=false,
                          bool analysis=false,
                          bool silent=false
                          ) const;

    // Richardson iteration with Block Gauss Seidel (BGS) 
    // and Smoother (S) or Vcycle (V) for individual fields
    double RichardsonBGS_SV(
                            const int myrank,
                            const int sweeps,
                            const double damp,
                            vector<int>* blocksweeps,
                            vector<double>* blockdamps,
                            const int level,
                            AnalyzeBest& sbest,
                            AnalyzeBest& fbest,
                            AnalyzeBest& abest,
                            MLAPI::MultiVector& sy,
                            MLAPI::MultiVector& fy,
                            MLAPI::MultiVector& ay,
                            const MLAPI::MultiVector& sf,
                            const MLAPI::MultiVector& ff,
                            const MLAPI::MultiVector& af,  
                            vector<MLAPI::Operator>& Ass,
                            vector<MLAPI::Operator>& Pss, vector<MLAPI::Operator>& Rss,
                            vector<MLAPI::Operator>& Aff,
                            vector<MLAPI::Operator>& Pff, vector<MLAPI::Operator>& Rff,
                            vector<MLAPI::Operator>& Aaa,
                            vector<MLAPI::Operator>& Paa, vector<MLAPI::Operator>& Raa,
                            MLAPI::Operator& Asf,
                            MLAPI::Operator& Afs,
                            MLAPI::Operator& Afa,
                            MLAPI::Operator& Aaf,
                            bool initiguesszero=false,
                            bool analysis=false,
                            bool silent=false
                            ) const;
    
    // AMG(BGS) Richardson iteration with block V cycle for FSI
    double Richardson_BlockV(
                             const int myrank,
                             const int sweeps,
                             const double damp, 
                             vector<int>& Vsweeps,
                             vector<double>& Vdamps,
                             vector<int>* blocksweeps,
                             vector<double>* blockdamps,
                             AnalyzeBest& sbest,
                             AnalyzeBest& fbest,
                             AnalyzeBest& abest,
                             MLAPI::MultiVector& sy,
                             MLAPI::MultiVector& fy,
                             MLAPI::MultiVector& ay,
                             const MLAPI::MultiVector& sf,
                             const MLAPI::MultiVector& ff,
                             const MLAPI::MultiVector& af,  
                             vector<MLAPI::Operator>& Ass,
                             vector<MLAPI::Operator>& Pss, 
                             vector<MLAPI::Operator>& Rss,
                             vector<MLAPI::Operator>& Aff,
                             vector<MLAPI::Operator>& Pff, 
                             vector<MLAPI::Operator>& Rff,
                             vector<MLAPI::Operator>& Aaa,
                             vector<MLAPI::Operator>& Paa, 
                             vector<MLAPI::Operator>& Raa,
                             vector<MLAPI::Operator>& Asf,
                             vector<MLAPI::Operator>& Afs,
                             vector<MLAPI::Operator>& Afa,
                             vector<MLAPI::Operator>& Aaf,
                             bool initiguesszero=false,
                             bool analysis=false,
                             bool silent=false) const;

    // AMG(BGS) block V cycle for FSI
    void BlockVcycle(
                     const int myrank,
                     vector<int>& Vsweeps,
                     vector<double>& Vdamps,
                     vector<int>* blocksweeps,
                     vector<double>* blockdamps,
                     const int level,
                     const int minnlevel,
                     AnalyzeBest& sbest,
                     AnalyzeBest& fbest,
                     AnalyzeBest& abest,
                     MLAPI::MultiVector& sy,
                     MLAPI::MultiVector& fy,
                     MLAPI::MultiVector& ay,
                     const MLAPI::MultiVector& sf,
                     const MLAPI::MultiVector& ff,
                     const MLAPI::MultiVector& af,  
                     vector<MLAPI::Operator>& Ass,
                     vector<MLAPI::Operator>& Pss, 
                     vector<MLAPI::Operator>& Rss,
                     vector<MLAPI::Operator>& Aff,
                     vector<MLAPI::Operator>& Pff, 
                     vector<MLAPI::Operator>& Rff,
                     vector<MLAPI::Operator>& Aaa,
                     vector<MLAPI::Operator>& Paa, 
                     vector<MLAPI::Operator>& Raa,
                     vector<MLAPI::Operator>& Asf,
                     vector<MLAPI::Operator>& Afs,
                     vector<MLAPI::Operator>& Afa,
                     vector<MLAPI::Operator>& Aaf) const;
    
    inline const bool& SisAMG() const { return sisml_; }
    inline const bool& FisAMG() const { return fisml_; }
    inline const bool& AisAMG() const { return aisml_; }
    
    bool                                 sisml_;
    bool                                 fisml_;
    bool                                 aisml_;
    int                                  srun_;
    int                                  frun_;
    int                                  arun_;
    
    int                                  minnlevel_; /// min of the below nlevel_
    int                                  maxnlevel_; /// max of the below nlevel_
    int                                  analyze_;   /// run analysis of FSIAMG
    INPAR::FSI::LinearBlockSolver        strategy_; /// type of preocnidtioner  to run: BGS(AMG) or AMG(BGS)
    vector<double>                       pcomega_;
    vector<int>                          pciter_;
    vector<double>                       somega_;
    vector<int>                          siterations_;
    vector<double>                       fomega_;
    vector<int>                          fiterations_;
    vector<double>                       aomega_;
    vector<int>                          aiterations_;

    int                                  snlevel_;   /// num level in structure AMG
    Teuchos::ParameterList               sparams_;
    mutable vector<MLAPI::Operator>      Ass_;
    vector<RCP<MLAPI::InverseOperator> > Sss_;
    vector<MLAPI::Operator>              Pss_;
    vector<MLAPI::Operator>              Rss_;

    int                                  fnlevel_;   /// num level in fluid AMG
    Teuchos::ParameterList               fparams_;
    mutable vector<MLAPI::Operator>      Aff_;
    vector<RCP<MLAPI::InverseOperator> > Sff_;
    vector<MLAPI::Operator>              Pff_;
    vector<MLAPI::Operator>              Rff_;

    int                                  anlevel_;   /// num level in ale AMG
    Teuchos::ParameterList               aparams_;
    mutable vector<MLAPI::Operator>      Aaa_;
    vector<RCP<MLAPI::InverseOperator> > Saa_;
    vector<MLAPI::Operator>              Paa_;
    vector<MLAPI::Operator>              Raa_;

    // these are for implicit off-diagonal in block smoother
    mutable MLAPI::Operator              Asf_;
    mutable MLAPI::Operator              Afs_;
    mutable MLAPI::Operator              Afa_;
    mutable MLAPI::Operator              Aaf_;

    // these are for explicit coarse off-diagonals in block smoother
    mutable vector<MLAPI::Operator>      ASF_;
    mutable vector<MLAPI::Operator>      AFS_;
    mutable vector<MLAPI::Operator>      AFA_;
    mutable vector<MLAPI::Operator>      AAF_;

  };


}

#endif
#endif
