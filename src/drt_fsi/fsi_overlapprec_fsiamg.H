/*----------------------------------------------------------------------*/
/*! \file

\brief Special version of block matrix that includes the FSI block preconditioner


\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FSI_OVERLAPPREC_FSIAMG_H
#define FSI_OVERLAPPREC_FSIAMG_H

#include "fsi_overlapprec.H"
#include "MLAPI_Operator.h"
#include "MLAPI_LoadBalanceInverseOperator.h"

#define FSIAMG_STRENGTH \
  0.85                     ///< emphasis on strength instead of speed on a 0...1 scale in analysis
#define FSIAMG_ANALYSIS 0  ///< level 0(off),1,2,3,4 switch on live analysis during actual run

namespace FSI
{
  class OverlappingBlockMatrixHybridSchwarz;
}

namespace FSI
{
  /*! \brief Special version of block matrix that includes the FSI block preconditioner
   *
   *  This one does coupling of coarse grids for fluid and structure
   */
  class OverlappingBlockMatrixFSIAMG : public OverlappingBlockMatrix
  {
   public:
    /// construction
    OverlappingBlockMatrixFSIAMG(const LINALG::MultiMapExtractor& maps,
        ADAPTER::FSIStructureWrapper& structure, ADAPTER::Fluid& fluid, ADAPTER::AleFsiWrapper& ale,
        bool structuresplit, int symmetric, std::vector<std::string>& blocksmoother,
        std::vector<double>& schuromega, std::vector<double>& omega, std::vector<int>& iterations,
        std::vector<double>& somega, std::vector<int>& siterations, std::vector<double>& fomega,
        std::vector<int>& fiterations, std::vector<double>& aomega, std::vector<int>& aiterations,
        int analyze, INPAR::FSI::LinearBlockSolver strategy, INPAR::FSI::Verbosity verbosity,
        FILE* err = NULL, OverlappingBlockMatrixHybridSchwarz* hybridPrec = NULL);

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    //@}

    /// setup of block preconditioners
    virtual void SetupPreconditioner();

   protected:
    /// symmetric Gauss-Seidel block preconditioner
    virtual void SGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /// do list for MLAPI smoother
    void SelectMLAPISmoother(std::string& type, const int level, Teuchos::ParameterList& subp,
        Teuchos::ParameterList& p, Teuchos::ParameterList& pushlist);

    /// wrap ILU smoother from ML
    void WrapILUSmoother(
        ML* ml, MLAPI::Operator& A, MLAPI::LoadBalanceInverseOperator& S, const int level);


    /// generic Vcycle that works on all fields
    virtual void Vcycle(const int level, const int nlevel, MLAPI::MultiVector& z,
        const MLAPI::MultiVector& b, const std::vector<MLAPI::Operator>& A,
        const std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S,
        const std::vector<MLAPI::Operator>& P, const std::vector<MLAPI::Operator>& R,
        const bool trigger = false) const;


    /// FSIAMG multigrid with explicit coarse off-diagonals
    virtual void ExplicitBlockVcycle(const int level, const int nlevel, MLAPI::MultiVector& mlsy,
        MLAPI::MultiVector& mlfy, MLAPI::MultiVector& mlay, const MLAPI::MultiVector& mlsx,
        const MLAPI::MultiVector& mlfx, const MLAPI::MultiVector& mlax) const;


    /// block Gauss-Seidel smoother within one level (explicit off-diagonals)
    virtual void ExplicitBlockGaussSeidelSmoother(const int level, MLAPI::MultiVector& mlsy,
        MLAPI::MultiVector& mlfy, MLAPI::MultiVector& mlay, const MLAPI::MultiVector& mlsx,
        const MLAPI::MultiVector& mlfx, const MLAPI::MultiVector& mlax, const bool amgsolve) const;


    /// iterate on the field individual blocks within the block Gauss Seidel smoother
    virtual void LocalBlockRichardson(const int iterations, const double omega, const int level,
        const bool amgsolve, const int nlevel, MLAPI::MultiVector& z, const MLAPI::MultiVector& b,
        const std::vector<MLAPI::Operator>& A,
        const std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S,
        const std::vector<MLAPI::Operator>& P, const std::vector<MLAPI::Operator>& R) const;


    /// triple matrix product with fine level operator only
    void RAPfine(MLAPI::Operator& RAP, const MLAPI::Operator& R, Teuchos::RCP<Epetra_CrsMatrix> A,
        const MLAPI::Operator& P);

    /// triple matrix product with coarse level operators
    void RAPcoarse(MLAPI::Operator& RAP, const MLAPI::Operator& R, const MLAPI::Operator& A,
        const MLAPI::Operator& P);

    /// build off-diagonal coupling blocks for FSIAMG
    void RAPoffdiagonals();

    /// build Schur Complement operator from fluid block
    void SchurComplementOperator(MLAPI::Operator& Schur, MLAPI::Operator& Ass, MLAPI::Operator& Aff,
        MLAPI::Operator& Aaa, MLAPI::Operator& Asf, MLAPI::Operator& Afs, MLAPI::Operator& Afa,
        MLAPI::Operator& Aaf, const double omega, const bool structuresplit);


    //! methods for performance analysis of FSIAMG preconditioner
    inline bool Analyze() { return (bool)analyze_; }

    //! Operator to analyze multigrid settings of a single field
    class AnalyzeBest
    {
     public:
      //! constructor
      AnalyzeBest(const int nlevel)
          : nlevel_(nlevel),
            besttype_(6, ""),
            bestdamp_(6, 1.0),
            bestpoly_(6, 1),
            bestsweeps_(6, 1),
            S_(7, Teuchos::null)
      {
        if (nlevel - 1 > 6) dserror("Can only analyze V cycles upto 7 levels");
        return;
      }

      //! Return number of levels
      inline int Nlevel() { return nlevel_; }

      //! Return set of optimal algorithm types
      inline std::vector<std::string>& Type() { return besttype_; }

      //! Return set of optimal damping paramters
      inline std::vector<double>& Damp() { return bestdamp_; }

      //! Return set of optimal polynomial degrees
      inline std::vector<int>& Poly() { return bestpoly_; }

      //! Return set of optimal number of sweeps
      inline std::vector<int>& Sweeps() { return bestsweeps_; }

      //! Access the underlying algorithmic operator
      inline std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S() { return S_; }

     private:
      int nlevel_;                         ///< number of levels
      std::vector<std::string> besttype_;  //< best set of algorithms
      std::vector<double> bestdamp_;       ///< best set of damping parameters
      std::vector<int> bestpoly_;          ///< best set of polynomial degrees
      std::vector<int> bestsweeps_;        ///< best set of number of sweeps
      std::vector<Teuchos::RCP<MLAPI::InverseOperator>> S_;
    };

    //! main routine analysis of FSIAMG
    void AnalyzeFSIAMG(const int myrank, const int snlevel, const int fnlevel, const int anlevel,
        Teuchos::ParameterList& sparams, Teuchos::ParameterList& fparams,
        Teuchos::ParameterList& aparams, std::vector<MLAPI::Operator>& Ass,
        std::vector<MLAPI::Operator>& Aff, std::vector<MLAPI::Operator>& Aaa,
        std::vector<MLAPI::Operator>& Pss, std::vector<MLAPI::Operator>& Rss,
        std::vector<MLAPI::Operator>& Pff, std::vector<MLAPI::Operator>& Rff,
        std::vector<MLAPI::Operator>& Paa, std::vector<MLAPI::Operator>& Raa,
        std::vector<MLAPI::Operator>& Asf, std::vector<MLAPI::Operator>& Afs,
        std::vector<MLAPI::Operator>& Afa, std::vector<MLAPI::Operator>& Aaf, ML* sml, ML* fml,
        ML* aml);

    //! analyze a single field
    void Analyse_SingleField(const std::string fieldname, const int myrank, std::string field,
        const int nlevel, Teuchos::ParameterList& params, std::vector<MLAPI::Operator>& A,
        std::vector<MLAPI::Operator>& P, std::vector<MLAPI::Operator>& R, AnalyzeBest& best);

    //! analyse BGS(AMG)
    void Analyse_BGSAMG(const int myrank, AnalyzeBest& sbest, AnalyzeBest& fbest,
        AnalyzeBest& abest, std::vector<MLAPI::Operator>& Ass, std::vector<MLAPI::Operator>& Pss,
        std::vector<MLAPI::Operator>& Rss, std::vector<MLAPI::Operator>& Aff,
        std::vector<MLAPI::Operator>& Pff, std::vector<MLAPI::Operator>& Rff,
        std::vector<MLAPI::Operator>& Aaa, std::vector<MLAPI::Operator>& Paa,
        std::vector<MLAPI::Operator>& Raa, std::vector<MLAPI::Operator>& Asf,
        std::vector<MLAPI::Operator>& Afs, std::vector<MLAPI::Operator>& Afa,
        std::vector<MLAPI::Operator>& Aaf);

    //! analyse AMG(BGS)
    void Analyse_AMGBGS(const int myrank, AnalyzeBest& sbest, AnalyzeBest& fbest,
        AnalyzeBest& abest, std::vector<MLAPI::Operator>& Ass, std::vector<MLAPI::Operator>& Pss,
        std::vector<MLAPI::Operator>& Rss, std::vector<MLAPI::Operator>& Aff,
        std::vector<MLAPI::Operator>& Pff, std::vector<MLAPI::Operator>& Rff,
        std::vector<MLAPI::Operator>& Aaa, std::vector<MLAPI::Operator>& Paa,
        std::vector<MLAPI::Operator>& Raa, std::vector<MLAPI::Operator>& Asf,
        std::vector<MLAPI::Operator>& Afs, std::vector<MLAPI::Operator>& Afa,
        std::vector<MLAPI::Operator>& Aaf);

    //! Compute rate of residual reduction
    double Rate(const int myrank, double t, double r, double initr, double l) const
    {
      // r /= sqrt(l); initr /= sqrt(l);
      double linrate = t * (r / initr);
      double lograte = 0.0;
      if (r / initr > 1.0)
        lograte = 2.0;
      else
        lograte = t / (-log(r / initr));
      double rate = FSIAMG_STRENGTH * linrate + (1. - FSIAMG_STRENGTH) * lograte;
#if (FSIAMG_ANALYSIS >= 4)
      if (r > initr)
        if (!myrank)
          printf("**residual increase in individual field** r0/r:     %8.4e/%8.4e\n", initr, r);
#endif
      return rate;
    }


    //! a single field single level Richardson iteration with smoother
    double RichardsonS(const std::string field, const int myrank, const int level, const int sweeps,
        const double damp, const MLAPI::Operator& A, const MLAPI::InverseOperator& S,
        MLAPI::MultiVector& x, const MLAPI::MultiVector& f, bool initiguesszero = false,
        bool analysis = false, bool silent = true) const;

    //! a non-Trilinos preconditioner Richardson
    double RichardsonMixed(const std::string field, const int myrank, const int level,
        const int sweeps, const double damp, const MLAPI::Operator& A,
        const LINALG::SparseMatrix& matrix, const Teuchos::RCP<LINALG::Preconditioner>& solver,
        MLAPI::MultiVector& x, const MLAPI::MultiVector& f, int& run, bool initiguesszero = false,
        bool analysis = false, bool silent = true) const;


    //! a single field V cycle analysis Richardson iteration
    double RichardsonV(const std::string field, const int myrank, int sweeps, const double damp,
        std::vector<int>& levelsweeps, std::vector<double>& leveldamps,
        std::vector<MLAPI::Operator>& A, std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S,
        std::vector<MLAPI::Operator>& P, std::vector<MLAPI::Operator>& R, const int level,
        const int nlevel, MLAPI::MultiVector& x, const MLAPI::MultiVector& f,
        bool initiguesszero = false, bool analysis = false, bool silent = false) const;

    //! a single field V cycle
    void Vcycle(const std::string field, const int myrank, std::vector<int>& sweeps,
        std::vector<double>& damps, const int level, const int nlevel, MLAPI::MultiVector& z,
        const MLAPI::MultiVector& b, const std::vector<MLAPI::Operator>& A,
        const std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S,
        const std::vector<MLAPI::Operator>& P, const std::vector<MLAPI::Operator>& R) const;



    /*! BGS(AMG) Richardson iteration with Block Gauss Seidel (BGS)
     *  and V cycle (V) for individual fields
     */
    double RichardsonBGS_V(const int myrank, const int sweeps, const double damp,
        std::vector<int>& blocksweeps, std::vector<double>& blockdamps, AnalyzeBest& sbest,
        AnalyzeBest& fbest, AnalyzeBest& abest, MLAPI::MultiVector& sy, MLAPI::MultiVector& fy,
        MLAPI::MultiVector& ay, const MLAPI::MultiVector& sf, const MLAPI::MultiVector& ff,
        const MLAPI::MultiVector& af, std::vector<MLAPI::Operator>& Ass,
        std::vector<MLAPI::Operator>& Pss, std::vector<MLAPI::Operator>& Rss,
        std::vector<MLAPI::Operator>& Aff, std::vector<MLAPI::Operator>& Pff,
        std::vector<MLAPI::Operator>& Rff, std::vector<MLAPI::Operator>& Aaa,
        std::vector<MLAPI::Operator>& Paa, std::vector<MLAPI::Operator>& Raa,
        std::vector<MLAPI::Operator>& Asf, std::vector<MLAPI::Operator>& Afs,
        std::vector<MLAPI::Operator>& Afa, std::vector<MLAPI::Operator>& Aaf,
        bool initiguesszero = false, bool analysis = false, bool silent = false) const;

    /*! BGS(AMG) Richardson iteration with Block Gauss Seidel (BGS)
     *  and V cycle (V) or other for individual fields
     */
    double RichardsonBGS_Mixed(const int myrank, const int sweeps, const double damp,
        std::vector<int>& blocksweeps, std::vector<double>& blockdamps, const bool sisamg,
        const bool fisamg, const bool aisamg, AnalyzeBest& sbest, AnalyzeBest& fbest,
        AnalyzeBest& abest, MLAPI::MultiVector& sy, MLAPI::MultiVector& fy, MLAPI::MultiVector& ay,
        const MLAPI::MultiVector& sf, const MLAPI::MultiVector& ff, const MLAPI::MultiVector& af,
        std::vector<MLAPI::Operator>& Ass, std::vector<MLAPI::Operator>& Pss,
        std::vector<MLAPI::Operator>& Rss, std::vector<MLAPI::Operator>& Aff,
        std::vector<MLAPI::Operator>& Pff, std::vector<MLAPI::Operator>& Rff,
        std::vector<MLAPI::Operator>& Aaa, std::vector<MLAPI::Operator>& Paa,
        std::vector<MLAPI::Operator>& Raa, std::vector<MLAPI::Operator>& Asf,
        std::vector<MLAPI::Operator>& Afs, std::vector<MLAPI::Operator>& Afa,
        std::vector<MLAPI::Operator>& Aaf, bool initiguesszero = false, bool analysis = false,
        bool silent = false) const;

    /*! Richardson iteration with Block Gauss Seidel (BGS)
     *  and Smoother (S) or Vcycle (V) for individual fields
     */
    double RichardsonBGS_SV(const int myrank, const int sweeps, const double damp,
        std::vector<int>* blocksweeps, std::vector<double>* blockdamps, const int level,
        AnalyzeBest& sbest, AnalyzeBest& fbest, AnalyzeBest& abest, MLAPI::MultiVector& sy,
        MLAPI::MultiVector& fy, MLAPI::MultiVector& ay, const MLAPI::MultiVector& sf,
        const MLAPI::MultiVector& ff, const MLAPI::MultiVector& af,
        std::vector<MLAPI::Operator>& Ass, std::vector<MLAPI::Operator>& Pss,
        std::vector<MLAPI::Operator>& Rss, std::vector<MLAPI::Operator>& Aff,
        std::vector<MLAPI::Operator>& Pff, std::vector<MLAPI::Operator>& Rff,
        std::vector<MLAPI::Operator>& Aaa, std::vector<MLAPI::Operator>& Paa,
        std::vector<MLAPI::Operator>& Raa, MLAPI::Operator& Asf, MLAPI::Operator& Afs,
        MLAPI::Operator& Afa, MLAPI::Operator& Aaf, bool initiguesszero = false,
        bool analysis = false, bool silent = false) const;

    /*! Richardson iteration with Schurcomplement preconditioner
     *  for structure split case
     */
    double RichardsonSchurStructureSplit_S(const int myrank, const int sweeps, const double damp,
        std::vector<int>* blocksweeps, std::vector<double>* blockdamps, const int level,
        AnalyzeBest& sbest, AnalyzeBest& fbest, AnalyzeBest& abest, MLAPI::MultiVector& sy,
        MLAPI::MultiVector& fy, MLAPI::MultiVector& ay, const MLAPI::MultiVector& sf,
        const MLAPI::MultiVector& ff, const MLAPI::MultiVector& af,
        std::vector<MLAPI::Operator>& Ass, std::vector<MLAPI::Operator>& Pss,
        std::vector<MLAPI::Operator>& Rss, std::vector<MLAPI::Operator>& Aff,
        std::vector<MLAPI::Operator>& Pff, std::vector<MLAPI::Operator>& Rff,
        std::vector<MLAPI::Operator>& Aaa, std::vector<MLAPI::Operator>& Paa,
        std::vector<MLAPI::Operator>& Raa, MLAPI::Operator& Asf, MLAPI::Operator& Afs,
        MLAPI::Operator& Afa, MLAPI::Operator& Aaf, bool initiguesszero = false,
        bool analysis = false, bool silent = false) const;

    /*! Richardson iteration with Schurcomplement preconditioner
     *  for fluid split case
     */
    double RichardsonSchurFluidSplit_S(const int myrank, const int sweeps, const double damp,
        std::vector<int>* blocksweeps, std::vector<double>* blockdamps, const int level,
        AnalyzeBest& sbest, AnalyzeBest& fbest, AnalyzeBest& abest, MLAPI::MultiVector& sy,
        MLAPI::MultiVector& fy, MLAPI::MultiVector& ay, const MLAPI::MultiVector& sf,
        const MLAPI::MultiVector& ff, const MLAPI::MultiVector& af,
        std::vector<MLAPI::Operator>& Ass, std::vector<MLAPI::Operator>& Pss,
        std::vector<MLAPI::Operator>& Rss, std::vector<MLAPI::Operator>& Aff,
        std::vector<MLAPI::Operator>& Pff, std::vector<MLAPI::Operator>& Rff,
        std::vector<MLAPI::Operator>& Aaa, std::vector<MLAPI::Operator>& Paa,
        std::vector<MLAPI::Operator>& Raa, MLAPI::Operator& Asf, MLAPI::Operator& Afs,
        MLAPI::Operator& Afa, MLAPI::Operator& Aaf, bool initiguesszero = false,
        bool analysis = false, bool silent = false) const;

    //! AMG(BGS) Richardson iteration with block V cycle for FSI
    double Richardson_BlockV(const int myrank, const int sweeps, const double damp,
        std::vector<int>& Vsweeps, std::vector<double>& Vdamps, std::vector<int>* blocksweeps,
        std::vector<double>* blockdamps, const std::vector<std::string>& blocksmoother,
        AnalyzeBest& sbest, AnalyzeBest& fbest, AnalyzeBest& abest, MLAPI::MultiVector& sy,
        MLAPI::MultiVector& fy, MLAPI::MultiVector& ay, const MLAPI::MultiVector& sf,
        const MLAPI::MultiVector& ff, const MLAPI::MultiVector& af,
        std::vector<MLAPI::Operator>& Ass, std::vector<MLAPI::Operator>& Pss,
        std::vector<MLAPI::Operator>& Rss, std::vector<MLAPI::Operator>& Aff,
        std::vector<MLAPI::Operator>& Pff, std::vector<MLAPI::Operator>& Rff,
        std::vector<MLAPI::Operator>& Aaa, std::vector<MLAPI::Operator>& Paa,
        std::vector<MLAPI::Operator>& Raa, std::vector<MLAPI::Operator>& Asf,
        std::vector<MLAPI::Operator>& Afs, std::vector<MLAPI::Operator>& Afa,
        std::vector<MLAPI::Operator>& Aaf, bool initiguesszero = false, bool analysis = false,
        bool silent = false) const;

    //! AMG(BGS) block V cycle for FSI
    void BlockVcycle(const int myrank, std::vector<int>& Vsweeps, std::vector<double>& Vdamps,
        std::vector<int>* blocksweeps, std::vector<double>* blockdamps,
        const std::vector<std::string>& blocksmoother, const int level, const int minnlevel,
        AnalyzeBest& sbest, AnalyzeBest& fbest, AnalyzeBest& abest, MLAPI::MultiVector& sy,
        MLAPI::MultiVector& fy, MLAPI::MultiVector& ay, const MLAPI::MultiVector& sf,
        const MLAPI::MultiVector& ff, const MLAPI::MultiVector& af,
        std::vector<MLAPI::Operator>& Ass, std::vector<MLAPI::Operator>& Pss,
        std::vector<MLAPI::Operator>& Rss, std::vector<MLAPI::Operator>& Aff,
        std::vector<MLAPI::Operator>& Pff, std::vector<MLAPI::Operator>& Rff,
        std::vector<MLAPI::Operator>& Aaa, std::vector<MLAPI::Operator>& Paa,
        std::vector<MLAPI::Operator>& Raa, std::vector<MLAPI::Operator>& Asf,
        std::vector<MLAPI::Operator>& Afs, std::vector<MLAPI::Operator>& Afa,
        std::vector<MLAPI::Operator>& Aaf) const;

    /// @name Access routines
    //!{

    //! Does solid use AMG?
    inline const bool& SisAMG() const { return sisml_; }

    //! Does fluid use AMG?
    inline const bool& FisAMG() const { return fisml_; }

    //! Does ALE use AMG?
    inline const bool& AisAMG() const { return aisml_; }
    //@}

    bool sisml_;  ///< solid uses AMG (true/false)
    bool fisml_;  ///< fluid uses AMG (true/false)
    bool aisml_;  ///< ALE uses AMG (true/false)
    int srun_;
    int frun_;
    int arun_;

    int minnlevel_;  ///< min of the below nlevel_
    int maxnlevel_;  ///< max of the below nlevel_
    int analyze_;    ///< run analysis of FSIAMG
    INPAR::FSI::LinearBlockSolver
        strategy_;  ///< type of preconditioner  to run: BGS(AMG) or AMG(BGS)
    std::vector<std::string> blocksmoother_;  ///< type of inter-field block smoother
    std::vector<double> schuromega_;  ///< damping factor for construction of Schur complement
    std::vector<double> pcomega_;
    std::vector<int> pciter_;
    std::vector<double> somega_;    ///< damping factors for solid AMG hierarchy
    std::vector<int> siterations_;  ///< number of sweeps for solid AMG hierarchy
    std::vector<double> fomega_;    ///< damping factors for fluid AMG hierarchy
    std::vector<int> fiterations_;  ///< number of sweeps for fluid AMG hierarchy
    std::vector<double> aomega_;    ///< damping factors for ALE AMG hierarchy
    std::vector<int> aiterations_;  ///< number of sweeps for ALE AMG hierarchy

    const INPAR::FSI::Verbosity verbosity_;  ///< verbosity level of FSI algorithm

    //! hybrid additive/multiplicative Schwarz preconditioner
    OverlappingBlockMatrixHybridSchwarz* hybridPrec_;

    /// @name AMG for structure
    int snlevel_;                     ///< num level in structure AMG
    Teuchos::ParameterList sparams_;  ///< parameter list
    mutable std::vector<MLAPI::Operator> Ass_;
    std::vector<Teuchos::RCP<MLAPI::InverseOperator>> Sss_;
    std::vector<MLAPI::Operator> Pss_;  ///< prolongation operators of solid AMG hierarchy
    std::vector<MLAPI::Operator> Rss_;  ///< restriction operators of solid AMG hierarchy
    //@}

    /// @name AMG for fluid
    int fnlevel_;  ///< num level in fluid AMG
    Teuchos::ParameterList fparams_;
    mutable std::vector<MLAPI::Operator> Aff_;
    std::vector<Teuchos::RCP<MLAPI::InverseOperator>> Sff_;
    mutable std::vector<MLAPI::Operator> Schurff_;
    std::vector<MLAPI::Operator> Pff_;  ///< prolongation operators of fluid AMG hierarchy
    std::vector<MLAPI::Operator> Rff_;  ///< restriction operators of fluid AMG hierarchy
    //@}

    /// @name AMG for ale
    int anlevel_;  ///< num level in ale AMG
    Teuchos::ParameterList aparams_;
    mutable std::vector<MLAPI::Operator> Aaa_;
    std::vector<Teuchos::RCP<MLAPI::InverseOperator>> Saa_;
    std::vector<MLAPI::Operator> Paa_;  ///< prolongation operators of ALE AMG hierarchy
    std::vector<MLAPI::Operator> Raa_;  ///< restriction operators of ALE AMG hierarchy
    //@}

    //! @name implicit off-diagonal blocks in block smoother
    //@{
    mutable MLAPI::Operator Asf_;
    mutable MLAPI::Operator Afs_;
    mutable MLAPI::Operator Afa_;
    mutable MLAPI::Operator Aaf_;
    //@}

    //! @name explicit coarse off-diagonal blocks in block smoother
    //@{
    mutable std::vector<MLAPI::Operator> ASF_;
    mutable std::vector<MLAPI::Operator> AFS_;
    mutable std::vector<MLAPI::Operator> AFA_;
    mutable std::vector<MLAPI::Operator> AAF_;
    //@}
  };


}  // namespace FSI

#endif
