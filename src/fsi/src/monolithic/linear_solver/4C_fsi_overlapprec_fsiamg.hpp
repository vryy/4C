/*----------------------------------------------------------------------*/
/*! \file

\brief Special version of block matrix that includes the FSI block preconditioner


\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_OVERLAPPREC_FSIAMG_HPP
#define FOUR_C_FSI_OVERLAPPREC_FSIAMG_HPP

#include "4C_config.hpp"

#include "4C_fsi_overlapprec.hpp"

#include <MLAPI_LoadBalanceInverseOperator.h>
#include <MLAPI_Operator.h>

#define FSIAMG_STRENGTH \
  0.85                     ///< emphasis on strength instead of speed on a 0...1 scale in analysis
#define FSIAMG_ANALYSIS 0  ///< level 0(off),1,2,3,4 switch on live analysis during actual run

FOUR_C_NAMESPACE_OPEN

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
    OverlappingBlockMatrixFSIAMG(const Core::LinAlg::MultiMapExtractor& maps,
        Adapter::FSIStructureWrapper& structure, Adapter::Fluid& fluid, Adapter::AleFsiWrapper& ale,
        bool structuresplit, int symmetric, std::vector<std::string>& blocksmoother,
        std::vector<double>& schuromega, std::vector<double>& omega, std::vector<int>& iterations,
        std::vector<double>& somega, std::vector<int>& siterations, std::vector<double>& fomega,
        std::vector<int>& fiterations, std::vector<double>& aomega, std::vector<int>& aiterations,
        int analyze, Inpar::FSI::LinearBlockSolver strategy, Inpar::FSI::Verbosity verbosity,
        OverlappingBlockMatrixHybridSchwarz* hybridPrec = nullptr);

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    const char* Label() const override;

    //@}

    /// setup of block preconditioners
    void SetupPreconditioner() override;

   protected:
    /// symmetric Gauss-Seidel block preconditioner
    void sgs(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    /// do list for MLAPI smoother
    void select_mlapi_smoother(std::string& type, const int level, Teuchos::ParameterList& subp,
        Teuchos::ParameterList& p, Teuchos::ParameterList& pushlist);

    /// wrap ILU smoother from ML
    void wrap_ilu_smoother(
        ML* ml, MLAPI::Operator& A, MLAPI::LoadBalanceInverseOperator& S, const int level);


    /// generic Vcycle that works on all fields
    virtual void vcycle(const int level, const int nlevel, MLAPI::MultiVector& z,
        const MLAPI::MultiVector& b, const std::vector<MLAPI::Operator>& A,
        const std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S,
        const std::vector<MLAPI::Operator>& P, const std::vector<MLAPI::Operator>& R,
        const bool trigger = false) const;


    /// FSIAMG multigrid with explicit coarse off-diagonals
    virtual void explicit_block_vcycle(const int level, const int nlevel, MLAPI::MultiVector& mlsy,
        MLAPI::MultiVector& mlfy, MLAPI::MultiVector& mlay, const MLAPI::MultiVector& mlsx,
        const MLAPI::MultiVector& mlfx, const MLAPI::MultiVector& mlax) const;


    /// block Gauss-Seidel smoother within one level (explicit off-diagonals)
    virtual void explicit_block_gauss_seidel_smoother(const int level, MLAPI::MultiVector& mlsy,
        MLAPI::MultiVector& mlfy, MLAPI::MultiVector& mlay, const MLAPI::MultiVector& mlsx,
        const MLAPI::MultiVector& mlfx, const MLAPI::MultiVector& mlax, const bool amgsolve) const;


    /// iterate on the field individual blocks within the block Gauss Seidel smoother
    virtual void local_block_richardson(const int iterations, const double omega, const int level,
        const bool amgsolve, const int nlevel, MLAPI::MultiVector& z, const MLAPI::MultiVector& b,
        const std::vector<MLAPI::Operator>& A,
        const std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S,
        const std::vector<MLAPI::Operator>& P, const std::vector<MLAPI::Operator>& R) const;


    /// triple matrix product with fine level operator only
    void ra_pfine(MLAPI::Operator& RAP, const MLAPI::Operator& R, Teuchos::RCP<Epetra_CrsMatrix> A,
        const MLAPI::Operator& P);

    /// triple matrix product with coarse level operators
    void ra_pcoarse(MLAPI::Operator& RAP, const MLAPI::Operator& R, const MLAPI::Operator& A,
        const MLAPI::Operator& P);

    /// build off-diagonal coupling blocks for FSIAMG
    void ra_poffdiagonals();

    /// build Schur Complement operator from fluid block
    void schur_complement_operator(MLAPI::Operator& Schur, MLAPI::Operator& Ass,
        MLAPI::Operator& Aff, MLAPI::Operator& Aaa, MLAPI::Operator& Asf, MLAPI::Operator& Afs,
        MLAPI::Operator& Afa, MLAPI::Operator& Aaf, const double omega, const bool structuresplit);

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
            s_(7, Teuchos::null)
      {
        if (nlevel - 1 > 6) FOUR_C_THROW("Can only analyze V cycles upto 7 levels");
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
      inline std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S() { return s_; }

     private:
      int nlevel_;                         ///< number of levels
      std::vector<std::string> besttype_;  //< best set of algorithms
      std::vector<double> bestdamp_;       ///< best set of damping parameters
      std::vector<int> bestpoly_;          ///< best set of polynomial degrees
      std::vector<int> bestsweeps_;        ///< best set of number of sweeps
      std::vector<Teuchos::RCP<MLAPI::InverseOperator>> s_;
    };

    //! Compute rate of residual reduction
    double rate(const int myrank, double t, double r, double initr, double l) const
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
    double richardson_s(const std::string field, const int myrank, const int level,
        const int sweeps, const double damp, const MLAPI::Operator& A,
        const MLAPI::InverseOperator& S, MLAPI::MultiVector& x, const MLAPI::MultiVector& f,
        bool initiguesszero = false, bool analysis = false, bool silent = true) const;

    //! a non-Trilinos preconditioner Richardson
    double richardson_mixed(const std::string field, const int myrank, const int level,
        const int sweeps, const double damp, const MLAPI::Operator& A,
        const Core::LinAlg::SparseMatrix& matrix,
        const Teuchos::RCP<Core::LinAlg::Preconditioner>& solver, MLAPI::MultiVector& x,
        const MLAPI::MultiVector& f, int& run, bool initiguesszero = false, bool analysis = false,
        bool silent = true) const;


    //! a single field V cycle analysis Richardson iteration
    double richardson_v(const std::string field, const int myrank, int sweeps, const double damp,
        std::vector<int>& levelsweeps, std::vector<double>& leveldamps,
        std::vector<MLAPI::Operator>& A, std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S,
        std::vector<MLAPI::Operator>& P, std::vector<MLAPI::Operator>& R, const int level,
        const int nlevel, MLAPI::MultiVector& x, const MLAPI::MultiVector& f,
        bool initiguesszero = false, bool analysis = false, bool silent = false) const;

    //! a single field V cycle
    void vcycle(const std::string field, const int myrank, std::vector<int>& sweeps,
        std::vector<double>& damps, const int level, const int nlevel, MLAPI::MultiVector& z,
        const MLAPI::MultiVector& b, const std::vector<MLAPI::Operator>& A,
        const std::vector<Teuchos::RCP<MLAPI::InverseOperator>>& S,
        const std::vector<MLAPI::Operator>& P, const std::vector<MLAPI::Operator>& R) const;



    /*! BGS(AMG) Richardson iteration with Block Gauss Seidel (BGS)
     *  and V cycle (V) for individual fields
     */
    double richardson_bgs_v(const int myrank, const int sweeps, const double damp,
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
    double richardson_bgs_mixed(const int myrank, const int sweeps, const double damp,
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

    /// @name Access routines
    //!{

    //! Does solid use AMG?
    inline const bool& sis_amg() const { return sisml_; }

    //! Does fluid use AMG?
    inline const bool& fis_amg() const { return fisml_; }

    //! Does ALE use AMG?
    inline const bool& ais_amg() const { return aisml_; }
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
    Inpar::FSI::LinearBlockSolver
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

    const Inpar::FSI::Verbosity verbosity_;  ///< verbosity level of FSI algorithm

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

FOUR_C_NAMESPACE_CLOSE

#endif
