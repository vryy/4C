/*---------------------------------------------------------------------*/
/*!
\file smoothing_strategy.H

\brief Add meaningful comment.

\level 3

\maintainer Alexander Popp

*/
/*---------------------------------------------------------------------*/
#ifndef SMOOTHING_STRATEGY_H_
#define SMOOTHING_STRATEGY_H_

#include "contact_abstract_strategy.H"
#include "../drt_mortar/mortar_utils.H"  // TODO: nötig?


// forward declaration
namespace MORTAR
{
  class MortarInterface;
}

namespace CONTACT
{
  class SmoothingStrategy : public CoAbstractStrategy
  {
   public:
    //! Standard constructor
    SmoothingStrategy(const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<CONTACT::CoInterface>> cinterface,
        std::vector<Teuchos::RCP<MORTAR::MortarInterface>> mtinterface, int dim,
        Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof);

    //! Shared data constructor
    SmoothingStrategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::CoInterface>> cinterface,
        std::vector<Teuchos::RCP<MORTAR::MortarInterface>> mtinterface, int dim,
        Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~SmoothingStrategy(){};

    //! @name Access methods

    /*!
    \brief Return convergence status of semi-smooth active set search

    If this Lagrange contact strategy is not based on a semi-smooth
    Newton approach, but on a fixed-point approach with two nested
    loops, then this method simply returns true, of course. Convergence
    of the active set is monitored with the flag activesetconv_ in
    this case and activesetssconv_ is meaningless.

    */
    bool ActiveSetSemiSmoothConverged() const
    {
      bool semismooth = DRT::INPUT::IntegralValue<int>(Params(), "SEMI_SMOOTH_NEWTON");
      if (semismooth)
        return activesetssconv_;
      else
        return true;
    }

    /*!
    \brief Return convergence status of fixed-point active set search

    If this Lagrange contact strategy is based on a semi-smooth
    Newton approach and not on a fixed-point approach with two nested
    loops, then this method simply returns true, of course. Convergence
    of the active set is monitored with the flag activesetssconv_ in
    this case and activesetconv_ is meaningless.

    */
    bool ActiveSetConverged()
    {
      bool semismooth = DRT::INPUT::IntegralValue<int>(Params(), "SEMI_SMOOTH_NEWTON");
      if (!semismooth)
        return activesetconv_;
      else
        return true;
    }

    /*!
    \brief Return no. of fixed-point active sets in this time step

    */
    int ActiveSetSteps() { return activesetsteps_; }

    //! Get the active node row map of the previous Newton step
    virtual Teuchos::RCP<const Epetra_Map> GetOldActiveRowNodes() const
    {
      dserror("Currently unsupported for the smoothing strategy!");
      return Teuchos::null;
    };

    //! Get the slip node row map of the previous Newton step
    virtual Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const
    {
      dserror("Currently unsupported for the smoothing strategy!");
      return Teuchos::null;
    };

    //@}

    //! @name Evaluation methods

    /*!
    \brief Initialize general contact variables for next Newton step

    For a lagrangian strategy this includes the global normal / tangent matrices N and T,
    the global derivative matrices S and P and Tresca friction matrix L + vector r.

    */
    void Initialize();

    /*!
    \brief Evaluate contact

    For a lagrangian strategy this involves heavy modification to the initial kteff and feff.
    Hence, they are infact build from scratch here.
    The application of modifications to groups of dofs (slave, master, active etc.)
    results in some matrix and vector splitting and a lot of matrix-vector calculation in here!

    */
    void EvaluateContact(
        Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff);

    /*!
    \brief Evaluate frictional contact

    */
    void EvaluateFriction(
        Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff);

    /*!
    \brief Build 2x2 saddle point system

    \param kdd (in): the displacement dof stiffness (upper left block)
    \param fd (in): the displacement dof r.h.s. (upper block)
    \param sold (in): the displacement dof solution increment
    \param dirichtoggle (in): toggle vector for dirichlet conditions
    \param blockMat (out): Epetra_Operator containing the 2x2 block sparse matrix object
    \param mergedsol (out): Epetra_Vector for merged solution vector
    \param mergedrhs (out): Epetra_Vector for merged right hand side vector
    */
    void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) override;

    /*!
    \brief Update internal member variables after solving the 2x2 saddle point contact system

    \param sold (out): the displacement dof solution increment (associated with displacement dofs)
    \param mergedsol (in): Epetra_Vector for merged solution vector (containing the new solution
    vector of the full merged linear system)
    */
    void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) override;

    /*!
    \brief The entries of the constraint right-hand side are
    calculated. This function is outsourced cause
    the vector is needed for the line search algorithm.

    */
    void EvalConstrRHS();

    /*!
    \brief evaluate mortar coupling (MT)

    */
    void MortarCoupling(const Teuchos::RCP<const Epetra_Vector>& dis);

    /*!
    \brief Recovery method

    We only recover the Lagrange multipliers here, which had been
    statically condensated during the setup of the global problem!
    Optionally satisfaction or violation of the contact boundary
    conditions can be checked, too.

    */
    void Recover(Teuchos::RCP<Epetra_Vector> disi);

    /*!
    \brief Update active set and check for convergence

    In this function we loop over all interfaces and then over all
    slave nodes to check, whether the assumption of them being active
    or inactive respectively has been correct. If a single node changes
    state, the active set is adapted accordingly and the convergence
    flag is kept on false. Here we have the "standard" case of two
    nested iteration loops, and as a consequence this method is
    called AFTER convergence of the inner Newton iteration. If there
    is a change in the active set, another full Newton iteration has
    to be performed for the current time / load step.

    */
    void UpdateActiveSet() { return; };  // TODO

    /*!
    \brief Update active set and check for convergence

    In this function we loop over all interfaces and then over all
    slave nodes to check, whether the assumption of them being active
    or inactive respectively has been correct. If a single node changes
    state, the active set is adapted accordingly and the convergence
    flag is kept on false.

    Here we have the semi-smooth Newton case
    with one combined iteration loop for active set search and large
    deformations. As a consequence this method is called AFTER each
    (not yet converged) Newton step. If there is a change in the active
    set or the residual and disp norm are still above their limits,
    another Newton step has to be performed.

    \note We use the flag \c firstStepPredictor to overwrite the active set status
    for each node in the predictor of the first time step

    \param[in] firstStepPredictor Boolean flag to indicate the predictor step in the first time step
    */
    void UpdateActiveSetSemiSmooth(const bool firstStepPredictor = false);

    /*!
    \brief Reset active set status for next time step

    */
    void ResetActiveSet()
    {
      activesetssconv_ = false;
      activesetconv_ = false;
      activesetsteps_ = 1;
      return;
    }

    /*!
    \brief Set current displ state [derived]

    */
    virtual void SetState(const enum MORTAR::StateType& statetype, const Epetra_Vector& vec);

    /*!
    \brief Return matrix T

    */
    Teuchos::RCP<LINALG::SparseMatrix> TMatrix() { return tmatrix_; }

    //@}

    //! @name Debugging and visualization methods

    /*!
    \brief Check linear and angular momentum conservation

    */
    void CheckConservationLaws(const Epetra_Vector& fs, const Epetra_Vector& fm) { return; };


    //@}

    //! @name Empty methods

    /*!
    \brief Empty methods only relevant for other strategies

    For a Lagrange strategy these are functions without functionality.
    Call them whenever you like.

    */
    double ConstraintNorm() const { return 0.0; }
    void EvaluateRelMovPredict() {}
    double InitialPenalty() { return 0.0; }
    void InitializeUzawa(
        Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
    {
    }
    void ResetPenalty() {}
    void ModifyPenalty() {}
    void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) {}
    void UpdateUzawaAugmentedLagrange() {}
    void UpdateConstraintNorm(int uzawaiter = 0) {}
    bool IsPenalty() const override { return false; };

   protected:
    //! derived
    virtual std::vector<Teuchos::RCP<CONTACT::CoInterface>>& Interfaces() { return interface_; }

    //! derived
    virtual const std::vector<Teuchos::RCP<CONTACT::CoInterface>>& Interfaces() const
    {
      return interface_;
    }

   private:
    /*!
    \brief create vector with node positions

    */
    void AssembleCoords(const std::string& sidename, bool ref, Teuchos::RCP<Epetra_Vector> vec);

    /*!
    \brief Setup mortar maps for strategy

    */
    void SetupMT(bool redistributed);

   private:
    // don't want = operator and cctor
    SmoothingStrategy operator=(const SmoothingStrategy& old) = delete;
    SmoothingStrategy(const SmoothingStrategy& old) = delete;

    // Store Coupling Matrices in case of Poro Lagrange Strategy ... here just ignore!
    virtual void SaveCouplingMatrices(Teuchos::RCP<LINALG::SparseMatrix> dhat,
        Teuchos::RCP<LINALG::SparseMatrix> mhataam, Teuchos::RCP<LINALG::SparseMatrix> invda)
    {
      return;
    }

    bool activesetssconv_;  // convergence flag for semi-smooth active set search
    bool activesetconv_;    // convergence flag for fixed-point active set search
    int activesetsteps_;    // number of fixed-point active set steps in this time step

    std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface_;  //!< contact interfaces

    Teuchos::RCP<LINALG::SparseMatrix> mhatmatrix_;  // product of global Mortar matrices inv(D)*M

    Teuchos::RCP<LINALG::SparseMatrix> tmatrix_;  // global Matrix T containing active node tangents
    Teuchos::RCP<LINALG::SparseMatrix>
        smatrix_;  // global Matrix S containing normal+D+M derivatives
    Teuchos::RCP<LINALG::SparseMatrix> smatrixW_;  // global Matrix S containing W derivatives
    Teuchos::RCP<LINALG::SparseMatrix>
        tderivmatrix_;  // global Matrix P containing tangent derivatives


    Teuchos::RCP<Epetra_Vector> fs_;           // slave side effective forces (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> invd_;  // inverse of Mortar matrix D (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> ksn_;   // stiffness block K_sn (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> ksm_;   // stiffness block K_sm (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> kss_;   // stiffness block K_ss (needed for LM)

    Teuchos::RCP<LINALG::SparseMatrix>
        linslipLM_;  // global matrix containing derivatives (LM) of slip condition
    Teuchos::RCP<LINALG::SparseMatrix>
        linslipDIS_;  // global matrix containing derivatives (DIS) of slip condition
    Teuchos::RCP<Epetra_Vector> linslipRHS_;  // r.h.s vector friction slip nodes
    Teuchos::RCP<LINALG::SparseMatrix>
        linstickLM_;  // global matrix containing derivatives (LM) of slip condition
    Teuchos::RCP<LINALG::SparseMatrix>
        linstickDIS_;  // global matrix containing derivatives (DIS) of stick condition
    Teuchos::RCP<Epetra_Vector> linstickRHS_;  // r.h.s vector for friction stick condition

    Teuchos::RCP<Epetra_Map> zigzagone_;    // active node set of last active set try
    Teuchos::RCP<Epetra_Map> zigzagtwo_;    // active node set of second-last active set try
    Teuchos::RCP<Epetra_Map> zigzagthree_;  // active node set of third-last active set try

    //******************************************
    // meshtying stuff
    std::vector<Teuchos::RCP<MORTAR::MortarInterface>> minterface_;  // mt interfaces

    Teuchos::RCP<LINALG::SparseMatrix> MT_dmatrix_;  // Mortar matrix D for global mt
    Teuchos::RCP<LINALG::SparseMatrix> MT_mmatrix_;  // Mortar matrix M for global mt
    Teuchos::RCP<Epetra_Vector> MT_g_;               // global weighted gap vector g for global mt

    Teuchos::RCP<LINALG::SparseMatrix> mt1_dmatrix_;  // Mortar matrix D for mt1
    Teuchos::RCP<LINALG::SparseMatrix> mt1_mmatrix_;  // Mortar matrix M for mt1

    Teuchos::RCP<LINALG::SparseMatrix> mt2_dmatrix_;  // Mortar matrix D for mt2
    Teuchos::RCP<LINALG::SparseMatrix> mt2_mmatrix_;  // Mortar matrix M for mt2

    Teuchos::RCP<Epetra_Vector> MT_z_;     // current vector of Lagrange multipliers at t_n+1
    Teuchos::RCP<Epetra_Vector> MT_zold_;  // old vector of Lagrange multipliers at t_n
    Teuchos::RCP<Epetra_Vector>
        MT_zincr_;  // Lagrange multiplier vector increment within SaddlePointSolve (this is NOT the
                    // increment of z_ between t_{n+1} and t_{n}!)

    Teuchos::RCP<Epetra_Vector> mt1_z_;     // current vector of Lagrange multipliers at t_n+1
    Teuchos::RCP<Epetra_Vector> mt1_zold_;  // old vector of Lagrange multipliers at t_n
    Teuchos::RCP<Epetra_Vector>
        mt1_zincr_;  // Lagrange multiplier vector increment within SaddlePointSolve (this is NOT
                     // the increment of z_ between t_{n+1} and t_{n}!)

    Teuchos::RCP<Epetra_Vector> mt2_z_;     // current vector of Lagrange multipliers at t_n+1
    Teuchos::RCP<Epetra_Vector> mt2_zold_;  // old vector of Lagrange multipliers at t_n
    Teuchos::RCP<Epetra_Vector>
        mt2_zincr_;  // Lagrange multiplier vector increment within SaddlePointSolve (this is NOT
                     // the increment of z_ between t_{n+1} and t_{n}!)

    int globaloffset_;  // contact lm gids

    Teuchos::RCP<Epetra_Map>
        MT_glmdofrowmap_;  // global MT Lagrange mult. dof row map (of all minterfaces)

    Teuchos::RCP<Epetra_Map>
        MT_gsdofrowmap_;  // global MT slave dof row map (of all minterfaces) = contact s+m
    Teuchos::RCP<Epetra_Map> MT_gmdofrowmap_;  // global MT master dof row map (of all minterfaces)

    Teuchos::RCP<Epetra_Map>
        MT_gsnoderowmap_;  // global MT slave node row map (of all minterfaces) = contact s+m
    Teuchos::RCP<Epetra_Map>
        MT_gmnoderowmap_;  // global MT master node row map (of all minterfaces)

    Teuchos::RCP<Epetra_Map>
        MT_gsmdofrowmap_;  // global MT slave + master dof row map (of all minterface



    // penalty
    Teuchos::RCP<LINALG::SparseMatrix> linzmatrix_;  // global matrix LinZ with derivatives of LM
    double constrnorm_;                              // L2-norm of normal contact constraints
    double constrnormtan_;                           // L2-norm of tangential contact constraints
    double initialpenalty_;                          // initial penalty parameter
    double initialpenaltytan_;                       // initial tangential penalty parameter

  };  // class SmoothingStrategy
}  // namespace CONTACT

#endif /* SMOOTHING_STRATEGY_H_ */
