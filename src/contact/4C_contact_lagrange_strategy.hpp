/*---------------------------------------------------------------------*/
/*! \file
\brief Contact solving strategy with (standard/dual) Lagrangian multipliers.

\level 2


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_LAGRANGE_STRATEGY_HPP
#define FOUR_C_CONTACT_LAGRANGE_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_contact_abstract_strategy.hpp"

FOUR_C_NAMESPACE_OPEN


namespace CONTACT
{
  /*! \brief Contact solving strategy with (standard/dual) Lagrangian multipliers.

   This is a specialization of the abstract contact algorithm as defined in AbstractStrategy.
   For a more general documentation of the involved functions refer to AbstractStrategy.
  */

  class LagrangeStrategy : public AbstractStrategy
  {
   public:
    /*!
    \brief Shared data constructor

    Creates the strategy base object and initializes all global variables.

    \param[in] data_ptr Data container object
    \param[in] DofRowMap Dof row map of underlying problem
    \param[in] NodeRowMap Node row map of underlying problem
    \param[in] params List of contact/parameters
    \param[in] interface All contact interface objects
    \param[in] spatialDim Spatial dimension of the problem
    \param[in] comm Communicator
    \param[in] alphaf Mid-point for Generalized-alpha time integration
    \param[in] maxdof Highest DOF number in global problem
    */
    LagrangeStrategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, const int spatialDim,
        Teuchos::RCP<const Epetra_Comm> comm, const double alphaf, const int maxdof);



    //! @name Access methods

    /*!
    \brief Return convergence status of semi-smooth active set search

    If this Lagrange contact strategy is not based on a semi-smooth
    Newton approach, but on a fixed-point approach with two nested
    loops, then this method simply returns true, of course. Convergence
    of the active set is monitored with the flag activesetconv_ in
    this case and activesetssconv_ is meaningless.

    */
    bool active_set_semi_smooth_converged() const override
    {
      bool semismooth = CORE::UTILS::IntegralValue<int>(Params(), "SEMI_SMOOTH_NEWTON");
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
    bool ActiveSetConverged() override
    {
      bool semismooth = CORE::UTILS::IntegralValue<int>(Params(), "SEMI_SMOOTH_NEWTON");
      if (!semismooth)
        return activesetconv_;
      else
        return true;
    }

    /*!
    \brief Return no. of fixed-point active sets in this time step

    */
    int ActiveSetSteps() override { return activesetsteps_; }

    /*! \brief Return the desired right-hand-side block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired vector block type, e.g. displ, constraint,*/
    Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum CONTACT::VecBlockType& bt) const override;


    /*! \brief Recover the current state
     *
     * The main task of this method is to recover the Lagrange multiplier solution.
     * The Lagrange multiplier solution will be stored inside the corresponding strategy
     * and is necessary for different internal evaluation methods. If the Lagrange multiplier
     * is condensed, this method is the right place to recover it from the displacement solution.
     * If it is not condensed (saddle-point system) use the ResetLagrangeMultiplier routine
     * instead.
     *
     * \param cparams (in): parameter interface between the contact objects and the structural time
     *                       integration
     * \param xold (in): old solution vector of the NOX solver
     * \param dir (in): current search direction
     * \param xnew (in): new solution vector of the NOX solver
     *
     * \note The search direction \c dir in general differs from the actual step since the step
     * length can be differ from 1.0.
     * */
    void RunPostComputeX(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
        const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

    /*! \brief Return the desired matrix block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired matrix block type, e.g. displ_displ, displ_lm, ...
     *  \param cparams (in): contact parameter interface (read-only) */
    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetMatrixBlockPtr(const enum CONTACT::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams = nullptr) const override;

    /*! \brief Apply modifications (e.g. condensation) directly before linear solve
     *
     * \todo Complete documentation of parameters.
     *
     * @param kteff Jacobian matrix
     * @param rhs ?? Right-hand side of linear system or residual?
     */
    void run_pre_apply_jacobian_inverse(
        Teuchos::RCP<CORE::LINALG::SparseMatrix> kteff, Epetra_Vector& rhs) override;

    /*! \brief Perform condensation of frictionless contact
     *
     * \todo Complete documentation of parameters.
     *
     * @param kteff Jacobian matrix
     * @param rhs ?? Right-hand side of linear system or residual?
     */
    virtual void condense_frictionless(
        Teuchos::RCP<CORE::LINALG::SparseMatrix> kteff, Epetra_Vector& rhs);

    /*! \brief Perform condensation of frictional contact
     *
     * \todo Complete documentation of parameters.
     *
     * @param kteff Jacobian matrix
     * @param rhs ?? Right-hand side of linear system or residual?
     */
    virtual void CondenseFriction(
        Teuchos::RCP<CORE::LINALG::SparseMatrix> kteff, Epetra_Vector& rhs);

    /*! Recover condensed Lagrange multiplier after linear solve
     *
     * \todo Complete documentation of parameters.
     *
     * @param cparams Parameter interface between contact and structural time integration
     * @param rhs ?? Right-hand side of linear system or residual?
     * @param result ??
     * @param xold ??
     * @param grp ??
     */
    void run_post_apply_jacobian_inverse(const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& rhs, Epetra_Vector& result, const Epetra_Vector& xold,
        const NOX::NLN::Group& grp) override;

    //! The contributions to the structural right-hand-side block are calculated.
    virtual void EvalStrContactRHS();

    //! Assemble contact contributions to the rhs
    virtual void AssembleContactRHS();

    //! Assemble all contact contributions
    virtual void assemble_all_contact_terms();

    //! Assemble all contact contributions (frictionless contact)
    virtual void assemble_all_contact_terms_frictionless();

    //! Assemble all contact contributions (frictional contact)
    virtual void assemble_all_contact_terms_friction();

    const Epetra_Map& SlNormalDoFRowMap(const bool& redist) const override
    {
      if ((not redist) and ParRedist())
        FOUR_C_THROW("The original / not redistributed slave normal row map is not available!");

      return *gsdofrowmap_;
    }

    //! Get the active node row map of the previous Newton step
    Teuchos::RCP<const Epetra_Map> get_old_active_row_nodes() const override
    {
      return gOldActiveSlaveNodes_;
    };

    //! Get the slip node row map of the previous Newton step
    Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const override { return gOldslipnodes_; };

    //@}

    //! @name Evaluation methods

    /*!
    \brief Build 2x2 saddle point system

    The saddle-point system reads
    \f[
    \left[\begin{array}{cc}
    K_{dd} & K_{dz}\\
    K_{zd} & K_{zz}
    \end{array}\right]
    \left[\begin{array}{c}\Delta d\\ \Delta z\end{array}\right]
    = - \left[\begin{array}{c}f_d\\ f_z\end{array}\right]
    \f]
    with \f$d\f$ and \f$z\f$ denoting the displacement and Lagrange multiplier degrees of freedom
    (DOFs), respectively.

    In terms of data structures and parallel layout, the Lagrange multiplier DOFs have the same
    parallel layout as their slave side displacement DOFs counterparts, but distinct global IDs
    (GIDs). For evaluation of mortar matrices, the salve side DofRowMap #gsdofrowmap_ is used.
    However, for assembling the saddle-point system, the constraint equations need their own unique
    map, #glmdofrowmap_. To also account for the possibility of a parallel redistribution of the
    Mortar interface discretizations, the assembly is performed in two steps:

    1. Transform maps related to constraint DOFs of matrix blocks \f$K_{dz}\f$, \f$K_{zd}\f$, and
       \f$K_{zz}\f$ to the Lagrange multiplier DofRowMap #glmdofrowmap_.
    2. If parallel redistribution is active, do another transformation to go from the redistributed
       constraint DofRowMap #glmdofrowmap_ to the original distribution #pglmdofrowmap_.

    \param[in] kdd Displacement dof stiffness (upper left block)
    \param[in] fd  Displacement dof r.h.s. (upper block)
    \param[in] sold Displacement dof solution increment
    \param[in] dbcmaps Map extractor to apply Dirichlet boundary conditions
    \param[out] blockMat Epetra_Operator containing the 2x2 block sparse matrix object
    \param[out] blocksol Epetra_Vector for merged solution vector
    \param[out] blockrhs Epetra_Vector for merged right hand side vector

    \todo Update documentation once class members have been removed and replaced by accessing the
    data container directly.
    \todo Argument \c numiter not used. Can it be removed?

    \sa update_displacements_and_l_mincrements

    \author Tobias Wiesner \date 11/2014
    */
    void build_saddle_point_system(Teuchos::RCP<CORE::LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) override;

    /*!
    \brief Update internal member variables after solving the 2x2 saddle point contact system

    We have to extract the displacement and Lagrange multiplier solution increment from the blocked
    solution vector. This requires to also transfer constraint GIDs back to #gsdofrowmap_ since
    constraints internally use the same GIDs as their displacement DOF counterparts.

    \note In case of parallel redistribution, we also have to transfer constraint data from the
    undredistributed layout (used to build the saddle-point system) to the redistributed layout
    (used to perform contact evaluations).

    \note More details on handling of maps w/ and w/o parallel redistribution are given in
    build_saddle_point_system().

    \param[out] sold Displacement dof solution increment (associated with displacement dofs)
    \param[in] blocksol Merged solution vector (containing the new solution vector of the full
                    merged linear system, i.e. displacement and Lagrange multiplier DOFs)

    \sa build_saddle_point_system

    \author Tobias Wiesner \date 11/2014
    */
    void update_displacements_and_l_mincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) override;

    /*!
    \brief The entries of the constraint right-hand side are calculated.

    This function is outsourced cause the vector is needed for the line search algorithm.

    */
    void EvalConstrRHS() override;

    /*! \brief Compute force and stiffness terms
     *
     * \param cparams (in): parameter interface between the contact objects and the structural time
     * integration*/
    void EvalForceStiff(CONTACT::ParamsInterface& cparams) override;

    /*! \brief Run at the beginning of the Evaluate() routine
     *         set force evaluation flag
     *
     */
    void PreEvaluate(CONTACT::ParamsInterface& cparams) override;

    /*! \brief Run in the end of the Evaluate() routine to reset
     *         force evaluation flag
     */
    void PostEvaluate(CONTACT::ParamsInterface& cparams) override;

    /*! \brief This is a postprocessing functionality for nonsmooth contact
     */
    void compute_contact_stresses() final;

    /*!
    \brief Recovery method

    We only recover the Lagrange multipliers here, which had been
    statically condensed during the setup of the global problem!
    Optionally satisfaction or violation of the contact boundary
    conditions can be checked, too.

    */
    void Recover(Teuchos::RCP<Epetra_Vector> disi) override;

    /*! \brief Reset the internal stored Lagrange multipliers
     *
     * \param cparams (in): parameter interface between the contact objects and the structural time
     *                      integration
     * \param xnew    (in): new solution vector of the NOX solver
     */
    void reset_lagrange_multipliers(
        const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew) override;

    /*! \brief Compute force terms
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     * integration*/
    void EvalForce(CONTACT::ParamsInterface& cparams) override;

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
    void UpdateActiveSet() override;

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
    set or the residual and displacement increment norm are still above their limits,
    another Newton step has to be performed.

    \note We use the flag \c firstStepPredictor to overwrite the active set status
    for each node in the predictor of the first time step (or the first time step
    after a restart).

    \param[in] firstStepPredictor Boolean flag to indicate the predictor step in the first time step
    */
    void update_active_set_semi_smooth(const bool firstStepPredictor = false) override;

    /*!
    \brief Reset active set status for next time step

    */
    void ResetActiveSet() override
    {
      activesetssconv_ = false;
      activesetconv_ = false;
      activesetsteps_ = 1;
      return;
    }

    /*!
    \brief Return matrix T

    */
    Teuchos::RCP<CORE::LINALG::SparseMatrix> TMatrix() override { return tmatrix_; }

    //@}

    //! @name Debugging and visualization methods

    /*!
    \brief Check linear and angular momentum conservation

    */
    void check_conservation_laws(const Epetra_Vector& fs, const Epetra_Vector& fm);

    /*!
    \brief do additional matrix manipulations for regularization scaling

    */
    void do_regularization_scaling(bool aset, bool iset,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& invda,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& kan,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& kam,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& kai,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& kaa, Teuchos::RCP<Epetra_Vector>& fa,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& kteffnew, Teuchos::RCP<Epetra_Vector>& feffnew);

    /*!
    \brief calculate regularization scaling and apply it to matrixes

    */
    void evaluate_regularization_scaling(Teuchos::RCP<Epetra_Vector> gact);

    /*!
    \brief Saving reference state is required for penalty support (LTL)

    */
    void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) override;
    //@}

    //! @name Empty methods

    /*!
    \brief Empty methods only relevant for other strategies

    For a Lagrange strategy these are functions without functionality.
    Call them whenever you like.

    */
    double ConstraintNorm() const override { return 0.0; }
    void evaluate_rel_mov_predict() override {}
    double InitialPenalty() override { return 0.0; }
    void InitializeUzawa(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override
    {
    }
    void ResetPenalty() override {}
    void ModifyPenalty() override {}
    void update_uzawa_augmented_lagrange() override {}
    void update_constraint_norm(int uzawaiter = 0) override {}
    bool IsPenalty() const override { return false; };

   protected:
    //! derived
    std::vector<Teuchos::RCP<CONTACT::Interface>>& Interfaces() override { return interface_; }

    //! derived
    const std::vector<Teuchos::RCP<CONTACT::Interface>>& Interfaces() const override
    {
      return interface_;
    }

    /*!
    \brief Initialize general contact variables for next Newton step

    For a lagrangian strategy this includes the global normal / tangent matrices N and T,
    the global derivative matrices S and P and Tresca friction matrix L + vector r.

    */
    void Initialize() override;

    /*!
    \brief Evaluate contact

    For a lagrangian strategy this involves heavy modification to the initial kteff and feff.
    Hence, they are in fact build from scratch here.
    The application of modifications to groups of dofs (slave, master, active etc.)
    results in some matrix and vector splitting and a lot of matrix-vector calculation in here!

    */
    void EvaluateContact(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override;

    /*!
    \brief Evaluate frictional contact

    */
    void EvaluateFriction(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override;

    void Update(Teuchos::RCP<const Epetra_Vector> dis) override;


   protected:
    /*!
    \brief Add penalty terms for LTL edge contact

    */
    void add_line_to_lin_contributions(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, bool add_time_integration = true);

    /*!
    \brief Add penalty terms for master contact

    */
    void add_master_contributions(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, bool add_time_integration = true);

    /*!
    \brief Add penalty terms for LTL edge contact

    */
    void add_line_to_lin_contributions_friction(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, bool add_time_integration = true);

    // don't want = operator and cctor
    LagrangeStrategy operator=(const LagrangeStrategy& old) = delete;
    LagrangeStrategy(const LagrangeStrategy& old) = delete;

    // Store Coupling Matrices in case of Poro Lagrange Strategy ... here just ignore!
    virtual void save_coupling_matrices(Teuchos::RCP<CORE::LINALG::SparseMatrix> dhat,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> mhataam,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> invda)
    {
      return;
    }

    std::vector<Teuchos::RCP<CONTACT::Interface>> interface_;

    bool evalForceCalled_;  //< flag for evaluate force call
    bool activesetssconv_;  //< convergence flag for semi-smooth active set search
    bool activesetconv_;    //< convergence flag for fixed-point active set search
    int activesetsteps_;    //< number of fixed-point active set steps in this time step

    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        mhatmatrix_;  //< product of global Mortar matrices inv(D)*M

    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        tmatrix_;  //< global Matrix T containing active node tangents
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        nmatrix_;  //< global Matrix N containing active node normals

    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        smatrix_;  //< global Matrix S containing normal+D+M derivatives
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        smatrixW_;  //< global Matrix S containing W derivatives
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        tderivmatrix_;  //< global Matrix containing tangent derivatives
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        nderivmatrix_;  //< global Matrix containing normal derivatives


    Teuchos::RCP<Epetra_Vector> fs_;                 //< slave side effective forces (needed for LM)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> invd_;  //< inverse of Mortar matrix D (needed for LM)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ksn_;   //< stiffness block K_sn (needed for LM)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ksm_;   //< stiffness block K_sm (needed for LM)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> kss_;   //< stiffness block K_ss (needed for LM)

    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        linslipLM_;  //< global matrix containing derivatives (LM) of slip condition
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        linslipDIS_;  //< global matrix containing derivatives (DIS) of slip condition
    Teuchos::RCP<Epetra_Vector> linslipRHS_;  //< r.h.s vector friction slip nodes
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        linstickLM_;  //< global matrix containing derivatives (LM) of slip condition
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        linstickDIS_;  //< global matrix containing derivatives (DIS) of stick condition
    Teuchos::RCP<Epetra_Vector> linstickRHS_;  //< r.h.s vector for friction stick condition

    Teuchos::RCP<Epetra_Map> zigzagone_;    //< active node set of last active set try
    Teuchos::RCP<Epetra_Map> zigzagtwo_;    //< active node set of second-last active set try
    Teuchos::RCP<Epetra_Map> zigzagthree_;  //< active node set of third-last active set try

    Teuchos::RCP<Epetra_Map> gOldActiveSlaveNodes_;  // active slave nodes from previous newton step

    Teuchos::RCP<Epetra_Map> gOldslipnodes_;  // active slave nodes from previous newton step

    Teuchos::RCP<Epetra_Vector> fLTLOld_;        //< old line to line forces
    Teuchos::RCP<Epetra_Vector> fLTL_;           //< current line to line forces combined
    Teuchos::RCP<Epetra_Vector> fLTLn_;          //< current line to line forces normal
    Teuchos::RCP<Epetra_Vector> fLTLt_;          //< current line to line forces tangent
    Teuchos::RCP<Epetra_Vector> fconservation_;  //< current line to line forces normal

    Teuchos::RCP<Epetra_Vector> nonsmooth_Penalty_force_;  //< penalty forces of non-smooth contact
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        nonsmooth_Penalty_stiff_;  //< tangent to penalty forces of non-smooth contact

  };  // class LagrangeStrategy
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
