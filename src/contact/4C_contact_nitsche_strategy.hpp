/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NITSCHE_STRATEGY_HPP
#define FOUR_C_CONTACT_NITSCHE_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_contact_abstract_strategy.hpp"
#include "4C_contact_utils.hpp"

#include <Epetra_FEVector.h>

#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  /*!
   \brief Contact solving strategy with Nitsche's method.

   This is a specialization of the abstract contact algorithm as defined in AbstractStrategy.
   For a more general documentation of the involved functions refer to CONTACT::AbstractStrategy.

   */
  class NitscheStrategy : public AbstractStrategy
  {
   public:
    //! Standard constructor
    NitscheStrategy(const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        const Teuchos::RCP<Epetra_Comm>& comm, double alphaf, int maxdof)
        : AbstractStrategy(Teuchos::rcp(new CONTACT::AbstractStratDataContainer()), dof_row_map,
              NodeRowMap, params, dim, comm, alphaf, maxdof),
          interface_(std::move(interface)),
          curr_state_eval_(false)
    { /* empty */
    }

    //! Shared data constructor
    NitscheStrategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        const Teuchos::RCP<const Epetra_Comm>& comm, double alphaf, int maxdof)
        : AbstractStrategy(data_ptr, dof_row_map, NodeRowMap, params, dim, comm, alphaf, maxdof),
          interface_(std::move(interface)),
          curr_state_eval_(false)
    { /* empty */
    }

    // don't want = operator and cctor
    NitscheStrategy operator=(const NitscheStrategy& old) = delete;
    NitscheStrategy(const NitscheStrategy& old) = delete;

    void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<CORE::LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, int step,
        int iter, bool predictor) override;

    void DoReadRestart(CORE::IO::DiscretizationReader& reader,
        Teuchos::RCP<const Epetra_Vector> dis,
        Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr) override;

    bool IsSaddlePointSystem() const override { return false; }

    bool IsCondensedSystem() const override { return false; }

    /*!
     * @brief Integrate all contact interfaces
     *
     * @note this method is called from the new structural time integration
     *
     * @param[in] cparams  contact data container
     */
    virtual void Integrate(const CONTACT::ParamsInterface& cparams);

    Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum CONTACT::VecBlockType& bt) const override;

    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetMatrixBlockPtr(const enum CONTACT::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams) const override;

    /*! \brief Setup this strategy object (maps, vectors, etc.)

     derived from contact abstract strategy.
     The Nitsche strategy does not have
      */
    void Setup(bool redistributed, bool init) override;

    virtual void update_trace_ineq_etimates();

    /*! \brief Get dirichlet B.C. status and store into Nodes

     This is called once at the beginning of the simulation
     to set the D.B.C. status in each CNode.

     \param dbcmaps (in): MapExtractor carrying global dbc map */
    void store_dirichlet_status(Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmaps) override{
        /* we don't care about dirichlet for now */
    };
    void Update(Teuchos::RCP<const Epetra_Vector> dis) override;
    void evaluate_reference_state() override;
    void DoWriteRestart(std::map<std::string, Teuchos::RCP<Epetra_Vector>>& restart_vectors,
        bool forcedrestart) const override{
        /* nothing stored in nitsche strategy that would need to be written */
    };
    void compute_contact_stresses() final{/* nothing stress output in nitsche strategy yet */};
    virtual void reconnect_parent_elements();
    void set_state(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) override;

    /*!
     * @brief  Set the parent state
     *
     * @param[in] statename  name of state to be set
     * @param[in] vec        corresponding state vector
     */
    virtual void SetParentState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec);

    Teuchos::RCP<const Epetra_Vector> GetLagrMultN(const bool& redist) const override
    {
      return Teuchos::null;
    }
    Teuchos::RCP<const Epetra_Vector> GetLagrMultNp(const bool& redist) const override
    {
      return Teuchos::null;
    }
    Teuchos::RCP<Epetra_Vector> LagrMultOld() override { return Teuchos::null; }
    Teuchos::RCP<const Epetra_Map> LMDoFRowMapPtr(const bool& redist) const override
    {
      return Teuchos::null;
    }
    // All these functions only have functionality in Lagrange contact simulations,
    // thus they are defined empty here in the case of Penalty contact.

    //! Get the active node row map of the previous Newton step
    Teuchos::RCP<const Epetra_Map> get_old_active_row_nodes() const override
    {
      return Teuchos::null;
    };
    Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const override { return Teuchos::null; };
    bool IsNitsche() const override { return true; }
    void PrintActiveSet() const override{};
    bool active_set_semi_smooth_converged() const override { return true; }
    bool ActiveSetConverged() override { return true; }
    int ActiveSetSteps() override { return 0; }
    void ResetActiveSet() override {}
    void Recover(Teuchos::RCP<Epetra_Vector> disi) override {}
    void build_saddle_point_system(Teuchos::RCP<CORE::LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) override
    {
      FOUR_C_THROW(
          "Nitsche does not have Lagrange multiplier DOFs. So, saddle point system makes no sense "
          "here.");
    }
    void update_displacements_and_l_mincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) override
    {
      FOUR_C_THROW(
          "Nitsche does not have Lagrange multiplier DOFs. So, saddle point system makes no sense "
          "here.");
    }
    void EvalConstrRHS() override {}
    void UpdateActiveSet() override {}
    void update_active_set_semi_smooth(const bool firstStepPredictor) override {}
    void evaluate_rel_mov_predict() override {}
    void ModifyPenalty() override {}
    void update_uzawa_augmented_lagrange() override {}
    void update_constraint_norm(int uzawaiter) override {}
    void Initialize() override{};
    void EvaluateContact(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override
    {
      FOUR_C_THROW("not supported in this strategy");
    }
    void EvaluateFriction(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override
    {
      FOUR_C_THROW("not supported in this strategy");
    }
    void InitializeUzawa(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override
    {
    }
    void ResetPenalty() override {}
    void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) override {}
    double InitialPenalty() override { return 0.0; }
    double ConstraintNorm() const override { return 0.0; }
    bool IsPenalty() const override { return false; }

   protected:
    std::vector<Teuchos::RCP<CONTACT::Interface>>& interfaces() override { return interface_; }

    const std::vector<Teuchos::RCP<CONTACT::Interface>>& interfaces() const override
    {
      return interface_;
    }

    void eval_force(CONTACT::ParamsInterface& cparams) override;

    void eval_force_stiff(CONTACT::ParamsInterface& cparams) override;

    void Reset(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dispnp,
        const Epetra_Vector& xnew) override;

    void run_post_compute_x(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
        const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

    /*!
     * @brief Fill RHS vector of vector block type
     *
     * @param[in] bt  vector block type
     * @return the filled RHS vector of given vector block type
     */
    virtual Teuchos::RCP<Epetra_FEVector> create_rhs_block_ptr(
        const enum CONTACT::VecBlockType& bt) const;

    /*!
     * @brief  Create an appropriate vector for the RHS
     *
     * @param[in] bt  block type
     * @return  vector for given vector block type
     */
    virtual Teuchos::RCP<Epetra_FEVector> setup_rhs_block_vec(
        const enum CONTACT::VecBlockType& bt) const;

    /*!
     * @brief Create appropriate matrix block
     *
     * @param[in] bt  matrix block type
     * @return matrix block for given matrix block type
     */
    virtual Teuchos::RCP<CORE::LINALG::SparseMatrix> setup_matrix_block_ptr(
        const enum MatBlockType& bt);

    /*!
     * @brief Complete the matrix block with correct maps
     *
     * @param[in] bt      matrix block type
     * @param[in,out] kc  matrix block of given matrix block type that has to be completed
     */
    virtual void complete_matrix_block_ptr(
        const enum MatBlockType& bt, Teuchos::RCP<CORE::LINALG::SparseMatrix> kc);

    /*!
     * @brief Fill block matrix of given matrix block type
     *
     * @param[in] bt  matrix block type
     * @return the filled block matrix of given matrix block type
     */
    virtual Teuchos::RCP<CORE::LINALG::SparseMatrix> create_matrix_block_ptr(
        const enum MatBlockType& bt);

    std::vector<Teuchos::RCP<CONTACT::Interface>> interface_;

    Teuchos::RCP<Epetra_Vector> curr_state_;
    bool curr_state_eval_;

    Teuchos::RCP<Epetra_FEVector> fc_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> kc_;
  };
}  // namespace CONTACT
FOUR_C_NAMESPACE_CLOSE

#endif
