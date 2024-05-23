/*----------------------------------------------------------------------*/
/*! \file
\brief wear strategy for finite wear modeling

\level 2

*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              farah 09/13 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_LAGRANGE_STRATEGY_WEAR_HPP
#define FOUR_C_CONTACT_LAGRANGE_STRATEGY_WEAR_HPP

/*----------------------------------------------------------------------*
 | header                                                   farah 09/13 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_contact_lagrange_strategy.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN

namespace WEAR
{
  // forward declarations
  class WearInterface;


  class LagrangeStrategyWear : public CONTACT::LagrangeStrategy
  {
   public:
    /*!
    \brief Standard Constructor

    */
    LagrangeStrategyWear(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interfaces, int dim,
        Teuchos::RCP<const Epetra_Comm> comm, double alphaf, int maxdof);


    /*!
    \brief Condense discr. wear and lm. for frictional contact

    */
    void CondenseWearDiscr(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector>& gact);

    /*!
    \brief Condense lm. for frictional contact with explicit/implicit wear algorithm

    */
    void condense_wear_impl_expl(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector>& gact);

    /*!
    \brief Prepare SaddlePointSystem

    */
    void prepare_saddle_point_system(
        Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff);

    /*!
    \brief Recovery method

    We only recover the Lagrange multipliers here, which had been
    statically condensed during the setup of the global problem!
    Optionally satisfaction or violation of the contact boundary
    conditions can be checked, too.

    */
    void Recover(Teuchos::RCP<Epetra_Vector> disi) override;

    /*!
    \brief Redistribute all contact interfaces in parallel

    We hand in the current global displacement state so that a contact search can be performed and
    set state called.

    The current velocity state is required in case of extedning the ghosting via binning to account
    for relative motion between interfaces.

    \param[in] dis Current displacement state
    \param[in] vel Current velocity state

    \return TRUE if the interface has been redistributed. Return FALSE otherwise.
    */
    bool RedistributeContact(
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel) final;

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
    void build_saddle_point_system(Teuchos::RCP<CORE::LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) override;

    /*!
    \brief Update internal member variables after solving the 2x2 saddle point contact system

    \param sold (out): the displacement dof solution increment (associated with displacement dofs)
    \param mergedsol (in): Epetra_Vector for merged solution vector (containing the new solution
    vector of the full merged linear system)
    */
    void update_displacements_and_l_mincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) override;

    /*!
    \brief Evaluate wear vector

    Evaluates the unweighted wear vector.
    Refer also to the Semesterarbeit of Karl Wichmann 2010

    \warning In case of weighted wear, this is only implemented for dual/Petrov-Galerkin shape
    functions.

    \note This requires to solve for the real wear where we need to invert the slave side mortar
    matrix. Since we restrict ourselves to dual/Petrov-Galerkin shape functions, we can exploit
    the diagonality of the #dmatrix_. Hence, instead of solving we just divide by the diagonal
    elements of #dmatrix_.

    */
    void OutputWear() override;


    /*!
    \brief Perform a write restart

    A write restart is initiated by the contact manager. However, the manager has no
    direct access to the nodal quantities. Hence, a portion of the restart has to be
    performed on the level of the contact algorithm, for short: here's the right place.

    */
    void DoWriteRestart(std::map<std::string, Teuchos::RCP<Epetra_Vector>>& restart_vectors,
        bool forcedrestart = false) const override;

    /*!
    \brief Perform a write restart

    A write restart is initiated by the contact manager. However, the manager has no
    direct access to the nodal quantities. Hence, all the restart action has to be
    performed on the level of the contact algorithm, for short: here's the right place.

    */
    void DoReadRestart(
        IO::DiscretizationReader& reader, Teuchos::RCP<const Epetra_Vector> dis) override;

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
    void update_active_set_semi_smooth(const bool firstStepPredictor = false) override;

    /*!
    \brief Store/Reset nodal wear quantities for pv approach

    */
    void update_wear_discret_iterate(bool store);

    /*!
    \brief Store wear for accumulation due to different pseudo time scales!

    */
    void update_wear_discret_accumulation();

    /*!
    \brief Update wear contact at end of time step

    */
    void Update(Teuchos::RCP<const Epetra_Vector> dis) override;

    /*!
    \brief Store wear data into wear data container

    */
    void store_nodal_quantities(MORTAR::StrategyBase::QuantityType type) override;

    /*!
    \brief Return vector of wear (t_n+1) - D^-1 \times weighted wear!

    */
    Teuchos::RCP<Epetra_Vector> ContactWear() override { return wearoutput_; }  // for slave side
    Teuchos::RCP<const Epetra_Vector> ContactWear() const override
    {
      return wearoutput_;
    }                                                                    // for slave side
    Teuchos::RCP<Epetra_Vector> ContactWear2() { return wearoutput2_; }  // for master side

    /*!
    \brief Return wear interfaces

    */
    std::vector<Teuchos::RCP<WEAR::WearInterface>> WearInterfaces() { return interface_; }

    /*!
    \brief Return master map for both sided wear (slip), mapped from slave side

    */
    Teuchos::RCP<const Epetra_Map> MasterSlipNodes() const override { return gmslipnodes_; };

    /*!
    \brief Return master map for both sided wear (active), mapped from slave side

    */
    Teuchos::RCP<const Epetra_Map> MasterActiveNodes() const override { return gmactivenodes_; };

    /*!
     \brief Return discrete wear vector (t_n+1)

     */
    Teuchos::RCP<Epetra_Vector> WearVar() { return w_; }

    /*!
     \brief Return discrete wear vector (t_n+1) Master

     */
    Teuchos::RCP<Epetra_Vector> WearVarM() { return wm_; }

    /*!
     \brief Return wear rhs vector (only in saddle-point formulation

     */
    Teuchos::RCP<Epetra_Vector> WearRhs() override { return wearrhs_; }

    /*!
     \brief Return wear-master rhs vector (only in saddle-point formulation

     */
    Teuchos::RCP<Epetra_Vector> WearMRhs() override { return wearmrhs_; }

    /*!
     \brief Returns increment of W solution vector in SaddlePointSolve routine

     */
    Teuchos::RCP<Epetra_Vector> WSolveIncr() override { return wincr_; }

    /*!
     \brief Returns increment of W-master solution vector in SaddlePointSolve routine

     */
    Teuchos::RCP<Epetra_Vector> WMSolveIncr() override { return wmincr_; }

    /*!
     \brief Return global both sided wear status

     */
    bool WearBothDiscrete() const override { return wbothpv_; }

    /*!
     \brief Return global wear status

     */
    bool WeightedWear() const override { return weightedwear_; }

   private:
    /*!
    \brief Evaluate frictional contact

    */
    void EvaluateFriction(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override;

    /*!
    \brief Initialize and evaluate Mortar stuff for the next Newton step

    This method first checks if we are dealing with self contact and updates
    the interface slave and master sets if so. Then it resets the global
    Mortar matrices D and M and the global gap vector g accordingly.

    The nodal quantites computed in InitEvalInterface() are then assembled
    to global matrices and vectors respectively. No setup of the global system
    is to be done here yet, so there is no need to pass in the effective
    stiffness K or the effective load vector f.

    */
    void InitMortar() override;
    void AssembleMortar() override;

    /*!
    \brief Initialize general contact variables for next Newton step

    For a lagrangian strategy this includes the global normal / tangent matrices N and T,
    the global derivative matrices S and P and Tresca friction matrix L + vector r.

    */
    void Initialize() override;

    /*!
    \brief Setup this strategy object (maps, vectors, etc.)

    All global maps and vectors are initialized by collecting
    the necessary information from all interfaces. In the case
    of a parallel redistribution, this method is called again
    to re-setup the above mentioned quantities. In this case
    we set the input parameter redistributed=TRUE. Moreover,
    when called for the first time (in the constructor) this
    method is given the input parameter init=TRUE to account
    for initialization of the active set.
      */
    void Setup(bool redistributed, bool init) override;

    /*!
    \brief Setup this strategy object (maps, vectors, etc.)

    All wear specific maps here
    */
    void SetupWear(bool redistributed, bool init);

   private:
    // don't want = operator and cctor
    LagrangeStrategyWear operator=(const LagrangeStrategyWear& old) = delete;
    LagrangeStrategyWear(const LagrangeStrategyWear& old) = delete;

    std::vector<Teuchos::RCP<WEAR::WearInterface>> interface_;

    // basic data
    bool weightedwear_;              // flag for contact with wear (is) --> weighted wear
    bool wbothpv_;                   // flag for both sided wear disrete
    Teuchos::RCP<Epetra_Vector> w_;  // current vector of pv wear at t_n+1 (slave)
    Teuchos::RCP<Epetra_Vector>
        wincr_;  // Wear variables vector increment within SaddlePointSolve (this is NOT the
                 // increment of w_ between t_{n+1} and t_{n}!)
    Teuchos::RCP<Epetra_Vector> wearrhs_;

    Teuchos::RCP<Epetra_Vector> wm_;  // current vector of pv wear at t_n+1 (master)
    Teuchos::RCP<Epetra_Vector>
        wmincr_;  // Wear variables vector increment within SaddlePointSolve (this is NOT the
                  // increment of w_ between t_{n+1} and t_{n}!)
    Teuchos::RCP<Epetra_Vector> wearmrhs_;

    // implicit wear algorithm
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        wlinmatrix_;  // global Matrix Wg containing wear-lm derivatives
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        wlinmatrixsl_;  // global Matrix Wsl containing wear-lm slip derivatives
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        wlinmatrixst_;  // global Matrix Wst containing wear-lm stick derivatives

    // both-sided wear weak dirich cond
    Teuchos::RCP<CORE::LINALG::SparseMatrix> d2matrix_;  // global Mortar matrix D2

    Teuchos::RCP<Epetra_Map>
        gminvolvednodes_;  // global involved master node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>
        gminvolveddofs_;  // global involved master dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>
        gslipn_;  // global row map of matrix N for slip dofs (of all interfaces)
    Teuchos::RCP<Epetra_Map>
        gwinact_;  // global row map of matrix N for slip dofs (of all interfaces)
    Teuchos::RCP<Epetra_Map>
        gmslipn_;  // global row map of matrix N for slip dofs (of all interfaces)
    Teuchos::RCP<Epetra_Map>
        gwminact_;  // global row map of matrix N for slip dofs (of all interfaces)

    Teuchos::RCP<Epetra_Map>
        gwmdofrowmap_;  // global master wear dof row map (of all interfaces) -active
    Teuchos::RCP<Epetra_Map>
        gwdofrowmap_;  // global slave wear dof row map (of all interfaces) -active
    Teuchos::RCP<Epetra_Map> gsdofnrowmap_;    // global slave wear dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gmdofnrowmap_;    // global master wear dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> galldofnrowmap_;  // global master wear dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gwalldofrowmap_;  // all
    Teuchos::RCP<Epetra_Map> gmslipnodes_;     // global master slip nodes
    Teuchos::RCP<Epetra_Map> gmactivenodes_;   // global master active nodes

    Teuchos::RCP<Epetra_Vector> wearoutput_;   // vector of unweighted wear at t_n+1  -- slave
    Teuchos::RCP<Epetra_Vector> wearoutput2_;  // vector of unweighted wear at t_n+1  -- master
    Teuchos::RCP<Epetra_Vector> wearvector_;   // global weighted wear vector w

    int maxdofwear_;  // highest dof number in problem discretization

    bool wearimpl_;        // weartype: implicit
    bool wearprimvar_;     // bool for wear with own discretization
    bool wearbothpv_;      // bool for both-sided discrete wear
    bool weartimescales_;  // bool for different time scales
    bool sswear_;          // bool steady state wear

    // discrete wear algorithm (SLAVE)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> twmatrix_;  // global Mortar wear matrix T
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ematrix_;   // global Mortar wear matrix E
    Teuchos::RCP<CORE::LINALG::SparseMatrix> eref_;      // global Mortar wear matrix E
    Teuchos::RCP<CORE::LINALG::SparseMatrix> lintdis_;   // Lin T w.r.t. displ: Lin(T*n*lm)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> lintlm_;    // Lin T w.r.t. lm: (T*n)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> linedis_;   // Lin E w.r.t. displ: Lin(E*w)
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        linslip_w_;  // global matrix containing derivatives (LM) of slip condition
    Teuchos::RCP<Epetra_Vector> inactive_wear_rhs_;  // inactive wear rhs: -w_i
    Teuchos::RCP<Epetra_Vector> wear_cond_rhs_;      // rhs wear condition: -E*w_i + k*T*n*lm_i

    // discrete wear algorithm (MASTER)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> twmatrix_m_;  // global Mortar wear matrix T
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ematrix_m_;   // global Mortar wear matrix E
    Teuchos::RCP<CORE::LINALG::SparseMatrix> lintdis_m_;   // Lin T w.r.t. displ: Lin(T*n*lm)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> lintlm_m_;    // Lin T w.r.t. lm: (T*n)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> linedis_m_;   // Lin E w.r.t. displ: Lin(E*w)
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        linslip_wm_;  // global matrix containing derivatives (LM) of slip condition
    Teuchos::RCP<Epetra_FEVector> inactive_wear_rhs_m_;  // inactive wear rhs: -w_i
    Teuchos::RCP<Epetra_FEVector> wear_cond_rhs_m_;      // rhs wear condition: -E*w_i + k*T*n*lm_i

    // matrix blocks for recovering
    Teuchos::RCP<CORE::LINALG::SparseMatrix> dnblock_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> dmblock_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> diblock_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> dablock_;
    Teuchos::RCP<Epetra_Vector> fw_;

    Teuchos::RCP<Epetra_Map> gidofs_;

  };  // class

}  // namespace WEAR

FOUR_C_NAMESPACE_CLOSE

#endif
