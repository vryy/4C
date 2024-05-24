/*---------------------------------------------------------------------*/
/*! \file
\brief Main abstract class for meshtying solution strategies

\level 2


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_MESHTYING_ABSTRACT_STRATEGY_HPP
#define FOUR_C_CONTACT_MESHTYING_ABSTRACT_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_contact_utils.hpp"
#include "4C_inpar_mortar.hpp"
#include "4C_mortar_strategy_base.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace MORTAR
{
  class Interface;
}

namespace CORE::LINALG
{
  class SparseMatrix;
}
namespace NOX
{
  namespace NLN
  {
    class Group;
  }  // namespace NLN
}  // namespace NOX

namespace CONTACT
{
  // forward declarations
  class MtNoxInterface;

  /*! \brief Main abstract class for meshtying solution strategies

  This is the templating abstract class for all meshtying solution algorithms.
  Every solution algorithm has to fit into the set of functions and calls defined herein
  and has to be specified in a corresponding subclass defining the concrete algorithmic steps.

  This class it itself derived from the MORTAR::StrategyBase class, which is an even
  more abstract framework for any solution strategies involving mortar coupling.

  */
  class MtAbstractStrategy : public MORTAR::StrategyBase
  {
   public:
    /*!
    \brief Standard Constructor

    Creates the strategy object and initializes all global variables, including
    all necessary Epetra_Maps and global vector and matrix quantities.

    \param[in] dof_row_map Dof row map of underlying problem
    \param[in] NodeRowMap Node row map of underlying problem
    \param[in] params List of contact/parameters
    \param[in] interface All contact interface objects
    \param[in] spatialDim Spatial dimension of the problem
    \param[in] comm Communicator
    \param[in] alphaf Mid-point for Generalized-alpha time integration
    \param[in] maxdof Highest DOF number in global problem
    */
    MtAbstractStrategy(const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<MORTAR::Interface>> interface,
        const int spatialDim, const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf,
        const int maxdof);



    //! @name Access methods

    /*!
    \brief Return Lagrange multiplier vector (t_n+1)

    */
    Teuchos::RCP<Epetra_Vector> LagrMult() override { return z_; }

    /*!
    \brief Return old Lagrange multiplier vector (t_n)

    */
    Teuchos::RCP<Epetra_Vector> LagrMultOld() override { return zold_; }

    /*!
    \brief Return Lagrange multiplier vector from last Uzawa step

    */
    Teuchos::RCP<Epetra_Vector> LagrMultUzawa() { return zuzawa_; }

    /*!
    \brief Return constraint rhs vector (only in saddle-point formulation

    */
    Teuchos::RCP<Epetra_Vector> ConstrRhs() override { return constrrhs_; }

    /*!
    \brief Returns increment of LagrangeMultiplier solution vector in SaddlePointSolve routine

    */
    Teuchos::RCP<Epetra_Vector> LagrMultSolveIncr() override { return zincr_; }

    /*!
    \brief Gather maps needed for contact/meshtying specific multigrid preconditioners

    @param MasterDofMap Dof row map of master interface
    @param SlaveDofMap Dof row map of slave interface
    @param InnerDofMap Dof row map of interior volume
    @param ActiveDofMap Dof row map of active slave contact interface
    */
    void collect_maps_for_preconditioner(Teuchos::RCP<Epetra_Map>& MasterDofMap,
        Teuchos::RCP<Epetra_Map>& SlaveDofMap, Teuchos::RCP<Epetra_Map>& InnerDofMap,
        Teuchos::RCP<Epetra_Map>& ActiveDofMap) const override;

    /*!
    \brief Return mortar matrix D

    */
    Teuchos::RCP<CORE::LINALG::SparseMatrix> DMatrix() override { return dmatrix_; }

    /*!
    \brief Return mortar matrix M

    */
    Teuchos::RCP<CORE::LINALG::SparseMatrix> MMatrix() override { return mmatrix_; }

    /*!
    \brief Get dual quadratic 3d slave element flag

    Returns TRUE if at least one higher-order 3d slave element with
    dual Lagrange multiplier shape functions in any interface.

    */
    virtual const bool& Dualquadslavetrafo() const { return dualquadslavetrafo_; };

    /*!
    \brief Return parallel redistribution status (yes or no)

    */
    bool ParRedist() const
    {
      INPAR::MORTAR::ParallelRedist partype =
          Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(
              Params().sublist("PARALLEL REDISTRIBUTION"), "PARALLEL_REDIST");
      if (partype != INPAR::MORTAR::ParallelRedist::redist_none)
        return true;
      else
        return false;
    }

    //@}

    //! @name Evaluation methods

    /*!
    \brief Redistribute all meshtying interfaces in parallel

    Here, we call each interface to perform redistribution for each interface individually. Since
    this changes maps and interface discretizations, we have to fill_complete() all interface
    discretizations and re-setup the strategy object afterwards by calling Setup(bool).

    If parallel redistribution is disabled in the input file or if this is a serial computation,
    i.e. only one MPI rank, then we just print the current parallel distribution to the screen, but
    do not change it.

    \pre Meshtying interface discretizations are distributed to multiple processors, but maybe not
    in an optimal way, i.e. with sub-optimal load balancing.

    \post If desired by the user, all meshtying interface discretizations have been re-distributed,
    such that load balancing among processors is closer to optimal.
    */
    void redistribute_meshtying() final;

    /*!
    \brief Global evaluation method called from time integrator

    */
    void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<CORE::LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f,
        const int step, const int iter, bool predictor = false) override;

    /*! \brief Reset call at the beginning of the ApplyForce(), ApplyStiff() and ApplyForceStiff()
     * [derived]
     *
     *  \date 02/2016
     *  \author hiermeier */
    virtual void Reset(const Epetra_Vector& dis)
    {
      FOUR_C_THROW("Not yet considered for meshtying!");
    };

    /*! \brief Global evaluation method called from STR::MODELEVALUATOR::Contact class [derived]
     *
     *  Evaluation of the right-hand-side only. Necessary and meaningful for line search strategies
     *  for example.
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual bool ApplyForce()
    {
      FOUR_C_THROW("Not yet considered for msht!");
      return false;
    };

    /*! \brief Global evaluation method called from STR::MODELEVALUATOR::Contact class [derived]
     *
     *  Evaluation of the mesh-tying right-hand-side and the mesh-tying jacobian. We call this
     * method also, when we are only interested in the jacobian, since the created overhead is
     * negligible.
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual bool ApplyForceStiff()
    {
      FOUR_C_THROW("Not yet considered for msht!");
      return false;
    };

    /*!
    \brief Set current deformation state

    All interfaces are called to set the current deformation state.

    \param statename (in): std::string defining which quantity to set (only "displacement"
    applicable) \param vec (in): current global state of the quantity defined by statename

    */
    void set_state(const enum MORTAR::StateType& statetype, const Epetra_Vector& vec) override;

    /*!
    \brief Do mortar coupling in reference configuration

    Only do this ONCE for meshtying upon initialization!
    This method calls Initialize() on all contact interfaces, which
    resets all kind of nodal quantities. It then calls Evaluate() on
    all meshtying interfaces, which does all the geometric coupling stuff.
    Concretely, this is an evaluation of all involved quantities at nodal
    level. It includes the nodal normal calculations, search, projection
    and overlap detection and integration of the Mortar terms D and M.

    Then - on global level - it resets the Mortar matrices D and M accordingly.
    The nodal quantities computed before are assembled to global matrices. No
    setup of the global system is to be done here yet, so there is no need to
    pass in the effective stiffness K or the effective load vector f.

    Note: Only quantities common to all subsequent solving strategies (Lagrange,
    Penalty) are computed here. In case they need additional mortar variables,
    use the overloaded function call in the derived class and refer back to this function.

    */
    void MortarCoupling(const Teuchos::RCP<const Epetra_Vector>& dis) override;

    //@}

    //! @name Quantity control methods

    /*!
    \brief Get some nodal quantity globally and store into MORTAR::Node(s)

    The enum input parameter defines, which quantity has to be updated.
    Currently, the possibilities "lmold", "lmcurrent", "lmupdate" and
    "lmuzawa" exist. Note that "lmold" means the converged value LM_n
    of the last time / load step, whereas "lmcurrent" adresses the current
    (not necessarily converged) value of the LM_n+1. "lmupdate" is a special
    option called only in Recover() after the update of the Lagr. multipliers.
    It basically does the same as "lmcurrent", but also checks for D.B.C.
    problems. Finally, "lmuzawa" addresses the LM update within an Uzawa
    augmented Lagrangian scheme.

    \param type (in): enum defining which quantity to store into MORTAR::Node(s)

    */
    void store_nodal_quantities(MORTAR::StrategyBase::QuantityType type) override;

    /*!
    \brief Get dirichlet B.C. status and store into MORTAR::Node(s)

    This is called once at the beginning of the simulation
    to set the D.B.C. status in each MORTAR::Node.

    \param dbcmaps (in): MapExtractor carrying global dbc map

    */
    void store_dirichlet_status(Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmaps) override;

    /*!
    \brief Update meshtying at end of time step

    \param dis (in):  current displacements (-> old displacements)

    */
    void Update(Teuchos::RCP<const Epetra_Vector> dis) override;

    /*!
    \brief Perform a write restart

    A write restart is initiated by the contact manager. However, the manager has no
    direct access to the nodal quantities. Different from writing a restart step, now
    all the restart action has to be performed on the level of the meshtying algorithm,
    for short: here's the right place.

    */
    void DoReadRestart(
        IO::DiscretizationReader& reader, Teuchos::RCP<const Epetra_Vector> dis) override;

    //@}

    //! @name Output

    /*!
    \brief Compute interface forces and moments

    Compute current interface forces and moments at n+1-alphaf using current
    Lagrange multiplier values and current Mortar matrices D and M at n+1. When
    doing dynamics with alpha_f > 0, this also uses the old LM and Mortar
    matrices of the last converged time / load step n (TR-like interpolation).

    \param fresm (in): residual / force vector at state n+1 of current Newton step
    \param output (in): flag indicating whether force output shall be written

    */
    void InterfaceForces(bool output = false) override;

    /*!
    \brief Print interfaces

    \param[in] os Output stream used for printing
    */
    void Print(std::ostream& os) const override;

    /*!
    \brief Print current active set to screen for debugging purposes

    */
    void PrintActiveSet() const override;

    /*!
    \brief Write results for visualization separately for each meshtying/contact interface

    Call each interface, such that each interface can handle its own output of results.

    \param[in] outputParams Parameter list with stuff required by interfaces to write output
    */
    void postprocess_quantities_per_interface(
        Teuchos::RCP<Teuchos::ParameterList> outputParams) final;

    //! @}

    //! @name Debugging methods
    //! @{

    /*!
    \brief Visualize contact stuff with gmsh

    \param step (in): current time step index
    \param iter (in): current nonlinear iteration index

    */
    void VisualizeGmsh(const int step, const int iter) override;

    //! @}

    //! @name Preconditioner methods
    //! @{

    /*! Derived method
     *
     * (see NOX::NLN::CONSTRAINT::Interface::Preconditioner for more information) */
    bool IsSaddlePointSystem() const override;

    /*! \brief Derived method
     *
     * (see NOX::NLN::CONSTRAINT::Interface::Preconditioner for more information) */
    bool IsCondensedSystem() const override;

    /*! Fill the maps vector for the linear solver preconditioner
     *
     * The following order is defined:
     * (0) masterDofMap
     * (1) slaveDofMap
     * (2) innerDofMap
     * (3) activeDofMap
     *
     * \author hiermeier */
    void fill_maps_for_preconditioner(std::vector<Teuchos::RCP<Epetra_Map>>& maps) const override;

    //! compute the preconditioner operator
    bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M,
        Teuchos::ParameterList* precParams = nullptr) override;

    //! @}

    /*! @name Purely virtual functions
     *
     * All these functions are defined in one or more specific derived classes,
     * i.e CONTACT::MeshtyingLagrangeStrategy or CONTACT::MeshtyingPenaltyStrategy.
     * As the base class MORTAR::StrategyBase is always called from the control routine
     * (time integrator), these functions need to be defined purely virtual here.
     */

    double ConstraintNorm() const override = 0;
    void EvaluateMeshtying(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) override = 0;
    void InitializeUzawa(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override = 0;
    double InitialPenalty() override = 0;
    void Recover(Teuchos::RCP<Epetra_Vector> disi) override = 0;
    void ResetPenalty() override = 0;
    void ModifyPenalty() override = 0;
    void build_saddle_point_system(Teuchos::RCP<CORE::LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) override = 0;
    void update_displacements_and_l_mincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) override = 0;
    void update_uzawa_augmented_lagrange() override = 0;
    void update_constraint_norm(int uzawaiter = 0) override = 0;

    //! @}

    /*! @name Empty functions (contact)
     *
     * All these functions only have functionality in contact simulations, thus they
     * are defined as empty here in the case of meshtying. They can be called from the
     * control routine (time integrator), whenever you like.
     */

    bool ActiveSetConverged() override { return true; }
    bool active_set_semi_smooth_converged() const override { return true; }
    bool Friction() const override { return false; }
    bool WearBothDiscrete() const override { return false; }
    bool IsInContact() const override { return true; }
    bool WasInContact() const override { return true; }
    bool was_in_contact_last_time_step() const override { return true; }
    Teuchos::RCP<Epetra_Vector> ContactNorStress() override { return Teuchos::null; }
    Teuchos::RCP<Epetra_Vector> ContactTanStress() override { return Teuchos::null; }
    Teuchos::RCP<Epetra_Vector> ContactNorForce() override { return Teuchos::null; }
    Teuchos::RCP<Epetra_Vector> ContactTanForce() override { return Teuchos::null; }
    void AssembleMortar() override {}
    void DoWriteRestart(std::map<std::string, Teuchos::RCP<Epetra_Vector>>& restart_vectors,
        bool forcedrestart = false) const override
    {
    }
    void InitEvalInterface() override {}
    void InitMortar() override {}
    void Initialize() override {}
    double Inttime() override { return inttime_; };
    void Inttime_init() override { inttime_ = 0.0; };
    int NumberOfActiveNodes() const override { return 0; }
    int NumberOfSlipNodes() const override { return 0; }
    void compute_contact_stresses() final{};
    void AugForces(Epetra_Vector& augfs_lm, Epetra_Vector& augfs_g, Epetra_Vector& augfm_lm,
        Epetra_Vector& augfm_g){};
    bool RedistributeContact(
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel) final
    {
      return false;
    }
    void ResetActiveSet() override {}
    void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) override {}
    void UpdateActiveSet() override {}
    void update_active_set_semi_smooth(const bool firstStepPredictor = false) override {}
    Teuchos::RCP<CORE::LINALG::SparseMatrix> EvaluateNormals(
        Teuchos::RCP<Epetra_Vector> dis) override
    {
      return Teuchos::null;
    }
    void evaluate_reference_state() override {}
    void EvaluateRelMov() override {}
    void evaluate_rel_mov_predict() override {}
    Teuchos::RCP<Epetra_Map> SlaveRowNodes() override { return gsnoderowmap_; }
    Teuchos::RCP<Epetra_Map> ActiveRowNodes() override { return Teuchos::null; }
    Teuchos::RCP<Epetra_Map> ActiveRowDofs() override { return Teuchos::null; }
    Teuchos::RCP<Epetra_Map> not_re_dist_slave_row_dofs() override { return pgsdofrowmap_; }
    Teuchos::RCP<Epetra_Map> not_re_dist_master_row_dofs() override { return pgmdofrowmap_; }
    Teuchos::RCP<Epetra_Map> SlipRowNodes() override { return Teuchos::null; }

    //! @}

    //! @name New time integration
    //! @{

    //! Return the NOX::NLN::CONSTRAINT::Interface::Required member object
    const Teuchos::RCP<CONTACT::MtNoxInterface>& NoxInterfacePtr() { return noxinterface_ptr_; };

    /*! \brief Return the desired right-hand-side block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired vector block type, e.g. displ, constraint, ...
     *
     *  \date 05/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum CONTACT::VecBlockType& bt) const = 0;

    /*! \brief Return the desired matrix block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired matrix block type, e.g. displ_displ, displ_lm, ...
     *
     *  \date 05/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<CORE::LINALG::SparseMatrix> GetMatrixBlockPtr(
        const enum CONTACT::MatBlockType& bt) const = 0;

    /*! \brief Return the Lagrange multiplier dof row map
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> LMDoFRowMapPtr(const bool& redist) const
    {
      if ((not redist) and ParRedist()) return pglmdofrowmap_;

      return glmdofrowmap_;
    };

    //! Modify system before linear solve
    virtual void run_pre_apply_jacobian_inverse(
        Teuchos::RCP<CORE::LINALG::SparseMatrix> kteff, Epetra_Vector& rhs)
    { /* do nothing */
    }

    //! modify result after linear solve
    virtual void run_post_apply_jacobian_inverse(Epetra_Vector& result)
    { /* do nothing */
    }

    //! evaluate force terms
    virtual bool evaluate_force(const Teuchos::RCP<const Epetra_Vector> dis) = 0;

    //! evaluate stiffness terms
    virtual bool evaluate_stiff(const Teuchos::RCP<const Epetra_Vector> dis) = 0;

    //! evaluate force and stiffness terms
    virtual bool evaluate_force_stiff(const Teuchos::RCP<const Epetra_Vector> dis) = 0;

    //! after applying Newton increment
    virtual void RunPostComputeX(
        const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew){};

    /*! \brief Get the correct RHS for convergence check
     *
     * \todo Is this really about the right-hand side vector or the residual?
     *
     * @param[in/out] rhs Right-hand side vector
     */
    virtual void remove_condensed_contributions_from_rhs(Epetra_Vector& rhs) const {};

    //!@}

   protected:
    /*!
    \brief Assemble global coordinate vector

    \param sidename (in): std::string indicating slave or master side
    \param ref (in): boolean indicating evaluation in reference configuration
    \param vec (in/out)):  empty global vetcor to be assembled to

    */
    void AssembleCoords(const std::string& sidename, bool ref, Teuchos::RCP<Epetra_Vector> vec);

    /*!
    \brief Do mesh initialization for rotational invariance

    Only do this ONCE for meshtying upon initialization!
    This method relocates the slave nodes such that the meshtying constraint
    is satisfied in the reference configuration, which is a prerequisite for
    ensuring both rotational invariance and absence of initial stresses at the
    same time. Basically the constraint equation needs to be solved for this,
    which is specific to the applied solving strategy (dual Lagrange or Penalty).
    In the dual LM, matrix D is diagonal, thus its inversion is trivial and no
    linear system needs to be solved. In the penalty case, matrix D is not diagonal
    and we apply a default CORE::LINALG::Solver to solve for the modified slave positions.
    Thus, this linear system solve is done in the derived method FIRST and then
    we refer back to this base class function.

    \param Xslavemod (in): modified slave reference configuration

    */
    virtual void MeshInitialization(Teuchos::RCP<Epetra_Vector> Xslavemod);

   private:
    /*!
    \brief Evaluate contact

    This is just a tiny control routine, deciding which Evaluate-routine
    of those listed below is to be called (based on input-file information)
    Note that into ALL derived Evaluate() routines, a REFERENCE to the pointer
    on the effective stiffness matrix is handed in. This way, after building the
    new effective stiffness matrix with contact, we can simply let the pointer
    kteff point onto the new object. The same is true for the effective force
    vector feff. Be careful: kteff is of type Teuchos::RCP<CORE::LINALG::SparseOperator>&.

    \param kteff (in/out): effective stiffness matrix (without -> with contact)
    \param feff (in/out): effective residual / force vector (without -> with contact)

    */
    void Evaluate(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) override;

    /*!
    \brief Restrict slave boundary to actual meshtying zone

    Only do this ONCE for meshtying upon initialization!
    This method first detects for each interface the actually tied part
    of the slave surface (i.e. the nodes that carry a D/M contribution).
    Then all slave maps on interface level and on global level are
    re-initialized and re-setup according to the the above defined
    actual slave meshtying zone. This is necessary for problems in which
    the slave surface does not fully project onto the master surface
    and thus the actual meshtying zone cannot be defined within the
    input file. Thus, it is computed here.

    */
    void restrict_meshtying_zone() override;

    /*!
    \brief Setup this strategy object (maps, vectors, etc.)

    All global maps and vectors are initialized by collecting
    the necessary information from all interfaces. In the case
    of a parallel redistribution, this method is called again
    to re-setup the above mentioned quantities. In this case
    the input parameter is set to TRUE.

    */
    void Setup(bool redistributed);

   protected:
    // don't want cctor (= operator impossible anyway for abstract class)
    MtAbstractStrategy(const MtAbstractStrategy& old) = delete;

    //! Vector with all meshtying interfaces
    std::vector<Teuchos::RCP<MORTAR::Interface>> interface_;

    //! Global Lagrange multiplier dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> glmdofrowmap_;

    //! Global slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gsdofrowmap_;

    //! Global master dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gmdofrowmap_;

    //! Global internal dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gndofrowmap_;

    //! Global slave and master dof row map (slave+master map)
    Teuchos::RCP<Epetra_Map> gsmdofrowmap_;

    //! Global displacement dof row map (s+m+n map)
    Teuchos::RCP<Epetra_Map> gdisprowmap_;

    //! Global slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gsnoderowmap_;

    //! Global master node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gmnoderowmap_;

    //! @name Parallel redistribution
    //!@{

    //! Global Lagrange multiplier dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pglmdofrowmap_;

    //! Global slave dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pgsdofrowmap_;

    //! Global master dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pgmdofrowmap_;

    //! Global slave and master dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pgsmdofrowmap_;

    //! Global dirichlet toggle of all slave dofs (before parallel redistribution)
    Teuchos::RCP<Epetra_Vector> pgsdirichtoggle_;

    //!@}

    //! @name Binning strategy
    //!@{

    //! Initial element column map for binning strategy (slave and master)
    std::vector<Teuchos::RCP<Epetra_Map>> initial_elecolmap_;

    //! Global Mortar matrix \f$D\f$
    Teuchos::RCP<CORE::LINALG::SparseMatrix> dmatrix_;

    //! Global Mortar matrix \f$M\f$
    Teuchos::RCP<CORE::LINALG::SparseMatrix> mmatrix_;

    //! Global weighted gap vector \f$g\f$
    Teuchos::RCP<Epetra_Vector> g_;

    //! Global constraint right-hand side vector (only for saddlepoint problems)
    Teuchos::RCP<Epetra_Vector> constrrhs_;

    //! Current vector of Lagrange multipliers at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> z_;

    //! Old vector of Lagrange multipliers at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> zold_;

    /*! \brief Lagrange multiplier vector increment within SaddlePointSolve
     *
     *  \remark This is \em not the increment of #z_ between \f$t_{n+1}\f$ and \f$t_{n}\f$!)
     */
    Teuchos::RCP<Epetra_Vector> zincr_;

    //! Vector of Lagrange multipliers from last Uzawa step
    Teuchos::RCP<Epetra_Vector> zuzawa_;

    //! @name Status flags
    //!@{

    /*! \brief Flag indicating whether transformation should be applied
     *
     * \todo What transformation?
     */
    bool dualquadslavetrafo_;

    //!@}


    /*! \brief Transformation matrix \f$T\f$ for dual quad 3D case
     *
     * \todo What is matrix \f$T\f$?
     * \todo What is quad? Quadratic shape functions or quadrilateral elements?
     * \todo What is the difference to #systrafo_?
     */
    Teuchos::RCP<CORE::LINALG::SparseMatrix> trafo_;

    /*! \brief Transformation matrix \f$T\f$ for dual quad 3D case
     *
     * \todo What is matrix \f$T\f$?
     * \todo What is quad? Quadratic shape functions or quadrilateral elements?
     * \todo What is the difference to #trafo_?
     */
    Teuchos::RCP<CORE::LINALG::SparseMatrix> systrafo_;

    /*! \brief Inverse trafo matrix \f$T^{-1}\f$ for dual quad 3D case
     *
     * \todo What is matrix \f$T\f$?
     * \todo What is quad? Quadratic shape functions or quadrilateral elements?
     */
    Teuchos::RCP<CORE::LINALG::SparseMatrix> invtrafo_;

    /*! \brief Integration time
     *
     * \todo Is this the wall clock time required to perform the mortar integration?
     */
    double inttime_;

    //! Structural force
    Teuchos::RCP<Epetra_Vector> f_;

    //! Structural force (slave)
    Teuchos::RCP<Epetra_Vector> fs_;

    //! Matrix containing \f$D\f$ and \f$-M\f$
    Teuchos::RCP<CORE::LINALG::SparseMatrix> dm_matrix_;

    /*! \brief Matrix containing D and -M. transposed
     *
     * \todo Is it \f$D\f$ and \f$-M^T\f$ or (D and -M)^T? In latter case, why store it explicitly?
     */
    Teuchos::RCP<CORE::LINALG::SparseMatrix> dm_matrix_t_;

    //! Lagrange multiplier diagonal block
    Teuchos::RCP<CORE::LINALG::SparseMatrix> lm_diag_matrix_;

   private:
    ///< pointer to the NOX::NLN::CONSTRAINT::Interface::Required object
    Teuchos::RCP<CONTACT::MtNoxInterface> noxinterface_ptr_;

  };  // class MtAbstractStrategy
}  // namespace CONTACT

// << operator
std::ostream& operator<<(std::ostream& os, const CONTACT::MtAbstractStrategy& strategy);

FOUR_C_NAMESPACE_CLOSE

#endif
