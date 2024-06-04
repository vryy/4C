/*----------------------------------------------------------------------*/
/*! \file
\brief Generic class for all mortar solution strategies

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MORTAR_STRATEGY_BASE_HPP
#define FOUR_C_MORTAR_STRATEGY_BASE_HPP

#include "4C_config.hpp"

#include "4C_inpar_contact.hpp"     // for the INPAR::CONTACT enums
#include "4C_mortar_interface.hpp"  // for the enum state type
#include "4C_solver_nonlin_nox_constraint_interface_preconditioner.hpp"  // interface specifications

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace INPAR
{
  namespace STR
  {
    enum DynamicType : int;
  }  // namespace STR
}  // namespace INPAR
namespace DRT
{
  class Discretization;
}

namespace CORE::LINALG
{
  class MapExtractor;
  class Solver;
  class SparseOperator;
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace CORE::IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace CORE::IO

namespace MORTAR
{
  /*! \brief Data container object for the strategy base
   *
   *  This object makes it possible to interchange and share the current state of the
   *  contact simulation between different strategy objects. By using this the
   *  actual strategy stays stateless!
   *
   *  \author  hiermeier
   *  \date 05/16 */
  class StratDataContainer
  {
   public:
    //! constructor
    StratDataContainer();

    //! destructor
    virtual ~StratDataContainer() = default;

    //! Return underlying problem dof row map (not only interfaces)
    Teuchos::RCP<Epetra_Map>& ProbDofsPtr() { return probdofs_; };
    Teuchos::RCP<const Epetra_Map> ProbDofsPtr() const { return probdofs_; };

    //! Return underlying problem node row map (not only interfaces)
    Teuchos::RCP<Epetra_Map>& ProbNodesPtr() { return probnodes_; };
    Teuchos::RCP<const Epetra_Map> ProbNodesPtr() const { return probnodes_; };

    //! Return communicator
    Teuchos::RCP<const Epetra_Comm>& CommPtr() { return comm_; };
    Teuchos::RCP<const Epetra_Comm> CommPtr() const { return comm_; };

    //! Return containing contact input parameters
    Teuchos::ParameterList& SContact() { return scontact_; };
    const Teuchos::ParameterList& SContact() const { return scontact_; };

    //! Return dimension of problem (2D or 3D)
    int& Dim() { return dim_; };
    const int& Dim() const { return dim_; };

    //! Return generalized-alpha parameter (0.0 for statics)
    double& AlphaF() { return alphaf_; };
    const double& AlphaF() const { return alphaf_; };

    /// get the (dynamic) time integration type
    inline INPAR::STR::DynamicType GetDynType() const { return dyntype_; };

    /// return dynamic time integration parameter
    inline double GetDynParameterN() const { return dynparam_n_; }

    /// set dynamic time integration parameter
    inline void SetDynParameterN(const double dynparamN) { dynparam_n_ = dynparamN; }

    /// set the (dynamic) time integration type
    inline void SetDynType(INPAR::STR::DynamicType dyntype) { dyntype_ = dyntype; }

    //! Return flag indicating parallel redistribution status
    bool& IsParRedist() { return parredist_; };
    const bool& IsParRedist() const { return parredist_; };

    //! Return highest dof number in problem discretization
    int& MaxDof() { return maxdof_; };
    const int& MaxDof() const { return maxdof_; };

    //! Return current used system type
    INPAR::CONTACT::SystemType& SysType() { return systype_; };
    const INPAR::CONTACT::SystemType& SysType() const { return systype_; };

   private:
    //! Underlying problem dof row map (not only interfaces)
    Teuchos::RCP<Epetra_Map> probdofs_;

    //! Underlying problem node row map (not only interfaces)
    Teuchos::RCP<Epetra_Map> probnodes_;

    //! Communicator
    Teuchos::RCP<const Epetra_Comm> comm_;

    //! Containing contact input parameters
    Teuchos::ParameterList scontact_;

    //! Dimension of problem (2D or 3D)
    int dim_;

    //! Generalized-alpha parameter (0.0 for statics)
    double alphaf_;

    //! Flag indicating parallel redistribution status
    bool parredist_;

    //! Highest dof number in problem discretization
    int maxdof_;

    //! Current used system type
    INPAR::CONTACT::SystemType systype_;

    //! time integration type
    INPAR::STR::DynamicType dyntype_;

    //! time integration parameter for the contributions of the old/previous time step
    double dynparam_n_;
  };

  /*!
  \brief Abstract base class for mortar solution strategies

  Every specific solution algorithm (e.g. mortar contact with Lagrange multipliers or
  mortar meshtying with penalty method) has to be specified in a corresponding derived
  subclass defining the concrete algorithmic steps.

  */
  class StrategyBase : public NOX::NLN::CONSTRAINT::Interface::Preconditioner
  {
   public:
    //! @name Enums and Friends
    //! @{
    // can be called by store_nodal_quantities() or StoreDMtoNodes()
    enum QuantityType
    {
      lmcurrent,  //!< current lagr. mult.
      lmold,      //!< lagr. mult. for last converged state
      lmupdate,   //!< update current lagr. mult. (same as for lmcurrent + DBC check)
      lmuzawa,    //!< lagr. mutl. from last Uzawa step
      activeold,  //!< contact status of last converged state
      slipold,    //!< slip for last converged state
      dm,
      pentrac,
      weightedwear,  //!< weighted wear (internal state var. approach)
      wupdate,       //!< update current pv wear for current step (slave)
      wmupdate,      //!< update current pv wear for current step (master)
      wold,          //!< pv wear for last converged state (slave)
      wmold,         //!< pv wear for last converged state (master)
      wupdateT,      //!< accumulated pv wear for different time scales
      lmThermo,      //!< thermal Lagrange multiplier
      n_old          //!< old normal
    };
    //! @}

    /*!
    \brief Standard Constructor

    Creates the strategy base object and initializes all global variables.

    \param data_ptr (in): data container object
    \param dofrowmap (in): dofrowmap of underlying problem
    \param noderowmap (in): noderowmap of underlying problem
    \param elementrowmap (in): elementrowmap of underlying problem
    \param params (in): List of meshtying/contact parameters
    \param spatialDim (in): Global problem dimension
    \param comm (in): A communicator object
    \param alphaf (in): Midpoint for Gen-alpha time integration
    \param maxdof (in): Highest dof number in global problem

    */
    StrategyBase(const Teuchos::RCP<MORTAR::StratDataContainer>& data_ptr,
        const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params, const int spatialDim,
        const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf, const int maxdof);


    //! @name Access methods
    //! @{
    //! Get parameter list
    Teuchos::ParameterList& Params() { return scontact_; }
    const Teuchos::ParameterList& Params() const { return scontact_; }

    //! return the current system type
    const INPAR::CONTACT::SystemType& SystemType() const { return systype_; };

    //! Get problem dimension
    int Dim() const { return dim_; }

    //! Get Epetra communicator
    const Epetra_Comm& Comm() const { return *comm_; }

    //! Get the underlying problem dof row map
    const Teuchos::RCP<Epetra_Map>& ProblemDofs() { return probdofs_; };
    Teuchos::RCP<const Epetra_Map> ProblemDofs() const { return probdofs_; };

    //! Get the underlying problem node row map
    const Teuchos::RCP<Epetra_Map>& ProblemNodes() { return probnodes_; };
    Teuchos::RCP<const Epetra_Map> ProblemNodes() const { return probnodes_; };

    //@}

    /// Set the time integration information
    void set_time_integration_info(const double time_fac, const INPAR::STR::DynamicType dyntype);

    //! @name Purely virtual functions

    // All these functions are defined in one or more specific derived classes,
    // such as CONTACT::ContactLagrangeStrategy or CONTACT::MeshtyingPenaltyStrategy.
    // As the base class MORTAR::StrategyBase is always called from the control routine
    // (time integrator), these functions need to be defined purely virtual here.

    virtual Teuchos::RCP<Epetra_Map> SlaveRowNodes() = 0;
    virtual Teuchos::RCP<Epetra_Map> ActiveRowNodes() = 0;
    virtual Teuchos::RCP<Epetra_Map> ActiveRowDofs() = 0;
    virtual Teuchos::RCP<Epetra_Map> not_re_dist_slave_row_dofs() = 0;
    virtual Teuchos::RCP<Epetra_Map> not_re_dist_master_row_dofs() = 0;
    virtual bool ActiveSetConverged() = 0;
    virtual bool active_set_semi_smooth_converged() const = 0;
    virtual void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<CORE::LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f,
        const int step, const int iter, bool predictor = false) = 0;
    virtual void AssembleMortar() = 0;
    virtual void collect_maps_for_preconditioner(Teuchos::RCP<Epetra_Map>& MasterDofMap,
        Teuchos::RCP<Epetra_Map>& SlaveDofMap, Teuchos::RCP<Epetra_Map>& InnerDofMap,
        Teuchos::RCP<Epetra_Map>& ActiveDofMap) const = 0;
    virtual double ConstraintNorm() const = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactNorStress() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactTanStress() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactNorForce() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactTanForce() = 0;
    virtual Teuchos::RCP<CORE::LINALG::SparseMatrix> DMatrix() = 0;
    virtual void DoReadRestart(
        CORE::IO::DiscretizationReader& reader, Teuchos::RCP<const Epetra_Vector> dis) = 0;
    virtual void DoWriteRestart(std::map<std::string, Teuchos::RCP<Epetra_Vector>>& restart_vectors,
        bool forcedrestart = false) const = 0;
    virtual void Evaluate(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void EvaluateMeshtying(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual Teuchos::RCP<CORE::LINALG::SparseMatrix> EvaluateNormals(
        Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void evaluate_reference_state() = 0;
    virtual void EvaluateRelMov() = 0;
    virtual void evaluate_rel_mov_predict() = 0;
    virtual bool Friction() const = 0;
    virtual void InitEvalInterface() = 0;
    virtual void InitMortar() = 0;
    virtual void Initialize() = 0;
    virtual void InitializeUzawa(
        Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
    virtual double InitialPenalty() = 0;
    virtual void InterfaceForces(bool output = false) = 0;
    virtual double Inttime() = 0;
    virtual void Inttime_init() = 0;
    virtual bool IsInContact() const = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMult() = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMultOld() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ConstrRhs() = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMultSolveIncr() = 0;
    virtual Teuchos::RCP<const Epetra_Vector> MeshInitialization() = 0;
    virtual Teuchos::RCP<CORE::LINALG::SparseMatrix> MMatrix() = 0;
    virtual void MortarCoupling(const Teuchos::RCP<const Epetra_Vector>& dis) = 0;
    virtual int NumberOfActiveNodes() const = 0;
    virtual int NumberOfSlipNodes() const = 0;
    virtual void compute_contact_stresses() = 0;

    /*!
    \brief Write results for visualization separately for each meshtying/contact interface

    Call each interface, such that each interface can handle its own output of results.

    \param[in] outputParams Parameter list with stuff required by interfaces to write output
    */
    virtual void postprocess_quantities_per_interface(
        Teuchos::RCP<Teuchos::ParameterList> outputParams) = 0;

    virtual void Print(std::ostream& os) const = 0;
    virtual void PrintActiveSet() const = 0;
    virtual void Recover(Teuchos::RCP<Epetra_Vector> disi) = 0;
    virtual bool RedistributeContact(
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel) = 0;
    virtual void redistribute_meshtying() = 0;
    virtual void ResetActiveSet() = 0;
    virtual void ResetPenalty() = 0;
    virtual void ModifyPenalty() = 0;
    virtual void restrict_meshtying_zone() = 0;
    virtual void build_saddle_point_system(Teuchos::RCP<CORE::LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) = 0;
    virtual void update_displacements_and_l_mincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) = 0;
    virtual void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) = 0;
    virtual void set_state(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) = 0;
    virtual Teuchos::RCP<Epetra_Map> SlipRowNodes() = 0;
    virtual void store_dirichlet_status(Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmaps) = 0;
    virtual void store_nodal_quantities(MORTAR::StrategyBase::QuantityType type) = 0;
    virtual void Update(Teuchos::RCP<const Epetra_Vector> dis) = 0;
    virtual void UpdateActiveSet() = 0;
    virtual void update_active_set_semi_smooth(const bool firstStepPredictor = false) = 0;
    virtual void update_uzawa_augmented_lagrange() = 0;
    virtual void update_constraint_norm(int uzawaiter = 0) = 0;
    virtual void VisualizeGmsh(const int step, const int iter) = 0;
    virtual bool WasInContact() const = 0;
    virtual bool was_in_contact_last_time_step() const = 0;

    // Flag for Poro No Penetration Condition (overloaded by LagrangeStrategyPoro)
    virtual bool has_poro_no_penetration() const { return false; }

    // Nitsche stuff
    virtual bool IsNitsche() const { return false; }

    // wear stuff
    virtual bool WeightedWear() const { return false; };
    virtual bool WearBothDiscrete() const { return false; };
    virtual Teuchos::RCP<Epetra_Vector> WearRhs() { return Teuchos::null; };
    virtual Teuchos::RCP<Epetra_Vector> WearMRhs() { return Teuchos::null; };
    virtual Teuchos::RCP<Epetra_Vector> WSolveIncr() { return Teuchos::null; };
    virtual Teuchos::RCP<Epetra_Vector> WMSolveIncr() { return Teuchos::null; };
    virtual Teuchos::RCP<Epetra_Vector> ContactWear() { return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> ContactWear() const { return Teuchos::null; };
    virtual void OutputWear() { ; };
    virtual Teuchos::RCP<const Epetra_Map> MasterSlipNodes() const { return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Map> MasterActiveNodes() const { return Teuchos::null; };

    // constraint preconditioner functions
    bool IsSaddlePointSystem() const override = 0;
    bool IsCondensedSystem() const override = 0;
    virtual bool IsPenalty() const = 0;
    void fill_maps_for_preconditioner(
        std::vector<Teuchos::RCP<Epetra_Map>>& maps) const override = 0;
    //@}

   private:
    /*! return the mutable mortar data container
     *
     * \remark This has to stay PRIVATE, otherwise this function becomes ambiguous.
     *
     * \author hiermeier
     * \date 05/16 */
    MORTAR::StratDataContainer& data() { return *data_ptr_; };

    /*! return the read-only mortar data container
     *
     * \remark This has to stay PRIVATE, otherwise this function becomes ambiguous.
     *
     * \author hiermeier
     * \date 05/16 */
    const MORTAR::StratDataContainer& data() const { return *data_ptr_; };

   protected:
    // don't want cctor (= operator impossible anyway for abstract class)
    StrategyBase(const StrategyBase& old);

    /*! @name References to the data container content
     *
     * \remark Please add no new member variables to the strategy base! Use
     *  the corresponding data container instead (--> MORTAR::StratDataContainer).
     *  If you have any questions concerning this, do not hesitate and ask me.
     *                                                          hiermeier 05/16 */
    //! @{
    Teuchos::RCP<Epetra_Map>&
        probdofs_;  //!< ref. to underlying problem dof row map (not only interfaces)
    Teuchos::RCP<Epetra_Map>&
        probnodes_;  //!< ref. to underlying problem node row map (not only interfaces)

    Teuchos::RCP<const Epetra_Comm>& comm_;  //!< ref. to communicator
    Teuchos::ParameterList& scontact_;       //!< ref. to containing contact input parameters
    int& dim_;                               //!< ref. to dimension of problem (2D or 3D)
    double& alphaf_;   //!< ref. to generalized-alpha parameter (0.0 for statics)
    bool& parredist_;  //!< ref. to flag indicating parallel redistribution status
    int& maxdof_;      //!< ref. to highest dof number in problem discretization
    INPAR::CONTACT::SystemType& systype_;  //!< ref. to current used system type
                                           //! @}

   private:
    //! pointer to the data container object
    Teuchos::RCP<MORTAR::StratDataContainer> data_ptr_;

  };  // class StrategyBase
}  // namespace MORTAR

FOUR_C_NAMESPACE_CLOSE

#endif
