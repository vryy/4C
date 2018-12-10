/*----------------------------------------------------------------------*/
/*!
\file mortar_strategy_base.H

\brief Generic class for all mortar solution strategies

\level 2

\maintainer Matthias Mayr

*/
/*----------------------------------------------------------------------*/
#ifndef MORTAR_STRATEGY_BASE_H
#define MORTAR_STRATEGY_BASE_H

#include "../solver_nonlin_nox/nox_nln_constraint_interface_preconditioner.H"  // interface specifications
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"

#include "mortar_interface.H"            // for the enum state type
#include "../drt_inpar/inpar_contact.H"  // for the INPAR::CONTACT enums


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

namespace LINALG
{
  class MapExtractor;
  class Solver;
  class SparseOperator;
  class SparseMatrix;
}  // namespace LINALG

namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace IO

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
    virtual ~StratDataContainer(){};

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
    inline double GetDynParameterN() const { return dynparamN_; }

    /// set dynamic time integration parameter
    inline void SetDynParameterN(const double dynparamN) { dynparamN_ = dynparamN; }

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
    double dynparamN_;
  };

  /*!
  \brief Abstract base class for mortar solution strategies

  Every specific solution algorithm (e.g. mortar contact with Lagrange multipliers or
  mortar meshtying with penalty method) has to be specified in a corresponding derived
  subclass defining the concrete algorithmic steps.

  \author popp (popp@lnm.mw.tum.de)
  */
  class StrategyBase : public NOX::NLN::CONSTRAINT::Interface::Preconditioner
  {
   public:
    //! @name Enums and Friends
    //! @{
    // can be called by StoreNodalQuantities() or StoreDMtoNodes()
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
    \param dim (in): Global problem dimension
    \param comm (in): A communicator object
    \param alphaf (in): Midpoint for Gen-alpha time integration
    \param maxdof (in): Highest dof number in global problem

    */
    StrategyBase(const Teuchos::RCP<MORTAR::StratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params, int dim, const Teuchos::RCP<const Epetra_Comm>& comm,
        double alphaf, int maxdof);

    //! Destructor
    virtual ~StrategyBase(){};

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
    void SetTimeIntegrationInfo(const double time_fac, const INPAR::STR::DynamicType dyntype);

    //! @name Purely virtual functions

    // All these functions are defined in one or more specific derived classes,
    // such as CONTACT::ContactLagrangeStrategy or CONTACT::MeshtyingPenaltyStrategy.
    // As the base class MORTAR::StrategyBase is always called from the control routine
    // (time integrator), these functions need to be defined purely virtual here.

    virtual const Teuchos::RCP<Epetra_Map> SlaveRowNodes() = 0;
    virtual const Teuchos::RCP<Epetra_Map> ActiveRowNodes() = 0;
    virtual const Teuchos::RCP<Epetra_Map> ActiveRowDofs() = 0;
    virtual const Teuchos::RCP<Epetra_Map> NotReDistSlaveRowDofs() = 0;
    virtual const Teuchos::RCP<Epetra_Map> NotReDistMasterRowDofs() = 0;
    virtual bool ActiveSetConverged() = 0;
    virtual bool ActiveSetSemiSmoothConverged() const = 0;
    virtual void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
        const int iter, bool predictor = false) = 0;
    virtual void AssembleMortar() = 0;
    virtual void InitBinStrategyforTimestep(Teuchos::RCP<const Epetra_Vector> vel) = 0;
    virtual void CollectMapsForPreconditioner(Teuchos::RCP<Epetra_Map>& MasterDofMap,
        Teuchos::RCP<Epetra_Map>& SlaveDofMap, Teuchos::RCP<Epetra_Map>& InnerDofMap,
        Teuchos::RCP<Epetra_Map>& ActiveDofMap) = 0;
    virtual double ConstraintNorm() const = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactNorStress() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactTanStress() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactNorForce() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactTanForce() = 0;
    virtual Teuchos::RCP<LINALG::SparseMatrix> DMatrix() = 0;
    virtual void DoReadRestart(
        IO::DiscretizationReader& reader, Teuchos::RCP<const Epetra_Vector> dis) = 0;
    virtual void DoWriteRestart(std::map<std::string, Teuchos::RCP<Epetra_Vector>>& restart_vectors,
        bool forcedrestart = false) const = 0;
    virtual void Evaluate(Teuchos::RCP<LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void EvaluateMeshtying(Teuchos::RCP<LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual Teuchos::RCP<LINALG::SparseMatrix> EvaluateNormals(Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void EvaluateReferenceState(Teuchos::RCP<const Epetra_Vector> vec) = 0;
    virtual void EvaluateRelMov() = 0;
    virtual void EvaluateRelMovPredict() = 0;
    virtual bool Friction() const = 0;
    virtual void InitEvalInterface() = 0;
    virtual void InitMortar() = 0;
    virtual void Initialize() = 0;
    virtual void InitializeUzawa(
        Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
    virtual double InitialPenalty() = 0;
    virtual void InterfaceForces(bool output = false) = 0;
    virtual double Inttime() = 0;
    virtual void Inttime_init() = 0;
    virtual bool IsInContact() const = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMult() = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMultOld() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ConstrRhs() = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMultSolveIncr() = 0;
    virtual Teuchos::RCP<Epetra_Vector> MeshInitialization() = 0;
    virtual Teuchos::RCP<LINALG::SparseMatrix> MMatrix() = 0;
    virtual void MortarCoupling(const Teuchos::RCP<const Epetra_Vector>& dis) = 0;
    virtual int NumberOfActiveNodes() const = 0;
    virtual int NumberOfSlipNodes() const = 0;
    virtual void OutputStresses() = 0;
    virtual void Print(std::ostream& os) const = 0;
    virtual void PrintActiveSet() const = 0;
    virtual void Recover(Teuchos::RCP<Epetra_Vector> disi) = 0;
    virtual bool RedistributeContact(Teuchos::RCP<const Epetra_Vector> dis) = 0;
    virtual void RedistributeMeshtying() = 0;
    virtual void ResetActiveSet() = 0;
    virtual void ResetPenalty() = 0;
    virtual void ModifyPenalty() = 0;
    virtual void RestrictMeshtyingZone() = 0;
    virtual void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) = 0;
    virtual void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) = 0;
    virtual void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) = 0;
    virtual void SetState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) = 0;
    virtual const Teuchos::RCP<Epetra_Map> SlipRowNodes() = 0;
    virtual void StoreDirichletStatus(Teuchos::RCP<const LINALG::MapExtractor> dbcmaps) = 0;
    virtual void StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type) = 0;
    virtual void Update(Teuchos::RCP<const Epetra_Vector> dis) = 0;
    virtual void UpdateActiveSet() = 0;
    virtual void UpdateActiveSetSemiSmooth(const bool firstStepPredictor = false) = 0;
    virtual void UpdateUzawaAugmentedLagrange() = 0;
    virtual void UpdateConstraintNorm(int uzawaiter = 0) = 0;
    virtual void VisualizeGmsh(const int step, const int iter) = 0;
    virtual bool WasInContact() const = 0;
    virtual bool WasInContactLastTimeStep() const = 0;

    // Flag for Poro No Penetration Condition (overloaded by PoroLagrangeStrategy)
    virtual bool HasPoroNoPenetration() const { return false; }

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
    virtual bool IsSaddlePointSystem() const = 0;
    virtual bool IsCondensedSystem() const = 0;
    virtual bool IsPenalty() const = 0;
    virtual void FillMapsForPreconditioner(std::vector<Teuchos::RCP<Epetra_Map>>& maps) const = 0;
    //@}

   private:
    /*! return the mutable mortar data container
     *
     * \remark This has to stay PRIVATE, otherwise this function becomes ambiguous.
     *
     * \author hiermeier
     * \date 05/16 */
    MORTAR::StratDataContainer& Data() { return *data_ptr_; };

    /*! return the read-only mortar data container
     *
     * \remark This has to stay PRIVATE, otherwise this function becomes ambiguous.
     *
     * \author hiermeier
     * \date 05/16 */
    const MORTAR::StratDataContainer& Data() const { return *data_ptr_; };

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

#endif  // #ifndef MORTAR_STRATEGY_BASE_H
