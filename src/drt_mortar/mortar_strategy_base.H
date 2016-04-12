/*!----------------------------------------------------------------------
\file mortar_strategy_base.H

\maintainer Philipp Farah, Alexander Seitz

*----------------------------------------------------------------------*/
#ifndef MORTAR_STRATEGY_BASE_H
#define MORTAR_STRATEGY_BASE_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"

// forward declarations
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
}

namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}

namespace MORTAR
{

/*!
\brief Abstract base class for mortar solution strategies

Every specific solution algorithm (e.g. mortar contact with Lagrange multipliers or
mortar meshtying with penalty method) has to be specified in a corresponding derived
subclass defining the concrete algorithmic steps.

\author popp (popp@lnm.mw.tum.de)
*/
class StrategyBase
{
  public:

    //! @name Enums and Friends
    enum QuantityType // can be called by StoreNodalQuantities() or StoreDMtoNodes()
    {
      lmcurrent,   // current lagr. mult.
      lmold,       // lagr. mult. for last converged state
      lmupdate,    // update current lagr. mult. (same as for lmcurrent + DBC check)
      lmuzawa,     // lagr. mutl. from last Uzawa step
      activeold,   // contact status of last converged state
      slipold,     // slip for last converged state
      dm,
      pentrac,
      weightedwear,// weighted wear (internal state var. approach)
      wupdate,     // update current pv wear for current step (slave)
      wmupdate,    // update current pv wear for current step (master)
      wold,        // pv wear for last converged state (slave)
      wmold,       // pv wear for last converged state (master)
      wupdateT,    // accumulated pv wear for different time scales
      lmThermo     // thermal Lagrange multiplier
    };
    //@}

    /*!
    \brief Standard Constructor

    Creates the strategy base object and initializes all global variables.

    \param dofrowmap (in): dofrowmap of underlying problem
    \param noderowmap (in): noderowmap of underlying problem
    \param elementrowmap (in): elementrowmap of underlying problem
    \param params (in): List of meshtying/contact parameters
    \param dim (in): Global problem dimension
    \param comm (in): A communicator object
    \param alphaf (in): Midpoint for Gen-alpha time integration
    \param maxdof (in): Highest dof number in global problem

    */
    StrategyBase(
        const Epetra_Map* DofRowMap,
        const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params,
        int dim,
        Teuchos::RCP<Epetra_Comm> comm,
        double alphaf,
        int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~StrategyBase() {};

    //! @name Access methods

    /*!
    \brief Get parameter list

    */
    Teuchos::ParameterList& Params() { return scontact_; }

    /*!
    \brief Get problem dimension

    */
    const int& Dim() { return dim_; }

    /*!
    \brief Get Epetra communicator

    */
    const Epetra_Comm& Comm() const { return *comm_; }

    /*!
    \brief Get the underlying problem dof row map

    */
    const Teuchos::RCP<Epetra_Map> ProblemDofs() { return probdofs_; }

    /*!
    \brief Get the underlying problem node row map

    */
    const Teuchos::RCP<Epetra_Map> ProblemNodes() { return probnodes_; }

    //@}

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
    virtual bool ActiveSetSemiSmoothConverged() = 0;
    virtual void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step, const int iter, bool predictor = false) = 0;
    virtual void AssembleMortar() = 0;
    virtual void InitBinStrategyforTimestep(Teuchos::RCP<Epetra_Vector> vel) = 0;
    virtual void CollectMapsForPreconditioner(Teuchos::RCP<Epetra_Map>& MasterDofMap, Teuchos::RCP<Epetra_Map>& SlaveDofMap, Teuchos::RCP<Epetra_Map>& InnerDofMap, Teuchos::RCP<Epetra_Map>& ActiveDofMap ) = 0;
    virtual double ConstraintNorm() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactNorStress() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactTanStress() = 0;
    virtual Teuchos::RCP<LINALG::SparseMatrix> DMatrix() = 0;
    virtual void DoReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void DoWriteRestart(std::map<std::string,Teuchos::RCP<Epetra_Vector> >& restart_vectors, bool forcedrestart = false) = 0;
    virtual void Evaluate(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void EvaluateMeshtying(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual Teuchos::RCP<LINALG::SparseMatrix> EvaluateNormals(Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void EvaluateReferenceState(const Teuchos::RCP<Epetra_Vector> vec) = 0;
    virtual void EvaluateRelMov() = 0;
    virtual void EvaluateRelMovPredict() = 0;
    virtual bool Friction() = 0;
    virtual void InitEvalInterface() = 0;
    virtual void InitMortar() = 0;
    virtual void Initialize() = 0;
    virtual void InitializeUzawa(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
    virtual double InitialPenalty() = 0;
    virtual void InterfaceForces(bool output = false) = 0;
    virtual double Inttime() = 0;
    virtual void Inttime_init()=0;
    virtual bool IsInContact() = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMult() = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMultOld() = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMultOldRescaled() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ConstrRhs() = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMultSolveIncr() = 0;
    virtual Teuchos::RCP<Epetra_Vector> MeshInitialization() = 0;
    virtual Teuchos::RCP<LINALG::SparseMatrix> MMatrix() = 0;
    virtual void MortarCoupling(const Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual int NumberOfActiveNodes() = 0;
    virtual int NumberOfSlipNodes() = 0;
    virtual void OutputStresses() = 0;
    virtual void AugForces(Epetra_Vector& augfs_lm,Epetra_Vector& augfs_g,Epetra_Vector& augfm_lm,Epetra_Vector& augfm_g) = 0;
    virtual void Print(std::ostream& os) const = 0;
    virtual void PrintActiveSet() = 0;
    virtual void Recover(Teuchos::RCP<Epetra_Vector> disi) = 0;
    virtual void RedistributeContact(Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void RedistributeMeshtying() = 0;
    virtual void ResetActiveSet() = 0;
    virtual void ResetPenalty() = 0;
    virtual void ModifyPenalty() = 0;
    virtual void RestrictMeshtyingZone() = 0;
    virtual void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd, Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<LINALG::MapExtractor> dbcmaps, int numiter,Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) = 0;
    virtual void UpdateDisplacementsAndLMincrements(Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Epetra_Vector> blocksol) = 0;
    virtual void SaveReferenceState(Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void SetState(const std::string& statename, const Teuchos::RCP<Epetra_Vector> vec) = 0;
    virtual const Teuchos::RCP<Epetra_Map> SlipRowNodes() = 0;
    virtual void StoreDirichletStatus(Teuchos::RCP<LINALG::MapExtractor> dbcmaps) = 0;
    virtual void StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type) = 0;
    virtual void Update(Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void UpdateActiveSet() = 0;
    virtual void UpdateActiveSetSemiSmooth() = 0;
    virtual void UpdateUzawaAugmentedLagrange() = 0;
    virtual void UpdateConstraintNorm(int uzawaiter = 0) = 0;
    virtual void VisualizeGmsh(const int step, const int iter) = 0;
    virtual bool WasInContact() = 0;
    virtual bool WasInContactLastTimeStep() = 0;

    //Flag for Poro No Penetration Condition (overloaded by PoroLagrangeStrategy)
    virtual bool HasPoroNoPenetration(){return false;}

    // wear stuff
    virtual bool WeightedWear(){return false;};
    virtual bool WearBothDiscrete(){return false;};
    virtual Teuchos::RCP<Epetra_Vector> WearRhs(){return Teuchos::null;};
    virtual Teuchos::RCP<Epetra_Vector> WearMRhs(){return Teuchos::null;};
    virtual Teuchos::RCP<Epetra_Vector> WSolveIncr(){return Teuchos::null;};
    virtual Teuchos::RCP<Epetra_Vector> WMSolveIncr(){return Teuchos::null;};
    virtual Teuchos::RCP<Epetra_Vector> ContactWear(){return Teuchos::null;};
    virtual void OutputWear(){;};
    virtual const Teuchos::RCP<Epetra_Map> MasterSlipNodes(){return Teuchos::null;};
    virtual const Teuchos::RCP<Epetra_Map> MasterActiveNodes(){return Teuchos::null;};
    //@}

  protected:

    // don't want cctor (= operator impossible anyway for abstract class)
    StrategyBase(const StrategyBase& old);

    Teuchos::RCP<Epetra_Map>    probdofs_;         // underlying problem dof row map (not only interfaces)
    Teuchos::RCP<Epetra_Map>    probnodes_;        // underlying problem node row map (not only interfaces)

    Teuchos::RCP<Epetra_Comm>   comm_;             // communicator
    Teuchos::ParameterList      scontact_;         // containing contact input parameters
    int                         dim_;              // dimension of problem (2D or 3D)
    double                      alphaf_;           // Generalized-alpha parameter (0.0 for statics)
    bool                        parredist_;        // flag indicating parallel redistribution status
    int                         maxdof_;           // highest dof number in problem discretization

}; // class StrategyBase
} // namespace MORTAR

#endif  // #ifndef MORTAR_STRATEGY_BASE_H
