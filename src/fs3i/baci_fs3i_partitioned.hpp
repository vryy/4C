/*----------------------------------------------------------------------*/
/*! \file
\brief H-file associated with general algorithmic routines for
       partitioned solution approaches to fluid-structure-scalar-scalar
       interaction (FS3I), that is, algorithmic routines not specifically
       related to partitioned solution approaches to one -or
       two-way-coupled problem configurations, respectively

\level 2



*----------------------------------------------------------------------*/


#ifndef BACI_FS3I_PARTITIONED_HPP
#define BACI_FS3I_PARTITIONED_HPP


#include "baci_config.hpp"

#include "baci_coupling_adapter_volmortar.hpp"
#include "baci_fs3i.hpp"
#include "baci_inpar_fs3i.hpp"

#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}
namespace ADAPTER
{
  class Coupling;
  class MortarVolCoupl;
  class ScaTraBaseAlgorithm;
}  // namespace ADAPTER

namespace FSI
{
  class Monolithic;

  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }  // namespace UTILS
}  // namespace FSI

namespace CORE::LINALG
{
  class MultiMapExtractor;
  class BlockSparseMatrixBase;
  class SparseMatrix;
  class Solver;
}  // namespace CORE::LINALG


namespace FS3I
{
  class PartFS3I : public FS3I_Base
  {
   public:
    //! constructor of base class for partitioned FS3I
    PartFS3I(const Epetra_Comm& comm);

    //! initialize this class
    void Init() override;

    //! setup this class
    void Setup() override;

    //! @name overall FS3I system
    //@{

    //! time loop to be defined in inherited classes (structure depends on
    //! considered coupling, i.e. one-way or two-way)
    void Timeloop() override = 0;

    //! flag whether time loop should be finished
    bool NotFinished() const { return ((step_ < numstep_) and ((time_ + 1e-14) < timemax_)); };

    //! read and set fields needed for restart
    void ReadRestart() override;

    /// redistribute the  FPSI interface, if running on parallel. Just needed in the case of FPS3I
    void RedistributeInterface() override { return; };

    /// create a volmortar object
    Teuchos::RCP<CORE::ADAPTER::MortarVolCoupl> CreateVolMortarObject(
        Teuchos::RCP<DRT::Discretization> masterdis, Teuchos::RCP<DRT::Discretization> slavedis);

    //! set-up of FSI and ScaTra systems
    void SetupSystem() override;

    //! test results for individual fields
    void TestResults(const Epetra_Comm& comm) override;

    //! information transfer FSI -> ScaTra
    void SetFSISolution();

    /// set scatra solution on structure field
    void SetStructScatraSolution() const;

    //! check convergence of monolithic ScaTra problem (depends on which
    // coupling is considered)
    virtual bool ScatraConvergenceCheck(int itnum) = 0;

    //! return communicator
    const Epetra_Comm& Comm() const { return comm_; }

    /// extract fluid convective and structure convective velocities
    void ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector>>& vel,
        std::vector<Teuchos::RCP<const Epetra_Vector>>& convel) const;

    void SetVelocityFields() const;

    //! routine for preparing time step to be defined in inherited classes
    //! (structure depends on coupling, that is, either one- or two-way)
    virtual void PrepareTimeStep() = 0;

    void SetMeshDisp() const;

    /// provide wall shear stresses from FS3I subproblem for scatra subproblem
    virtual void SetWallShearStresses() const;

    /// extract Wall Shear Stresses at the interface
    void ExtractWSS(std::vector<Teuchos::RCP<const Epetra_Vector>>& wss) const;

    Teuchos::ParameterList& ManipulateFsiTimeParams(const Teuchos::ParameterList& fs3idyn) const;

    //@}

    /// transport quantity from fluid to fluid-scalar
    Teuchos::RCP<const Epetra_Vector> FluidToFluidScalar(
        const Teuchos::RCP<const Epetra_Vector> fluidvector) const;

    /// transport quantity from fluid-scalar to fluid
    Teuchos::RCP<const Epetra_Vector> FluidScalarToFluid(
        const Teuchos::RCP<const Epetra_Vector> fluidscalarvector) const;

    /// transport quantity from structure to structure-scalar
    Teuchos::RCP<const Epetra_Vector> StructureToStructureScalar(
        const Teuchos::RCP<const Epetra_Vector> structurevector) const;

    /// transport quantity from structure-scalar to structure
    Teuchos::RCP<const Epetra_Vector> StructureScalarToStructure(
        const Teuchos::RCP<const Epetra_Vector> structurescalavector) const;

   private:
    /// transport quantity from i-th volmortar master to i-th volmortar slave
    Teuchos::RCP<const Epetra_Vector> VolMortarMasterToSlavei(
        const int i, const Teuchos::RCP<const Epetra_Vector> mastervector) const;

    /// transport quantity from i-th volmortar slave to i-th volmortar master
    Teuchos::RCP<const Epetra_Vector> VolMortarSlaveToMasteri(
        const int i, const Teuchos::RCP<const Epetra_Vector> slavevector) const;

   protected:
    /// fsi algorithm
    Teuchos::RCP<FSI::Monolithic> fsi_;

    /// vector of scatra volume couplings (i.e. fluid to fluid-scalar and structure to
    /// structure-scalar)
    std::vector<INPAR::FS3I::VolumeCoupling> volume_fieldcouplings_;

   private:
    //! volume coupling (using mortar) adapter
    std::vector<Teuchos::RCP<CORE::ADAPTER::MortarVolCoupl>> volume_coupling_objects_;

    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra_;

    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra_;
  };
}  // namespace FS3I

BACI_NAMESPACE_CLOSE

#endif
