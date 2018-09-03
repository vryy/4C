/*!----------------------------------------------------------------------
\file fs3i_partitioned.H
\brief H-file associated with general algorithmic routines for
       partitioned solution approaches to fluid-structure-scalar-scalar
       interaction (FS3I), that is, algorithmic routines not specifically
       related to partitioned solution approaches to one -or
       two-way-coupled problem configurations, respectively

\level 2

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10364


*----------------------------------------------------------------------*/


#ifndef PARTFS3I_H
#define PARTFS3I_H


#include "fs3i.H"
#include "../drt_inpar/inpar_fs3i.H"
#include <Teuchos_ParameterList.hpp>


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

namespace LINALG
{
  class MultiMapExtractor;
  class BlockSparseMatrixBase;
  class SparseMatrix;
  class Solver;
}  // namespace LINALG


namespace FS3I
{
  class PartFS3I : public FS3I_Base
  {
   public:
    //! constructor of base class for partitioned FS3I
    PartFS3I(const Epetra_Comm& comm);

    //! destructor of base class for partitioned FS3I
    virtual ~PartFS3I(){};

    //! initialize this class
    virtual void Init();

    //! setup this class
    virtual void Setup();

    //! @name overall FS3I system
    //@{

    //! time loop to be defined in inherited classes (structure depends on
    //! considered coupling, i.e. one-way or two-way)
    virtual void Timeloop() = 0;

    //! flag whether time loop should be finished
    bool NotFinished() const { return ((step_ < numstep_) and ((time_ + 1e-14) < timemax_)); };

    //! read and set fields needed for restart
    virtual void ReadRestart();

    /// redistribute the  FPSI interface, if running on parallel. Just needed in the case of FPS3I
    void RedistributeInterface() { return; };

    /// create a volmortar object
    Teuchos::RCP<::ADAPTER::MortarVolCoupl> CreateVolMortarObject(
        Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<DRT::Discretization> structscatradis);

    //! set-up of FSI and ScaTra systems
    virtual void SetupSystem();

    //! test results for individual fields
    virtual void TestResults(const Epetra_Comm& comm);

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
    const Teuchos::RCP<const Epetra_Vector> FluidToFluidScalar(
        const Teuchos::RCP<const Epetra_Vector> fluidvector) const;

    /// transport quantity from fluid-scalar to fluid
    const Teuchos::RCP<const Epetra_Vector> FluidScalarToFluid(
        const Teuchos::RCP<const Epetra_Vector> fluidscalarvector) const;

    /// transport quantity from structure to structure-scalar
    const Teuchos::RCP<const Epetra_Vector> StructureToStructureScalar(
        const Teuchos::RCP<const Epetra_Vector> structurevector) const;

    /// transport quantity from structure-scalar to structure
    const Teuchos::RCP<const Epetra_Vector> StructureScalarToStructure(
        const Teuchos::RCP<const Epetra_Vector> structurescalavector) const;

   private:
    /// transport quantity from i-th volmortar master to i-th volmortar slave
    const Teuchos::RCP<const Epetra_Vector> VolMortarMasterToSlavei(
        const int i, const Teuchos::RCP<const Epetra_Vector> mastervector) const;

    /// transport quantity from i-th volmortar slave to i-th volmortar master
    const Teuchos::RCP<const Epetra_Vector> VolMortarSlaveToMasteri(
        const int i, const Teuchos::RCP<const Epetra_Vector> slavevector) const;

   protected:
    /// fsi algorithm
    Teuchos::RCP<FSI::Monolithic> fsi_;

    /// vector of scatra volume couplings (i.e. fluid to fluid-scalar and structure to
    /// structure-scalar)
    std::vector<INPAR::FS3I::VolumeCoupling> volume_fieldcouplings_;

   private:
    //! volume coupling (using mortar) adapter
    std::vector<Teuchos::RCP<::ADAPTER::MortarVolCoupl>> volume_coupling_objects_;

    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra_;

    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra_;
  };
}  // namespace FS3I

#endif
