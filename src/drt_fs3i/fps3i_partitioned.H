/*!----------------------------------------------------------------------
\file fps3i_partitioned.H
\brief General algorithmic routines for partitioned solution approaches
       to fluid-porous-structure-scalar-scalar interaction (FPS3I), that is,
       algorithmic routines not specifically related to partitioned
       solution approaches to one -or two-way-coupled problem
       configurations, respectively.

\level 3

\maintainer Christoph Schmidt


*----------------------------------------------------------------------*/


#ifndef PARTFPS3I_H
#define PARTFPS3I_H


#include "fs3i.H"


// forward declarations
namespace FPSI
{
  class Monolithic_Plain;
}

namespace ADAPTER
{
  class Coupling;
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
  class PartFPS3I : public FS3I_Base
  {
   public:
    //! constructor of base class for partitioned FPS3I
    PartFPS3I(const Epetra_Comm& comm);

    //! destructor of base class for partitioned FPS3I
    virtual ~PartFPS3I(){};

    //! initialize this class
    virtual void Init();

    //! setup this class
    virtual void Setup();

    //! time loop to be defined in inherited classes (structure depends on
    //! considered coupling, i.e. one-way or two-way)
    virtual void Timeloop() = 0;

    //! flag whether time loop should be finished
    bool NotFinished() { return step_ < numstep_ and time_ <= timemax_; };

    //! read and set fields needed for restart
    virtual void ReadRestart();

    /// redistribute FPS3I interface, if running on parallel
    void RedistributeInterface();

    //! set-up of FPSI and ScaTra systems
    virtual void SetupSystem();

    //! test results for individual fields
    virtual void TestResults(const Epetra_Comm& comm);

    //! evaluate ScaTra fields
    virtual void EvaluateScatraFields();

    //! information transfer FPSI -> ScaTra
    void SetFPSISolution();

    /// set scatra solution on structure field
    void SetStructScatraSolution();

    //! return communicator
    const Epetra_Comm& Comm() const { return comm_; }


    /// extract fluid convective and structure convective velocities
    void ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector>>& vel,
        std::vector<Teuchos::RCP<const Epetra_Vector>>& convel);

    /// extract Wall Shear Stresses at the interface
    void ExtractWSS(std::vector<Teuchos::RCP<const Epetra_Vector>>& wss);

    /// extracts pressures at the interface
    void ExtractPressure(std::vector<Teuchos::RCP<const Epetra_Vector>>& pressure);

    /// provide velocities from FPSI subproblem for scatra subproblem
    void SetVelocityFields();

    /// provide wall shear stresses from FPSI subproblem for scatra subproblem
    void SetWallShearStresses();

    /// provide pressures from FPSI subproblem for scatra subproblem
    void SetPressureFields();

    /// provide displacements from FPSI subproblem for scatra subproblem
    void SetMeshDisp();

   protected:
    /// fpsi algorithm
    Teuchos::RCP<FPSI::Monolithic_Plain> fpsi_;


   private:
    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

    /// scatra field on fluid
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra_;

    /// scatra field on structure
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra_;
  };
}  // namespace FS3I

#endif
