/*!----------------------------------------------------------------------
\file fps3i_partitioned.H
\brief General algorithmic routines for partitioned solution approaches
       to fluid-porous-structure-scalar-scalar interaction (FPS3I), that is,
       algorithmic routines not specifically related to partitioned
       solution approaches to one -or two-way-coupled problem
       configurations, respectively.

\level 3

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10364


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
}

namespace FSI
{
  class Monolithic;

  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }
}

namespace LINALG
{
  class MultiMapExtractor;
  class BlockSparseMatrixBase;
  class SparseMatrix;
  class Solver;
}


namespace FS3I
{
  class PartFPS3I : public FS3I_Base
  {

  public:

    //! constructor of base class for partitioned FPS3I
    PartFPS3I(const Epetra_Comm& comm);

    //! destructor of base class for partitioned FPS3I
    virtual ~PartFPS3I(){};

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
    void SetStructScatraSolution() ;

    //! return communicator
    const Epetra_Comm& Comm() const { return comm_; }


    /// extract fluid convective and structure convective velocities
    void ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector> >& vel,
                    std::vector<Teuchos::RCP<const Epetra_Vector> >& convel);

    /// extract Wall Shear Stresses at the interface
    void ExtractWSS(std::vector<Teuchos::RCP<const Epetra_Vector> >& wss);

    ///extracts mean concentration in membrane (interface)
    void ExtractMeanConcentration(std::vector<Teuchos::RCP<Epetra_Vector> >& MeanConc);

    ///extracts pressures at the interface
    void ExtractPressure(std::vector<Teuchos::RCP<const Epetra_Vector> >& pressure);

    ///provide velocities from FPSI subproblem for scatra subproblem
    void SetVelocityFields();

    ///provide wall shear stresses from FPSI subproblem for scatra subproblem
    void SetWallShearStresses();

    ///provide pressures from FPSI subproblem for scatra subproblem
    void SetPressureFields();

    ///provide information of mean concentration in membrane (interface)
    void SetMeanConcentration();

    ///provide displacements from FPSI subproblem for scatra subproblem
    void SetMeshDisp();

    ///Calculation of mean concentration in the membrane between fluid and poro
    Teuchos::RCP<Epetra_Vector> CalcMeanConcentration();


  protected:

    /// fpsi algorithm
    Teuchos::RCP<FPSI::Monolithic_Plain> fpsi_;


  private:

    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

  };
}

#endif
