/*!----------------------------------------------------------------------

\file ac_fsi.H

\brief H-file associated with algorithmic routines for two-way coupled partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I). Specifically related version for multiscale approches

\level 3

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10364

----------------------------------------------------------------------*/


#ifndef AC_FSI_H
#define AC_FSI_H

#include "fs3i_partitioned_1wc.H"
#include "../drt_lib/drt_dserror.H"

// forward declarations
namespace LINALG
{
  class MapExtractor;
}

namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace IO

namespace FS3I
{
  class MeanManager;

  /*!
   * What does the problem type Atherosclerosis_Fluid_Structure_Interaction?
   * Short answer: cool stuff!
   * And here is the long answer:
   * Its doing a multiscale (in time) approach with an full FS3I simulation at the small time scale
   * (seconds) and a scalar transport simulation at the larger time scale (days). The solving
   * strategy is as follows: We start with the full small time scale FS3i simulation (including
   * fluid Windkessel and WSS permeability). After each FSI cycle we check if the FSI problem is
   * periodic, by looking if the Windkessel does produce peridodic results. Afterwards we continue
   * the small time scale but do not solve the fsi subproblem anymore. Instead we peridically repeat
   * it by calling suitable restarts. When the fluid scatra subproblem gets periodic at the FS3I
   * interface we stop the small time scale and switch to the large time scale Now we higher dt_ and
   * only solve the structure scatra problem. We thereby use the WSS and interface concentrations of
   * the small time scale. Each time when there as been 'created' enough growth inducing mass we do
   * a growth update. If we have finally grew to much, we go back to the small time scale. And so
   * on, and so on,...
   */
  class ACFSI : public PartFS3I
  {
   public:
    /// constructor
    ACFSI(const Epetra_Comm& comm);

    /// destructor
    virtual ~ACFSI(){};

    /// initialize this class
    virtual void Init();

    /// setup this class
    virtual void Setup();

    /// Read restart
    void ReadRestart();

    /// timeloop
    void Timeloop();

    /// timeloop for small time scales
    void SmallTimeScaleLoop();

    /// flag whether small time scale time loop should be finished
    bool SmallTimeScaleLoopNotFinished();

    /// Prepare small time scale time step
    void SmallTimeScalePrepareTimeStep();

    /// Prepare time step
    void PrepareTimeStep()
    {
      dserror(
          "This function is not implemented! Use SmallTimeScalePrepareTimeStep() or "
          "LargeTimeScalePrepareTimeStep() instead!");
    };

    /// OuterLoop
    void SmallTimeScaleOuterLoop();

    /// OuterLoop for sequentially staggered FS3I scheme
    void SmallTimeScaleOuterLoopSequStagg();

    /// OuterLoop for iterative staggered FS3I scheme
    void SmallTimeScaleOuterLoopIterStagg();

    /// Do a single fsi step (including all subcycles)
    void DoFSIStep();

    void IsSmallTimeScalePeriodic();

    /// Decide if fsi problem is already periodic
    void IsFsiPeriodic();

    /// provide wall shear stresses from FS3I subproblem for scatra subproblem
    void SetWallShearStresses() const;

    /// Decide if fluid scatra problem is periodic
    void IsScatraPeriodic();

    /// Do a standard fsi step
    void DoFSIStepStandard();

    /// Do a fsi step with subcycling
    void DoFSIStepSubcycled(const int subcyclingsteps);

    /// Get fsi solution from one period before
    void DoFSIStepPeriodic();

    /// Get step number of on cycle ago
    double GetStepOfOnePeriodAgoAndPrepareReading(const int step_, const double time);

    /// Get step number of the beginning of this cycle
    double GetStepOfBeginnOfThisPeriodAndPrepareReading(
        const int actstep, const double acttime, const double dt);

    /// Get filename in which the equivalent step of the last period is written
    std::string GetFileName(const int step);

    /// Set time and step in FSI and all subfields
    void SetTimeAndStepInFSI(const double time, const int step);

    /// Do a single scatra step
    void SmallTimeScaleDoScatraStep();

    /// Update and output the small time scale
    void SmallTimeScaleUpdateAndOutput();

    /// Write FSI output
    void FsiOutput();

    /// check convergence of scatra fields
    bool ScatraConvergenceCheck(const int itnum);

    /// Convergence check for iterative staggered FS3I scheme
    bool PartFs3iConvergenceCkeck(const int itnum);

    //--------------------------------------------------------------------
    /// @name  control routines for the large time scale
    //--------------------------------------------------------------------

    /// timeloop for large time scales
    void LargeTimeScaleLoop();

    /// Prepare the large time scale loop
    void PrepareLargeTimeScaleLoop();

    /// Set mean wall shear stresses in scatra fields
    void SetMeanWallShearStresses() const;

    /// Set mean concentration of the fluid scatra field
    void SetMeanFluidScatraConcentration();

    /// Set zero velocity field in scatra fields
    void SetZeroVelocityField();

    /// Evaluate surface permeability condition for structure scatra field
    void EvaluateithScatraSurfacePermeability(const int i  ///< id of scalar to evaluate
    );

    /// Finish the large time scale loop
    void FinishLargeTimeScaleLoop();

    /// flag whether large time scale time loop should be finished
    bool LargeTimeScaleLoopNotFinished();

    /// Prepare large time scale time step
    void LargeTimeScalePrepareTimeStep();

    /// OuterLoop for sequentially staggered FS3I scheme
    void LargeTimeScaleOuterLoop();

    /// Do a large time scale structe scatra step
    void DoStructScatraStep();

    /// evaluate, solver and iteratively update structure scalar problem
    void StructScatraEvaluateSolveIterUpdate();

    /// check convergence of structure scatra field
    bool StructScatraConvergenceCheck(const int itnum  ///< current iteration number
    );

    /// Do the structure scatra displacments need to update
    bool DoesGrowthNeedsUpdate();

    /// update the structure scatra displacments due to growth
    void LargeTimeScaleDoGrowthUpdate();

    /// OuterLoop for large time scale iterative staggered FS3I scheme
    void LargeTimeScaleOuterLoopIterStagg();

    /// set mean FSI values in scatra fields (only to be used in large time scale!!)
    void LargeTimeScaleSetFSISolution();

    /// Update and output the large time scale
    void LargeTimeScaleUpdateAndOutput();
    //@}

    /// Build map extractor which extracts the j-th dof
    std::vector<Teuchos::RCP<LINALG::MapExtractor>> BuildMapExtractor();

    /// optional safety check for times and dt's of all fields
    void CheckIfTimesAndStepsAndDtsMatch();

    /// Compare if two doubles are relatively zero
    bool IsRealtiveEqualTo(const double A,  ///< first value
        const double B,                     ///< second value
        const double Ref = 1.0              ///< reference value
    );

    /// Compare if A mod B is relatively equal to zero
    bool ModuloIsRealtiveZero(const double A,  ///< value to mod
        const double B,                        ///< mod value
        const double Ref = 1.0                 ///< reference value
    );

   private:
    /// structure increment vector
    Teuchos::RCP<Epetra_Vector> structureincrement_;
    /// fluid increment vector
    Teuchos::RCP<Epetra_Vector> fluidincrement_;
    /// ale increment vector
    Teuchos::RCP<Epetra_Vector> aleincrement_;
    /// mean fluid phinp vector of the last period
    Teuchos::RCP<Epetra_Vector> fluidphinp_lp_;
    /// structurephinp vector at the beginning of the large time scale loop
    Teuchos::RCP<Epetra_Vector> structurephinp_blts_;
    /// growth update counter
    int growth_updates_counter_;

    /// mean WSS vector of the last period
    Teuchos::RCP<Epetra_Vector> WallShearStress_lp_;

    /// time of one fsi period, e.g. time of a heart cycle
    const double fsiperiod_;

    /// time step for the large time scale problem
    const double dt_large_;

    /// flag iff fsi subproblem is periodic
    bool fsiisperiodic_;

    /// flag iff fluid scatra subproblem is periodic
    bool scatraisperiodic_;

    /// flag iff fluid scatra subproblem is periodic
    bool fsineedsupdate_;

    /// Extract the j-th out of numscal_ dof
    std::vector<Teuchos::RCP<LINALG::MapExtractor>> extractjthstructscalar_;

    /// pointer to mean manager object
    Teuchos::RCP<FS3I::MeanManager> meanmanager_;
  };

  class MeanManager
  {
   public:
    /// constructor
    MeanManager(const Epetra_Map& wssmap, const Epetra_Map& phimap, const Epetra_Map& pressuremap);

    /// destructor
    virtual ~MeanManager(){};

    /// add value into the mean manager
    void AddValue(
        const std::string type, const Teuchos::RCP<const Epetra_Vector> value, const double dt);

    /// reset mean manager
    void Reset();

    /// get some mean value
    Teuchos::RCP<const Epetra_Vector> GetMeanValue(const std::string type) const;

    /// Write restart of mean manager
    void WriteRestart(Teuchos::RCP<IO::DiscretizationWriter> fluidwriter) const;

    /// Read restart of mean manager
    void ReadRestart(IO::DiscretizationReader& fluidreader);

   private:
    /// weighted sum of all prior wall shear stresses
    Teuchos::RCP<Epetra_Vector> SumWss_;
    /// weighted sum of all prior concentrations
    Teuchos::RCP<Epetra_Vector> SumPhi_;
    /// weighted sum of all prior pressures
    Teuchos::RCP<Epetra_Vector> SumPres_;

    double SumDtWss_;
    double SumDtPhi_;
    double SumDtPres_;
  };
}  // namespace FS3I

#endif
