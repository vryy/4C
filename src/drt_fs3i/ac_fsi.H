/*!----------------------------------------------------------------------
\file ac_fsi.H

\brief H-file associated with algorithmic routines for two-way coupled partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I). Specifically related version for multiscale approches

\date 2015-07-29

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            089/289-10364

\level 3
----------------------------------------------------------------------*/


#ifndef AC_FSI_H
#define AC_FSI_H

#include "fs3i_partitioned_1wc.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_inpar/inpar_fs3i.H"


namespace FS3I
{
  class ACFSI : public PartFS3I
  {
  public:

    ///constructor
    ACFSI(const Epetra_Comm& comm);

    virtual ~ACFSI(){};

    ///Read restart
    void ReadRestart();

    ///timeloop
    void Timeloop();

    ///timeloop for small time scales
    void SmallTimeScaleLoop();

    ///flag whether small time scale time loop should be finished
    bool SmallTimeScaleLoopNotFinished();

    ///Prepare small time scale time step
    void SmallTimeScalePrepareTimeStep();

    ///Prepare time step
    void PrepareTimeStep(){dserror("This function is not implemented! Use SmallTimeScalePrepareTimeStep() or LargeTimeScalePrepareTimeStep() instead!");};

    ///OuterLoop
    void SmallTimeScaleOuterLoop();

    ///OuterLoop for sequentially staggered FS3I scheme
    void SmallTimeScaleOuterLoopSequStagg();

    ///OuterLoop for iterative staggered FS3I scheme
    void SmallTimeScaleOuterLoopIterStagg();

    ///Do a single fsi step (including all subcycles)
    void DoFSIStep();

    void IsSmallTimeScalePeriodic();

    ///Decide if fsi problem is already periodic
    void IsFsiPeriodic();

    ///Do a standard fsi step
    void DoFSIStepStandard();

    ///Do a fsi step with subcycling
    void DoFSIStepSubcycled( const int subcyclingsteps);

    ///Get fsi solution from one period before
    void DoFSIStepPeriodic();

    ///Get step number of on cycle ago
    double GetStepOfPreviousPeriodAndPrepareReading(const int step_,const double time);

    ///Get filename in which the equivalent step of the last period is written
    std::string GetFileName(const int actstep);

    ///Set time and step in FSI and all subfields
    void SetTimeStepInFSI(const double time, const int step);

    ///Do a single scatra step
    void SmallTimeScaleDoScatraStep();

    ///Decide if fluid scatra problem is steady state
    void IsScatraPeriodic();

    ///Update and output the small time scale
    void SmallTimeScaleUpdateAndOutput();

    ///Write FSI output
    void FsiOutput();

    ///check convergence of scatra fields
    bool ScatraConvergenceCheck(const int itnum);

    ///Convergence check for iterative staggered FS3I scheme
    bool PartFs3iConvergenceCkeck( const int itnum);

    /// extract Wall Shear Stresses at the interface
    void ExtractWSS(std::vector<Teuchos::RCP<const Epetra_Vector> >& wss);

    //--------------------------------------------------------------------
    /// @name  control routines for the large time scale
    //--------------------------------------------------------------------

    ///timeloop for large time scales
    void LargeTimeScaleLoop();

    ///Prepare the large time scale loop
    void PrepareLargeTimeScaleLoop();

    ///Finish the large time scale loop
    void FinishLargeTimeScaleLoop();

    ///flag whether large time scale time loop should be finished
    bool LargeTimeScaleLoopNotFinished();

    ///Prepare large time scale time step
    void LargeTimeScalePrepareTimeStep();

    ///OuterLoop for sequentially staggered FS3I scheme
    void LargeTimeScaleOuterLoop();

    ///Do a large time scale structe scatra step
    void DoStructScatraStep();

    ///evaluate, solver and iteratively update structure scalar problem
    void StructScatraEvaluateSolveIterUpdate();

    ///check convergence of structure scatra field
    bool StructScatraConvergenceCheck(const int itnum);

    ///Do the structure scatra displacments need to update
    bool DoesGrowthNeedsUpdate();

    ///update the structure scatra displacments due to growth
    void LargeTimeScaleDoGrowthUpdate();

    ///Update and output the large time scale
    void LargeTimeScaleUpdateAndOutput();
    //@}

    ///Build map extractor which extracts the j-th dof
    std::vector<Teuchos::RCP<LINALG::MapExtractor> > BuildMapExtractor();

    ///Compare if two doubles are relatively zero
    bool IsRealtiveEqualTo(const double A, const double B, const double Ref=1.0);

    ///Compare if A mod B is relatively equal to zero
    bool ModuloIsRealtiveZero(const double A, const double B, const double Ref=1.0);

  private:

    /// structure increment vector
    const Teuchos::RCP<Epetra_Vector> structureincrement_;
    /// fluid increment vector
    const Teuchos::RCP<Epetra_Vector> fluidincrement_;
    /// ale increment vector
    const Teuchos::RCP<Epetra_Vector> aleincrement_;
    /// fluid phinp vector of the last period
    const Teuchos::RCP<Epetra_Vector> fluidphinp_lp_;
    /// structurephinp vector at the beginning of the large time scale loop
    const Teuchos::RCP<Epetra_Vector> structurephinp_bltsl_;
    /// structurephinp vector after the last growth update
    const Teuchos::RCP<Epetra_Vector> structurephinp_blgu_;

    ///WSS vector of the last period
    const Teuchos::RCP<Epetra_Vector> WallShearStress_lp_;

    ///time of one fsi period, e.g. time of a heart cycle
    const double fsiperiod_;

    ///time step for the large time scale problem
    const double dt_large_;

    ///flag iff fsi subproblem is periodic
    bool fsiisperiodic_;

    ///flag iff fluid scatra subproblem is periodic
    bool scatraisperiodic_;

    ///flag iff fluid scatra subproblem is periodic
    bool fsineedsupdate_;

    ///Extract the j-th out of numscal_ dof
    std::vector<Teuchos::RCP<LINALG::MapExtractor> > extractjthscalar_;
  };

}

#endif
