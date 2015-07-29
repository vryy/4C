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

    ///Prepare time step
    void PrepareTimeStep();

    ///OuterLoop
    void OuterLoop();

    ///OuterLoop for sequentially staggered FS3I scheme
    void OuterLoopSequStagg();

    ///OuterLoop for iterative staggered FS3I scheme
    void OuterLoopIterStagg();

    ///Do a single fsi step (including all subcycles)
    void DoFSIStep();

    void UpdateSmallTimeScale();

    ///Decide if fsi problem is already periodic
    void IsFsiPeriodic();

    ///Do a standard fsi step
    void DoFSIStepStandard();

    ///Do a fsi step with subcycling
    void DoFSIStepSubcycled();

    ///Get fsi solution from one period before
    void DoFSIStepPeriodic();

    ///Get step number of on cycle ago
    double GetStepOfPreviousPeriodAndPrepareReading(const int step_,const double time);

    ///Get filename in which the equivalent step of the last period is written
    std::string GetFileName(const int actstep);

    ///Set time and step in FSI and all subfields
    void SetTimeStepInFSI(const double time, const int step);

    ///Do a single scatra step
    void DoScatraStep();

    ///Decide if fluid scatra problem is steady state
    void IsScatraSteady();

    ///Update and output everything
    void UpdateAndOutput();

    ///Write FSI output
    void FsiOutput();

    ///check convergence of scatra fields
    bool ScatraConvergenceCheck(const int itnum);

    ///Convergence check for iterative staggered FS3I scheme
    bool PartFs3iConvergenceCkeck( const int itnum);

    /// extract Wall Shear Stresses at the interface
    void ExtractWSS(std::vector<Teuchos::RCP<const Epetra_Vector> >& wss);

    ///timeloop for large time scales
    void LargeTimeScaleLoop();

    ///flag whether large time scale time loop should be finished
    bool LargeTimeScaleLoopNotFinished();

  private:

    /// structure increment vector
    const Teuchos::RCP<Epetra_Vector> structureincrement_;
    /// fluid increment vector
    const Teuchos::RCP<Epetra_Vector> fluidincrement_;
    /// ale increment vector
    const Teuchos::RCP<Epetra_Vector> aleincrement_;
    /// phinp vector of the last period
    const Teuchos::RCP<Epetra_Vector> phinp_lp_;
    ///WSS vector of the last period
    const Teuchos::RCP<Epetra_Vector> WallShearStress_lp_;

    ///do how many fsi time steps in comparison to the scalars time steps
    const int fsiperssisteps_;

    ///time of one fsi period, e.g. time of a heart cycle
    const double fsiperiod_;

    ///flag iff fsi subproblem is periodic
    bool fsiisperiodic_;

    ///flag iff fluid scatra subproblem is periodic
    bool scatraisperiodic_;
  };

}

#endif
