/*!----------------------------------------------------------------------
\file ac_fsi.H
\brief H-file associated with algorithmic routines for two-way coupled partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I). Specifically related version for multiscale approches

<pre>
Maintainers: Moritz Thon
             thon@mhpc.mw.tum.de
             089/289-10364
</pre>

*----------------------------------------------------------------------*/


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

    //timeloop
    void Timeloop();

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

    ///Decide if fsi problem is already periodic
    bool IsFsiPeriodic();

    ///Do a standard fsi step
    void DoFSIStepStandard();

    ///Do a fsi step with subcycling
    void DoFSIStepSubcycled();

    ///Get fsi solution from one period before
    void DoFSIStepPeriodic();

    ///Get filename in which the equivalent step of the last period is written
    std::string GetFileName();

    ///Set time and step in FSI and all subfields
    void SetTimeStepInFSI(const double time, const int step);

    ///Do a single scatra step
    void DoScatraStep();

    ///Update and output everything
    void UpdateAndOutput();

    ///Write FSI output
    void FsiOutput();

    ///check convergence of scatra fields
    bool ScatraConvergenceCheck(int itnum);

    ///Convergence check for iterative staggered FS3I scheme
    bool PartFs3iConvergenceCkeck( const int itnum);

    /// extract Wall Shear Stresses at the interface
    void ExtractWSS(std::vector<Teuchos::RCP<const Epetra_Vector> >& wss);

  private:

    /// scatra increment vector
    Teuchos::RCP<Epetra_Vector> structureincrement_;
    /// scatra increment vector
    Teuchos::RCP<Epetra_Vector> fluidincrement_;
    /// scatra increment vector
    Teuchos::RCP<Epetra_Vector> aleincrement_;
    ///WSS vector
    Teuchos::RCP<Epetra_Vector> WallShearStress_;

    ///do how many fsi time steps in comparison to the scalars time steps
    int fsiperssisteps_;

    ///time of one fsi period, e.g. time of a heart cycle
    double fsiperiod_;

    ///from which fsi period is the fsi problem assumed to be periodic
    int periodstillfsiisperiodic_;

    ///amount of fsi periods till we update the fsi problem
    int periodstoupdatefsi_;

    ///tolerance to define the the time when the fsi problem needs to be updated
    double fsiupdatetol_;
  };

}

#endif
