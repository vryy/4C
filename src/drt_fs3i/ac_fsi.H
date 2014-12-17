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

    ///Update and output everything
    void UpdateAndOutput();

    ///Do a single fsi step (including all subcycles)
    void DoFSIStep();

    ///Do a single scatra step
    void DoScatraStep();

    ///Convergence check for iterative staggered FS3I scheme
    bool PartFs3iConvergenceCkeck( const int itnum);

    ///check convergence of scatra fields
    bool ScatraConvergenceCheck(int itnum);

    ///Read restart
    void ReadRestart();

  private:

    ///maximum number of iterationsfor outer iteration
    int     itmax_;
    ///tolerance for outer iteration
    double  ittol_;


    /// scatra increment vector
    Teuchos::RCP<Epetra_Vector> structureincrement_;
    /// scatra increment vector
    Teuchos::RCP<Epetra_Vector> fluidincrement_;
    /// scatra increment vector
    Teuchos::RCP<Epetra_Vector> aleincrement_;

    // coupling of fluid-structure interaction and scalar-scalar interaction problems
    INPAR::FS3I::SolutionSchemeOverFields couplingalgo_;

    //do how many fsi time steps in comparison to the scalars time steps
    int fsiperssisteps_;
  };

}

#endif
