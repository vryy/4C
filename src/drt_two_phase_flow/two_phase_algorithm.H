/*----------------------------------------------------------------------*/
/*!

\brief Basis of all two phase flow (TPF) algorithms

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915236

*/
/*----------------------------------------------------------------------*/

#ifndef TWO_PHASE_ALGORITHM_H
#define TWO_PHASE_ALGORITHM_H

#include <Epetra_Comm.h>
#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"
#include "../linalg/linalg_mapextractor.H"

// Needed for input parameter
#include "../drt_inpar/inpar_twophase.H"

namespace TWOPHASEFLOW
{
  /// Two phase flow base algorithm
  /*!

    Base class of two phase flow algorithms. Derives from FluidBaseAlgorithm
    and ScatraBaseAlgorithm.
    There can (and will) be different subclasses that implement
    different coupling schemes.

    \author winter
    \date 05/14
   */
  class Algorithm : public ADAPTER::ScaTraFluidCouplingAlgorithm
  {
   public:
    /// constructor
    explicit Algorithm(const Epetra_Comm& comm,     ///< communicator
        const Teuchos::ParameterList& prbdyn,       ///< problem-specific parameter list
        const Teuchos::ParameterList& solverparams  ///< solver-specific parameter list
    );

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in Init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in Init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref Setup().

    \warning none
    \return void
    \date 09/16
    \author rauch  */
    virtual void Init(const Teuchos::ParameterList& prbdyn,  ///< parameter list for global problem
        const Teuchos::ParameterList&
            scatradyn,  ///< parameter list for scalar transport subproblem
        const Teuchos::ParameterList& solverparams,  ///< parameter list for scalar transport solver
        const std::string& disname = "scatra",       ///< name of scalar transport discretization
        const bool isale = false                     ///< ALE flag
    );

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 09/16
    \author rauch  */
    virtual void Setup();

    /// Two phase flow time loop
    void TimeLoop();

    /// Two phase flow stationary problem
    void SolveStationaryProblem();

    /// Restart TPF problem
    void Restart(int step, const bool restartscatrainput);

    /// Restart TPF problem with turbulent inflow
    void ReadInflowRestart(int restart);

    /// Test the results aquired
    void TestResults();

   protected:
    /// do (partitioned) outer iteration loop
    void OuterLoop();
    /// set values from ScaTra and solve the fluid field
    void DoFluidField();
    /// set values from Fluid and solve the scatra field
    void DoScaTraField();

    /// convergence Check for present iteration step
    bool ConvergenceCheck(int itnum);
    /// get increment of Fluid solution between time step i and i+1
    void GetOuterLoopIncFluid(double& fsvelincnorm, double& fspressincnorm, int itnum);
    /// get increment of ScaTra solution between time step i and i+1
    void GetOuterLoopIncScaTra(double& fsphiincnorm, int itnum);

    /// update for next time step
    void TimeUpdate();
    /// prepare time step
    void PrepareTimeStep();
    /// prepare outer iteration
    void PrepareOuterIteration();

    /// write output
    void Output();
    /// Write Output for initial fields
    void OutputInitialField();

    /// get convective velocity at correct time
    void SetFluidValuesInScaTra(bool init);
    /// set scatra values required in fluid
    void SetScaTraValuesInFluid();

    /// Get smoothed LevelSet gradient
    Teuchos::RCP<Epetra_MultiVector> GetSmoothedLevelSetGradient(
        const Teuchos::RCP<const Epetra_Vector>& phinp);

   private:
    /// time-step length, maximum time and maximum number of steps
    double dt_;
    double maxtime_;
    int stepmax_;

    /// (preliminary) maximum number of iterations and tolerance for outer iteration
    int itmax_;
    double ittol_;

    // flag for output specific options

    /// Masscenter
    bool write_center_of_mass_;

    /// flag for turbulent inflow
    bool turbinflow_;

    /// number of inflow steps
    int numinflowsteps_;


    /// Outer-loop convergence variables:
    Teuchos::RCP<Epetra_Vector> velnpi_;
    Teuchos::RCP<Epetra_Vector> phinpi_;
    /// vectors for visualization of outerloop iteration progress.
    std::vector<double> fsvelincnorm_;
    std::vector<double> fspressincnorm_;
    std::vector<double> fsphiincnorm_;

    Teuchos::ParameterList prbdyn_;

    /// Two Phase Specific parameter
    INPAR::TWOPHASE::SmoothGradPhi smoothedgradphitype_;
    bool scalesmoothedgradients_;
  };

}  // namespace TWOPHASEFLOW

#endif  // TWO_PHASE_ALGORITHM_H
