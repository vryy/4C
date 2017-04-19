/*----------------------------------------------------------------------*/
/*!
\file xfluid_levelset_coupling_algorithm.H

\brief Base algorithm of Xfluid utilizing level set for cut

\level 3

<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915236
</pre>
 */
/*----------------------------------------------------------------------*/

#ifndef XFLUID_LEVELSET_ALGORITHM_H
#define XFLUID_LEVELSET_ALGORITHM_H

#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"
#include "../linalg/linalg_mapextractor.H"

//Needed for input parameter
#include "../drt_inpar/inpar_twophase.H"


namespace XFLUIDLEVELSET
{

  /// Xfluid level set coupling
  /*!

    Base class of xfluid and level set coupling. Derives from FluidBaseAlgorithm
    and ScatraBaseAlgorithm.


    \author winter
    \date 09/14
   */
  class Algorithm : public ADAPTER::ScaTraFluidCouplingAlgorithm
  {
  public:

    /// constructor
    explicit Algorithm(
        const Epetra_Comm& comm,             ///< communicator
        const Teuchos::ParameterList& prbdyn,///< problem-specific parameter list
        const Teuchos::ParameterList& solverparams, ///<solver-specific parameter list
        const std::map<std::string, int>& dofset_couling_map ///< key=meaning of a dofset in a dis, nodal dofset number
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
    \date 08/16
    \author rauch  */
    virtual void Init(
        const Teuchos::ParameterList&   prbdyn,              ///< parameter list for global problem
        const Teuchos::ParameterList&   scatradyn,           ///< parameter list for scalar transport subproblem
        const Teuchos::ParameterList&   solverparams,        ///< parameter list for scalar transport solver
        const std::string&              disname = "scatra",  ///< name of scalar transport discretization
        const bool                      isale = false        ///< ALE flag
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
    \date 08/16
    \author rauch  */
    virtual void Setup();

    /// Two phase flow time loop
    void TimeLoop();

    /// Two phase flow stationary problem
    void SolveStationaryProblem();

    /// Test the results aquired
    void TestResults();

    /// Restart
    void Restart(int restart);

  protected:

    virtual void DoAlgorithmSpecificInit();

    /// convergence Check for present iteration step
    bool ConvergenceCheck(int itnum);
    /// get increment of Fluid solution between time step i and i+1
    void GetOuterLoopIncFluid(double& fsvelincnorm, double& fspressincnorm, int itnum);
    /// get increment of ScaTra solution between time step i and i+1
    void GetOuterLoopIncScaTra(double& fsphiincnorm, int itnum);

    /// fluid-level set loop
    void OuterLoop();
    /// set values from ScaTra and solve the fluid field
    void DoFluidField();
    /// set values from Fluid and solve the scatra field
    void DoScaTraField();
    /// print string to std cout
    void PrintToScreen(const std::string str);

    //Communication between fields:
    /// set scatra values required in fluid
    void SetScaTraValuesInFluid();
    /// set fluid values required in scatra
    void SetFluidValuesInScaTra(bool init);

    /// compute geometric quantities (normalvector field and curvatures) from level set field
    void ComputeGeometricQuantities(
        const bool require_smoothedgradphi,
        const bool require_nodalcurvature,
        const Teuchos::RCP<const Epetra_Vector> & phinp,
        const INPAR::TWOPHASE::SmoothGradPhi smoothedgradphitype,
        Teuchos::RCP<Epetra_MultiVector> & smoothedgradphi,
        Teuchos::RCP<Epetra_Vector> & nodalcurvature
    );

    /// copy geometric quantities from scatra to cutter discretization stored in the conditionmanager
    void CopyGeometricQuantities(
        const bool require_smoothedgradphi,
        const bool require_nodalcurvature,
        const Teuchos::RCP<const Epetra_Vector> & scatra_phinp,
        const Teuchos::RCP<const Epetra_MultiVector> & scatra_smoothedgradphi,
        const Teuchos::RCP<const Epetra_Vector>      & scatra_nodalcurvature,
        const Teuchos::RCP<Epetra_Vector> & fluid_phinp,
        const Teuchos::RCP<Epetra_MultiVector> & fluid_smoothedgradphi,
        const Teuchos::RCP<Epetra_Vector>      & fluid_nodalcurvature
    );


    /// update time
    void TimeUpdate();
    /// prepare time step
    void PrepareTimeStep();
    /// prepare outer iteration
    void PrepareOuterIteration();

    /// synchronicity check between algorithm and fields
    void SynchronicityTimeCheck();

    /// write output
    void Output();
    /// write output of initial field
    void OutputInitialField();


    /// Get smoothed LevelSet gradient
    Teuchos::RCP<Epetra_MultiVector> GetSmoothedLevelSetGradient(const Teuchos::RCP<const Epetra_Vector> & phinp,
                                                                 INPAR::TWOPHASE::SmoothGradPhi        smoothedgradphi);
  private:

    /// time-step length, maximum time and maximum number of steps
    double  dt_;
    double  maxtime_;
    int     stepmax_;

    /// (preliminary) maximum number of iterations and tolerance for outer iteration
    int     itmax_;
    double  ittol_;

    /// flag for output specific options
    int upres_;

    /// masscenter output
    bool    write_center_of_mass_;

    /// Two Phase Specific parameter
    INPAR::TWOPHASE::SmoothGradPhi        smoothedgradphitype_;
    bool                                  scalesmoothedgradients_;

    /// Outer-loop convergence variables:
    Teuchos::RCP<Epetra_Vector> velnpi_;
    Teuchos::RCP<Epetra_Vector> phinpi_;
    /// vectors for visualization of outerloop iteration progress.
    std::vector<double> fsvelincnorm_;
    std::vector<double> fspressincnorm_;
    std::vector<double> fsphiincnorm_;

    /// problem dynamic parameters
    Teuchos::ParameterList prbdyn_;

    std::map<std::string, int> dofset_coupling_map_;

  };

} // namespace XFLUIDLEVELSET

#endif // XFLUID_LEVELSET_ALGORITHM_H
