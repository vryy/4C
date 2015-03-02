/*----------------------------------------------------------------------*/
/*!
\file xfluid_levelset_coupling_algorithm.H

\brief Base algorithm of Xfluid utilizing level set for cut

<pre>
Maintainer: Magnus Winter
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
        const Teuchos::ParameterList& solverparams ///<solver-specific parameter list
        );

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// Two phase flow time loop
    void TimeLoop();

    /// Two phase flow stationary problem
    void SolveStationaryProblem();

    /// Test the results aquired
    void TestResults();

    /// Restart
    void Restart(int restart);

  protected:

    /// convergence Check for present iteration step
    bool ConvergenceCheck(int itnum);


    /// fluid-level set loop
    void OuterLoop();
    /// set values from ScaTra and solve the fluid field
    void DoFluidField();
    /// set values from Fluid and solve the scatra field
    void DoScaTraField();


    //Communication between fields:
    /// set scatra values required in fluid
    void SetScaTraValuesInFluid();
    /// set scatra values required in fluid
    void SetFluidValuesInScaTra(bool init);


    /// update time
    void TimeUpdate();
    /// prepare time step
    void PrepareTimeStep();
    /// prepare outer iteration
    void PrepareOuterIteration();


    /// write output
    void Output();
    /// write output of initial field
    void OutputInitialField();


    /// set problem specific parameters
    void SetProblemSpecificParameters(const Teuchos::ParameterList& prbdyn);


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
    INPAR::TWOPHASE::SurfaceTensionApprox surftensapprox_;

    /// Outer-loop convergence variables:
    Teuchos::RCP<Epetra_Vector> velnpi_;
    Teuchos::RCP<Epetra_Vector> phinpi_;
    /// vectors for visualization of outerloop iteration progress.
    std::vector<double> fsvelincnorm_;
    std::vector<double> fspressincnorm_;
    std::vector<double> fsphiincnorm_;

  };

} // namespace XFLUIDLEVELSET

#endif // XFLUID_LEVELSET_ALGORITHM_H
