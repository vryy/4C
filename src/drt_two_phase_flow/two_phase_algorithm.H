/*----------------------------------------------------------------------*/
/*!
\file two_phase_algorithm.H

\brief Basis of all two phase flow (TPF) algorithms

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef TWO_PHASE_ALGORITHM_H
#define TWO_PHASE_ALGORITHM_H

#include <Epetra_Comm.h>
#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"
#include "../linalg/linalg_mapextractor.H"


/*!
\brief SCATRA namespace for two phase flow problems in BACI utilizing the levelset module

*/
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
    explicit Algorithm(
        const Epetra_Comm& comm,             ///< communicator
        const Teuchos::ParameterList& prbdyn,///< problem-specific parameter list
        const Teuchos::ParameterList& solverparams ///<solver-specific parameter list
        );

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// Two phase flow time loop
    void TimeLoop();

    ///Two phase flow stationary problem
    void SolveStationaryProblem();

    /// Restart TPF problem
    void Restart(int step, const bool restartscatrainput);

    /// Restart TPF problem with turbulent inflow
    void ReadInflowRestart(int restart);

    /// Test the results aquired
    void TestResults();

  protected:

    /// prepare time step
    void PrepareTimeStep();

    /// do (partitioned) outer iteration loop
    void OuterLoop();

    /// set scatra values required in fluid
    void SetScaTraValuesInFluid();

    /// convergence Check for present iteration step
    bool ConvergenceCheck(int itnum);

    /// update for next time step
    void TimeUpdate();

    /// write output
    void Output();

    /// get convective velocity at correct time
    void SetFluidValuesInScaTra(bool init);

    ///Write Output for initial fields
    void OutputInitialField();

    /// time-step length, maximum time and maximum number of steps
    double  dt_;
    double  maxtime_;
    int     stepmax_;

    /// (preliminary) maximum number of iterations and tolerance for outer iteration
    int     itmax_;
    double  ittol_;

    //flag for output specific options

    ///Masscenter
    bool    write_center_of_mass_;

    /// start of sampling period
    int     samstart_;

    /// flag for turbulent inflow
    bool    turbinflow_;

    /// number of inflow steps
    int     numinflowsteps_;

  private:

  };

} // namespace TWOPHASEFLOW

#endif // TWO_PHASE_ALGORITHM_H
