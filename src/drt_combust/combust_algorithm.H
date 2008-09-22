/*!----------------------------------------------------------------------*
\file combust_algorithm.H

\brief base combustion algorithm

	This algorithm handles the interaction between the different fields 
	involved in a combustion problem.
	So far there is only one algorithm. Others might follow in the future.
	This is why there are no further sub-routines of this algorithm.
	
	The interaction between the fluid field and the convection-diffusion 
	field (level-set/G-function) is controlled by this 	FGI algorithm 
	(Fluid-G-funcion-Interaction). It uses an implicit time integration 
	scheme.
	
<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef COMBUST_ALGORITHM_H
#define COMBUST_ALGORITHM_H

#include <Epetra_Comm.h>

#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"

/*!
\brief COMBUST: namespace of the combustion module in BACI

*/
namespace COMBUST
{

  /// COMBUST algorithm base
  /*!

    Base class of COMBUST algorithms. Derives from the Scatra-Fluid Coupling Algorithm. Therefore it
    has a Fluid Base Algorithm and a Scatra Base Algorithm.
    There can (and will) be different subclasses that implement different coupling schemes.

    \author henke
    \date 06/08
   */
  class Algorithm : public ADAPTER::ScaTraFluidCouplingAlgorithm
  {
  public:

    /// constructor using a Epetra_Comm
    explicit Algorithm(Epetra_Comm& comm, Teuchos::ParameterList& combustdyn);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// time loop of algorithm for dynamic combustion problem
    void TimeLoop();

    /// algorithm for static combustion problem
    void SolveStationaryProblem();

  protected:

    //! Options to reinitialize the G-function
    enum ReinitializationAction
    {
      initialize,
      compute_signeddistancefunction
    };

    /// FG iteration number
    int FGiter() const { return fgiter_; }

    /// reinitialize G-function
    void ReinitializeGfunc(ReinitializationAction action);

    /// initialize g-function by scalar function field
    void InitializeGFunc();

    /// build signed distance function
    void SignedDistFunc();

    /// tests if the FGI iteration has converged
    bool NotConvergedFGI();

    /// prepare a new time step
    void PrepareTimeStep();

    /// start a new Fluid-G-function iteration
    void PrepareFGIteration();

    /// solve Navier-Stokes equations for current time step
    void DoFluidField();

    /// solve G-function equation for current time step
    void DoGfuncField();

    /// update fields in FGI iteration;
    void UpdateFGIteration();

    /// take current results for converged and save for next time step
    void UpdateTimeStep();

    /// write output
    void Output();

  private:

    //! @name Time stepping variables
    //int step_;
    //int stepmax_;
    //double time_;
    //double timemax_;
    //double dt_;
    //@}
    
    //! @name Fluid-G-function iteration variables
    // Fluid-G-function iteration counter
    int fgiter_;
    // maximum number of Fluid-G-function iterations
    int fgitermax_;
    // L2-norm of velocity in Fluid-G-function iteration
    double fgvelnormL2_;
    // L2-norm of g-function in Fluid-G-function iteration
    double fggfuncnormL2_;
    //@}

    /// options for reinitialization of G-function;
    ReinitializationAction reinitializationaction_;

  };

} // namespace COMBUST

#endif // #ifndef COMBUST_ALGORITHM_H
#endif // #ifdef CCADISCRET
