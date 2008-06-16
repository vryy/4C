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

#include "../drt_adapter/adapter_fluid_base_algorithm.H"
#include "../drt_adapter/adapter_condif_base_algorithm.H"


/*!
\brief COMBUST: namespace of the combustion module in BACI

*/
namespace COMBUST
{

  /// COMBUST algorithm base
  /*!

    Base class of COMBUST algorithms. Derives from FluidBaseAlgorithm
    and ConDifBaseAlgorithm.
    There can (and will) be different subclasses that implement
    different coupling schemes.

    \author henke
    \date 06/08
   */
  class Algorithm : public ADAPTER::FluidBaseAlgorithm,
                    public ADAPTER::ConDifBaseAlgorithm
  {
  public:

    /// constructor using a Epetra_Comm
    explicit Algorithm(Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// outer level COMBUST time loop
    void TimeLoop();

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

  protected:

    /// time step size
    double Dt() const { return dt_; }

    /// time step number
    int Step() const { return step_; }

    /// tests if there are more time steps to do
    bool NotFinished() { return ((step_ < numstep_) and (time_ <= maxtime_)); }
    
    /// reinitialize G-function
    void ReinitializeGfunc();

    /// tests if the FGI iteration has converged
    bool NotConverged();
    
    /// start a new time step
    void PrepareTimeStep();
    
    /// start a new Fluid-G-function iteration
    void PrepareFGIteration();
    
    /// solve Navier-Stokes equations for current time step
    void DoFluidStep();
    
    /// solve G-function equation for current time step
    void DoGfuncStep();
    
    /// update fields in FGI iteration;
    void UpdateFGI();

    /// take current results for converged and save for next time step
    void UpdateTime();

    /// write output
    void Output();

  private:

    /// communication (mainly for screen output)
    Epetra_Comm& comm_;

    /// actual velocity field at time n+1
    RCP<const Epetra_Vector>    velocitynp_;

    //! @name Time stepping variables
    int step_;
    int numstep_;
    double time_;
    double maxtime_;
    double dt_;
    //@}
    
    // variables for FGI iteration
    int fgiter;
    int maxfgiter;
    bool notconverged;
    // variables for Fluid iterations
    int fiter;
    // variables for G-function iterations
    int giter;

  };

} // namespace COMBUST

#endif // #ifndef COMBUST_ALGORITHM_H
#endif // #ifdef CCADISCRET
