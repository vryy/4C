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
    
    /// FG iteration number
    int FGiter() const { return fgiter_; }

    /// tests if there are more time steps to do
    bool NotFinished() { return ((step_ < stepmax_) and (time_ <= timemax_)); }
    
    /// reinitialize G-function
    void ReinitializeGfunc();

    /// build signed distance function
    void SignedDistFunc();
    
    /// tests if the FGI iteration has converged
    bool NotConvergedFGI();
    
    /// start a new time step
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

    /// communication (mainly for screen output)
    Epetra_Comm& comm_;

    /// actual velocity field at time n+1
    RCP<const Epetra_Vector>    velocitynp_;

    //! @name Time stepping variables
    int step_;
    int stepmax_;
    double time_;
    double timemax_;
    double dt_;
    //@}
    
    //! @name FG iteration variables
    // FG iteration counter
    int fgiter_;
    // maximum number of FG iterations
    int fgitermax_;
    // L2-norm of velocity in FG iteration
    double fgvelnormL2_;
    // L2-norm of g-function in FG iteration
    double fggfuncnormL2_;
    // ?
    bool notconverged_;
    //@}
    
    // this is a parameter that should be read from the parameter list
    int reinitialize_action;
    // this is a parameter that should be contained in the validparameters list as an enum.
    int signeddistfunc;
  };

} // namespace COMBUST

#endif // #ifndef COMBUST_ALGORITHM_H
#endif // #ifdef CCADISCRET
