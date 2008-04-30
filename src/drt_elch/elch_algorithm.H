/*----------------------------------------------------------------------*/
/*!
\file elch_algorithm.H

\brief Basis of all ELCH algorithms

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ELCH_ALGORITHM_H
#define ELCH_ALGORITHM_H

#include <Epetra_Comm.h>

#include "../drt_adapter/adapter_fluid_base_algorithm.H"
#include "../drt_adapter/adapter_condif_base_algorithm.H"


/*!
\brief ELCH: namespace of the electrochemistry module in BACI

*/
namespace ELCH
{

  /// ELCH algorithm base
  /*!

    Base class of ELCH algorithms. Derives from FluidBaseAlgorithm
    and ConDifBaseAlgorithm.
    There can (and will) be different subclasses that implement
    different coupling schemes.

    \author gjb
    \date 03/08
   */
  class Algorithm : public ADAPTER::FluidBaseAlgorithm,
                    public ADAPTER::ConDifBaseAlgorithm
  {
  public:

    /// constructor using a Epetra_Comm
    explicit Algorithm(Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// outer level ELCH time loop
    void TimeLoop();

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

  protected:

    /// time step size
    double Dt() const { return dt_; }

    /// time step number
    int Step() const { return step_; }

    /// tests if there are more time steps to do
    bool NotFinished() {return ((step_ < nstep_) and (time_ <= maxtime_));}

    /// start a new time step
    void PrepareTimeStep();
    
    /// solve Navier-Stokes equations for current time step
    void DoFluidStep();
    
    /// solve transport equations for current time step
    void DoTransportStep();

    /// take current results for converged and save for next time step
    void Update();

    /// write output
    void Output();

  private:

    /// communication (mainly for screen output)
    Epetra_Comm& comm_;

    /// actual velocity field at time n+1
    RCP<const Epetra_Vector>    velocitynp_;

    //! @name Time stepping variables
    int step_;
    int nstep_;
    double time_;
    double maxtime_;
    double dt_;
    //@}

  };

} // namespace ELCH

#endif // ELCH_ALGORITHM_H
#endif // CCADISCRET
