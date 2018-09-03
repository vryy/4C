/*----------------------------------------------------------------------*/
/*!
\file adapter_algorithmbase.H

\brief Base algorithm for all kinds of coupled problems

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_ALGORITHMBASE_H
#define ADAPTER_ALGORITHMBASE_H

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

namespace ADAPTER
{
  /// base class for algorithms
  /*!

    Here we just keep the time variables as it is always the same.

    \author u.kue
    \date 07/09
   */
  class AlgorithmBase
  {
   public:
    /// create using a Epetra_Comm
    explicit AlgorithmBase(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /// virtual destruction
    virtual ~AlgorithmBase() {}

    /// read restart data
    virtual void ReadRestart(int step) = 0;

    /// read restart data from time
    virtual void ReadRestartfromTime(double time);

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

    /// tests if there are more time steps to do
    bool NotFinished() const { return step_ < nstep_ and time_ + 1e-8 * dt_ < maxtime_; }

    /// time step size
    double Dt() const { return dt_; }

    /// set time step size
    void SetDt(const double stepsize) { dt_ = stepsize; }

    /// current time
    double Time() const { return time_; }

    /// current time step number
    int Step() const { return step_; }

    /// set new time step in ReadRestart() or in a potential outer control
    void SetTimeStep(const double time,  ///< physical time to set
        const int step                   ///< time step number to set
    );

   protected:
    /// total number of time steps
    int NStep() const { return nstep_; }

    /// maximum simulation time
    double MaxTime() const { return maxtime_; }

    //! @name Time loop building blocks

    /// increment time and step value
    void IncrementTimeAndStep()
    {
      step_ += 1;
      time_ += dt_;
    }

    /// start a new time step
    virtual void PrepareTimeStep() {}

    /// take current results for converged and save for next time step
    virtual void Update() {}

    /// write output
    virtual void Output() {}

    /// set method name for screen output
    void SetMethod(std::string method) { method_ = method; }

    //@}

    /// print time step header
    void PrintHeader();

    /// return printscreen_
    int PrintScreenEvry() { return printscreen_; }

   private:
    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

    /// method name prepared for output
    std::string method_;

    /// print infos to standard out every printscreen_ steps
    int printscreen_;

    //! @name Time stepping variables

    /// current time step number
    int step_;

    /// total number of time steps
    int nstep_;

    /// current physical time at end of time step, i.e. \f$t_{n+1}\f$
    double time_;

    /// maximum simulation time
    double maxtime_;

    /// current time step size
    double dt_;

    //@}
  };
}  // namespace ADAPTER

#endif
