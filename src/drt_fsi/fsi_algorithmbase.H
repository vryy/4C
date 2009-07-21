#ifndef FSI_ALGORITHMBASE_H
#define FSI_ALGORITHMBASE_H

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

namespace FSI
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
    explicit AlgorithmBase(Epetra_Comm& comm,
                           const Teuchos::ParameterList& timeparams);

    /// virtual destruction
    virtual ~AlgorithmBase() {}

    /// read restart data
    virtual void ReadRestart(int step) = 0;

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

  protected:

    /// time step size
    double Dt() const { return dt_; }

    /// step number
    int Step() const { return step_; }

    /// total number of time steps
    int NStep() const { return nstep_;}

    /// maximum simulation time
    double MaxTime() const { return maxtime_;}

    /// current time
    double Time() const {return time_;}

    //! @name Time loop building blocks

    /// tests if there are more time steps to do
    bool NotFinished() { return step_ < nstep_ and time_+1e-10*dt_ < maxtime_; }

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

    /// set new time step in ReadRestart()
    void SetTimeStep(double time, int step);

    /// set method name for screen output
    void SetMethod(std::string method) { method_ = method; }

    //@}

    /// print time step header
    void PrintHeader();

  private:

    /// comunication (mainly for screen output)
    Epetra_Comm& comm_;

    /// method name prepared for output
    std::string method_;

    //! @name Time stepping variables
    int step_;
    int nstep_;
    double time_;
    double maxtime_;
    double dt_;
    //@}
  };
}

#endif
