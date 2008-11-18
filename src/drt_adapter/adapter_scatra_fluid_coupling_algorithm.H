/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_fluid_coupling_algorithm.H

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and (active or passive) scalar transport equations

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_SCATRA_FLUID_COUPLING_ALGORITHM_H
#define ADAPTER_SCATRA_FLUID_COUPLING_ALGORITHM_H

#include "../drt_adapter/adapter_fluid_base_algorithm.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

namespace ADAPTER
{

  /// basis coupling algorithm for scalar transport with fluid velocity field
  /*!

    Base class for scalar transport problems coupled to Navier-Stokes velocity field. 
    Derives from FluidBaseAlgorithm and ScaTraBaseAlgorithm.
    There are different subclasses that implement different coupling schemes 
    (one-way coupling for the transport of passive scalars and 
     fully coupled schemes for other multiphysics applications like combustion or electrochemistry).

    \author gjb
    \date 07/08
   */
  class ScaTraFluidCouplingAlgorithm : 
        public ADAPTER::FluidBaseAlgorithm,
        public ADAPTER::ScaTraBaseAlgorithm
  {
  public:

    /// constructor using a Epetra_Comm
    explicit ScaTraFluidCouplingAlgorithm(
        Epetra_Comm& comm,              ///< communicator
        const Teuchos::ParameterList& prbdyn  ///< problem-specific parameters
        );

    /// virtual destructor to support polymorph destruction
    virtual ~ScaTraFluidCouplingAlgorithm();

    /// outer level time loop (to be implemented by deriving classes)
    virtual void TimeLoop() = 0;

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

    /// read restart data
    void ReadRestart(
        int step ///< step number where the calculation is continued
        );

  protected:

    /// time step size
    double Dt() const { return dt_;}

    /// time step number
    int Step() const { return step_;}

    /// total number of time steps
    int NStep() const { return nstep_;}

    /// maximum simulation time
    double MaxTime() const { return maxtime_;}

    /// current time
    double Time() const {return time_;}

    /// tests whether there are more time steps to do or not
    bool NotFinished() { return ((step_ < nstep_) and ((time_ + EPS12) < maxtime_)); }

    /// get current velocity from fluid solver
    void GetCurrentFluidVelocity(){ velocitynp_=FluidField().ExtractVelocityPart(FluidField().Velnp()); }

    /// current velocity from fluid solver
    RCP<const Epetra_Vector> ConvectiveVelocity(){ return velocitynp_; }

    /// get current velocity+pressure from fluid solver
    void GetCurrentFluidVelPress(){ velpressnp_=FluidField().Velnp(); }

    /// current velocity from fluid solver
    RCP<const Epetra_Vector> VelocityPressure(){ return velpressnp_; }

    /// get current density from scalar transport solver
    void GetCurrentDensity(){ densnp_=ScaTraField().Densnp(); }

    /// current density from scalar transport solver
    RCP<const Epetra_Vector> Density(){ return densnp_; }

    /// get density at time step n from scalar transport solver
    void GetNDensity(){ densn_=ScaTraField().Densn(); }

    /// density at time step n from scalar transport solver
    RCP<const Epetra_Vector> NDensity(){ return densn_; }

    /// get density at time step n-1 from scalar transport solver
    void GetNmDensity(){ densnm_=ScaTraField().Densnm(); }

    /// density at time step n-1 from scalar transport solver
    RCP<const Epetra_Vector> NmDensity(){ return densnm_; }

    /// increment time and step value
    void IncrementTimeAndStep()
    {
      step_ += 1;
      time_ += dt_;
    }

  private:

    /// communication (mainly for screen output)
    Epetra_Comm& comm_;

    /// current velocity field at time n+1
    RCP<const Epetra_Vector>    velocitynp_;

    /// current velocity and pressure field at time n+1
    RCP<const Epetra_Vector>    velpressnp_;

    /// density fields at times n+1, n and n-1
    RCP<const Epetra_Vector>    densnp_;
    RCP<const Epetra_Vector>    densn_;
    RCP<const Epetra_Vector>    densnm_;

    /// problem-specific parameter list
    const Teuchos::ParameterList& params_;

    //! @name Time stepping variables
    /// time step counter
    int step_;
    /// maximum number of time steps
    int nstep_;
    /// current simulation time
    double time_;
    /// maximum simulation time
    double maxtime_;
    /// time step size
    double dt_;
    //@}

  };

} // namespace ADAPTER

#endif // ADAPTER_SCATRA_FLUID_COUPLING_ALGORITHM_H
#endif // CCADISCRET
