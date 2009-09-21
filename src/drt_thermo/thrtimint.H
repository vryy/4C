/*----------------------------------------------------------------------*/
/*!
\file thrtimint.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                             bborn 06/08 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef THRTIMINT_H
#define THRTIMINT_H

/*----------------------------------------------------------------------*
 |  headers                                                 bborn 06/08 |
 *----------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <fstream>

#include "Teuchos_RefCountPtr.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include "../drt_io/io_control.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_thermo.H"
#include "../drt_io/io.H"

#include "thrtimint_mstep.H"

/*----------------------------------------------------------------------*/
//! Thermal Dynamics
namespace THR
{

  /*====================================================================*/
  /*!
   * \brief Front-end for thermal dynamics by integrating in time.
   *
   * <h3> Intention </h3>
   * This front-end for thermal dynamics defines an interface to call
   * several derived time integrators. Thus it describes a plethora of pure
   * virtual methods which have to be implemented at the derived integrators.
   * However, it also offers a few non-empty methods and stores associated
   * data. The most important method of this base time integrator object
   * is #Integrate().
   *
   * #Integrate() performs a time integration (time loop) with constant
   * time steps and other parameters as set by the user.
   *
   * Although #Integrate is the main interface, this base time integrator
   * allows the public to access a few of its datum objects, for instance
   * the tangent system matrix #tang_ by #Tang(). This selective access
   * is needed in environments in which a independent time loop is provided.
   * This happens e.g. in fluid-structure-interaction.
   *
   * <h3> Responsibilties </h3>
   * Most importantly the base integrator manages the system state vectors and
   * matrices. It also deals with the output to files and offers method to
   * determine forces and stiffnesses (tangents).
   *
   * \author bborn
   * \date 06/08
   */
  class TimInt
  {

  public:

    //! @name Life
    //@{

    //! Print tea time logo
    void Logo();

    //! Constructor
    TimInt
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& tdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Empty constructor
    TimInt() { ; }

    //! Copy constructor
    TimInt(const TimInt& old) { ; }

    //! Destructor
    virtual ~TimInt()
    {
      DetachEnergyFile();
    }

    //! Resize #TimIntMStep<T> multi-step quantities
    virtual void ResizeMStep() = 0;

    //@}

    //! @name Actions
    //@{

    //! Equilibrate the initial state by identifying the consistent
    //! initial accelerations and (if applicable) internal variables
    //! Make capacity matrix
    void DetermineCapaConsistTempRate();

    //! Apply Dirichlet boundary conditions on provided state vectors
    void ApplyDirichletBC(
      const double time,  //!< at time
      Teuchos::RCP<Epetra_Vector> temp,  //!< temperatures
                                        //!< (may be Teuchos::null)
      Teuchos::RCP<Epetra_Vector> rate,  //!< temperature rate
                                        //!< (may be Teuchos::null)
      bool recreatemap  //!< recreate mapextractor/toggle-vector
                        //!< which stores the DOF IDs subjected
                        //!< to Dirichlet BCs
                        //!< This needs to be true if the bounded DOFs
                        //!< have been changed.
    );

    //! Do time integration of multiple steps
    void Integrate();

    //! Do time integration of single step
    virtual void IntegrateStep() = 0;

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awating the next
    //! time step.
    virtual void UpdateStepState() = 0;

    //! Update time and step counter
    void UpdateStepTime();

    //! Reset configuration after time step
    //!
    //! Thus the last converged state is copied back on the predictor
    //! for current time step. This applies only to elemet-wise
    //! quantities
    void ResetStep();

    //@}

    //! @name Output
    //@{

    //! print summary after step
    virtual void PrintStep() = 0;

    //! Output to file
    //! This routine always prints the last converged state, i.e.
    //! \f$T_{n}, R_{n}\f$. So, #UpdateIncrement should be called
    //! upon object prior to writing stuff here.
    //! \author mwgee (originally) \date 03/07
    void OutputStep();

    //! Write restart
    //! \author mwgee (originally) \date 03/07
    void OutputRestart(
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );

    //! Output temperatures, temperature rates
    //! and more system vectors
    //! \author mwgee (originally) \date 03/07
    void OutputState(
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );

    //! Heatflux & temperature gradient output
    //! \author lw (originally)
    void OutputHeatfluxTempgrad(
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );

    //! Energy output
    void OutputEnergy();

    //! Check wether energy output file is attached
    bool AttachedEnergyFile()
    {
      if (energyfile_) return true;
      else return false;
    }

    //! Attach file handle for energy file #energyfile_
    void AttachEnergyFile()
    {
      if (not energyfile_)
      {
        std::string energyname
          = DRT::Problem::Instance()->OutputControlFile()->FileName()
          + ".energy";
        energyfile_ = new std::ofstream(energyname.c_str());
        *energyfile_ << "# timestep time total_energy"
                     << " kinetic_energy internal_energy external_energy"
                     << std::endl;
      }
    }

    //! Detach file handle for energy file #energyfile_
    void DetachEnergyFile()
    {
      if (energyfile_) delete energyfile_;
    }

    //@}

    //! @name Forces and Tangent
    //@{

    //! Apply external force
    void ApplyForceExternal(
      const double time,  //!< evaluation time
      const Teuchos::RCP<Epetra_Vector> temp,  //!< temperature state
      Teuchos::RCP<Epetra_Vector>& fext  //!< external force
    );

    //! Evaluate ordinary internal force, its stiffness at state
    void ApplyForceTangInternal(
      Teuchos::ParameterList& p,  //!< parameter list handed down to elements
      const double time,  //!< evaluation time
      const double dt,  //!< step size
      const Teuchos::RCP<Epetra_Vector> temp,  //!< temperature state
      const Teuchos::RCP<Epetra_Vector> tempi,  //!< residual temperatures
      Teuchos::RCP<Epetra_Vector> fint,  //!< internal force
      Teuchos::RCP<LINALG::SparseMatrix> tang  //!< stiffness matrix
    );

    //! Evaluate ordinary internal force
    //!
    //! We need incremental temperatures, because the internal variables,
    //! chiefly EAS parameters with an algebraic constraint, are treated
    //! as well. They are not treated perfectly, ie they are not iteratively
    //! equilibriated according to their (non-linear) constraint and
    //! the pre-determined temperatures -- we talk explicit time integration
    //! here, but they are applied in linearised manner. The linearised
    //! manner means the static condensation is applied once with
    //! residual temperatures replaced by the full-step temperature
    //! increment \f$T_{n+1}-T_{n}\f$.
    void ApplyForceInternal
    (
      Teuchos::ParameterList& p,  //!< parameter list handed down to elements
      const double time,  //!< evaluation time
      const double dt,  //!< step size
      const Teuchos::RCP<Epetra_Vector> temp,  //!< temperature state
      const Teuchos::RCP<Epetra_Vector> tempi,  //!< incremental temperatures
      Teuchos::RCP<Epetra_Vector> fint  //!< internal force
    );

    //@}

    //! @name Attributes
    //@{

    //! Provide Name
    virtual enum INPAR::THR::DynamicType MethodName() const = 0;

    //! Provide title
    std::string MethodTitle() const
    {
      return INPAR::THR::DynamicTypeString(MethodName());
    }

    //! Return true, if time integrator is implicit
    virtual bool MethodImplicit() = 0;

    //! Return true, if time integrator is explicit
    bool MethodExplicit() { return (not MethodImplicit()); }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a \f$m\f$-multistep method returns \f$m\f$
    virtual int MethodSteps() = 0;

    //! Give order of accuracy
    virtual int MethodOrderOfAccuracy() = 0;

    //! Return linear error coefficient of temperatures
    virtual double MethodLinErrCoeff() = 0;

    //@}

    //! @name Access methods
    //@{

    //! Access dicretisation
    Teuchos::RCP<DRT::Discretization> Discretization()
    {
      return discret_;
    }

    //! Access solver
    Teuchos::RCP<LINALG::Solver> Solver()
    {
      return solver_;
    }

    //! Access output object
    Teuchos::RCP<IO::DiscretizationWriter> DiscWriter()
    {
      return output_;
    }

    //! Read restart values
    void ReadRestart(
      const int step //!< restart step
    );

    //! Read and set restart state
    void ReadRestartState();

    //! Read and set restart forces
    virtual void ReadRestartForce() = 0;

    //! Return temperatures \f$T_{n}\f$
    Teuchos::RCP<Epetra_Vector> Temp() { return (*temp_)(0); }

    //! Return temperatures \f$T_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> TempNew() { return tempn_; }

    //! Return temperature rates \f$V_{n}\f$
    Teuchos::RCP<Epetra_Vector> Rate() { return (*rate_)(0); }

    //! Return temperature rates \f$V_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> RateNew() { return raten_; }

    //! Return external force \f$F_{ext,n}\f$
    virtual Teuchos::RCP<Epetra_Vector> Fext() = 0;

    //! Return reaction forces
    virtual Teuchos::RCP<Epetra_Vector> Freact() = 0;

    //! Return tangent,
    //! i.e. force residual differentiated by temperatures
    Teuchos::RCP<LINALG::SparseMatrix> Tang() { return tang_; }

    //! Return domain map
    const Epetra_Map& GetDomainMap() { return tang_->DomainMap(); }

    //! Return current step number map
    double GetTime() const { return (*time_)[0]; }

    //! Get upper limit of time range of interest
    double GetTimeEnd() const { return timemax_; }

    //! Get time step size \f$\Delta t_n\f$
    double GetTimeStepSize() const { return (*dt_)[0]; }

    //! Return current step number map
    int GetStep() const { return step_; }

    //! Get number of time steps
    int GetTimeNumStep() const { return stepmax_; }


    //! Return MapExtractor for Dirichlet boundary conditions
    const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() const
    {
      return dbcmaps_;
    }

    //@}

  protected:

    //! @name General purpose algorithm members
    //@{
    Teuchos::RCP<DRT::Discretization> discret_;  //!< attached discretisation
    int myrank_;  //!< ID of actual processor in parallel
    const Epetra_Map* dofrowmap_;  //!< vector layout from the discretization
                                   //!< to construct matching
                                   //!< vectors and matrices
    Teuchos::RCP<LINALG::Solver> solver_;  //!< linear algebraic solver
    bool solveradapttol_;  //!< adapt solver tolerance
    double solveradaptolbetter_;  //!< tolerance to which is adpated ????
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;  //!< map extractor object
                                                  //!< containing non-overlapping
                                                  //!< map of global DOFs on Dirichlet
                                                  //!< boundary conditions
    //@}

    //! @name Printing and output
    //@{
    Teuchos::RCP<IO::DiscretizationWriter> output_;  //!< binary output
    bool printlogo_;  //!< true: enjoy your cuppa
    bool printscreen_;  //!< print infos to standard out
    FILE* errfile_;  //!< error file handle
    bool printerrfile_;  //!< print infos to error file #errfile_
    bool printiter_;  //!< print intermediate iterations during solution
    int writerestartevery_;  //!< write restart every given step;
                             //!< if 0, restart is not written
    bool writeglob_;  //!< write state on/off
    int writeglobevery_;  //!< write state every given step
    int writeelemevery_;  //!< write stress/strain every given step
    INPAR::THR::HeatFluxType writeheatflux_;
    INPAR::THR::TempGradType writetempgrad_;
    int writeenergyevery_;  //!< write system energy every given step
    std::ofstream* energyfile_;  //!< outputfile for energy
    //@}

    //! @name General control parameters
    //@{
    Teuchos::RCP<TimIntMStep<double> > time_;  //!< time \f$t_{n}\f$ of last converged step
    double timen_;  //!< target time \f$t_{n+1}\f$
    Teuchos::RCP<TimIntMStep<double> > dt_;  //!< time step size \f$\Delta t\f$
    double timemax_;  //!< final time \f$t_\text{fin}\f$
    int stepmax_;  //!< final step \f$N\f$
    int step_;  //!< time step index \f$n\f$
    int stepn_;  //!< time step index \f$n+1\f$
    //@}

    //! @name Global vectors
    //@{
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
    //@}


    //! @name Global state vectors
    //@{

    //! global temperatures \f${T}_{n}, T_{n-1}, ...\f$
    Teuchos::RCP<TimIntMStep<Epetra_Vector> > temp_;
    //! global temperature rates \f${R}_{n}, R_{n-1}, ...\f$
    Teuchos::RCP<TimIntMStep<Epetra_Vector> > rate_;
    Teuchos::RCP<Epetra_Vector> tempn_;  //!< global temperatures
                                        //!< \f${T}_{n+1}\f$
                                        //!< at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> raten_;  //!< global temperature rates
                                        //!< \f${R}_{n+1}\f$
                                        //!< at \f$t_{n+1}\f$

    //@}

    //! @name System matrices
    //@{
    //! holds eventually effective tangent
    Teuchos::RCP<LINALG::SparseMatrix> tang_;
    //! capacity matrix (constant)
    Teuchos::RCP<LINALG::SparseMatrix> capa_;

    //@}

  };  // TimInt

}  // namespace THR


#endif  // #ifndef THRTIMINT_H
#endif  // #ifdef CCADISCRET
