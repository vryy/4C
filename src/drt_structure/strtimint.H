/*----------------------------------------------------------------------*/
/*!
\file strtimint.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifdef CCADISCRET
#ifndef STRTIMINT_H
#define STRTIMINT_H

/*----------------------------------------------------------------------*/
/* headers */
#include <iostream>
#include <fstream>

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include <NOX.H>
#include <NOX_Epetra.H>
#include "ml_common.h"
#include "ml_include.h"
#include "nlnml_preconditioner.H"
#include "nlnml_linearsystem.H"
#include "nlnml_finelevelnoxinterface.H"

#include "../drt_io/io_control.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include "../drt_io/io.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_constraint/constraintsolver.H"
#include "../drt_surfstress/drt_surfstress_manager.H"
#include "../drt_surfstress/drt_potential_manager.H"
#include "../drt_contact/drt_contact_manager.H"



#include "strtimint_vector.H"
#include "strtimint_mstep.H"

/*----------------------------------------------------------------------*/
//! Structural Dynamics
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Front-end for structural dynamics by integrating in time.
   *
   * <h3> Intention </h3>
   * This front-end for structural dynamics defines an interface to call
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
   * the tangent system matrix #stiff_ by #Stiff(). This selective access
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

    //! Structural time adaptivity is friend
    friend class TimAda;

    //! Joint auxiliar schemes are friends
    template <typename T> friend class TimAdaJoint;

  public:

    //! Type of damping
    enum DampEnum
    {
      damp_none,  //!< damping off
      damp_rayleigh,  //!< globally applied Rayleigh damping
      damp_material  //!< element-wise applied damping
    };

    //! Map inpu string to damping enum
    static enum DampEnum MapDampStringToEnum
    (
      const std::string name  //!< identifier
    );

    //! Type of stress
    enum StressEnum
    {
      stress_none,  //!< undetermined
      stress_cauchy,  //!< Cauchy stress
      stress_pk1,  //!< 1st Piola-Kirchhoff stress
      stress_pk2  //! 2nd Piola-Kirchhoff stress
    };

    //! Map stress input string to enum
    static enum StressEnum MapStressStringToEnum
    (
      const std::string name  //!< identifier
    );

    //! Type of strain
    enum StrainEnum
    {
      strain_none,  //!< undetermined
      strain_ea,  //!< Euler-Almansi strain
      strain_gl  //!< Green-Lagrange strain
    };

    //! Map strain input string to enum
    static enum StrainEnum MapStrainStringToEnum
    (
      const std::string name  //!< identifier
    );

    //! List of all available time integration schemes
    enum NameEnum
    {
      name_statics,  //!< static analysis
      name_genalpha,  //!< generalised-alpha
      name_onesteptheta,  //!< one-step-theta
      name_gemm,  //!< generalised energy-momentum method
      name_ab2  //!< Adams-Bashforth 2nd order
    };

    //! Map name enum term to string
    static std::string MapNameEnumToString
    (
      const enum NameEnum term  //!< the enum
    );

    //! @name Life
    //@{

    //! Print tea time logo
    void Logo();

    //! Constructor
    TimInt
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimInt()
    {
      DetachEnergyFile();
    }

    //! Empty constructor
    TimInt() { ; }

    //! Copy constructor
    TimInt(const TimInt& old) { ; }

    //! Resize #TimIntMStep<T> multi-step quantities
    virtual void ResizeMStep() = 0;

    //! Merge
    //!
    //! Merge basically duplicates the base object content of time
    //! integrator #tis onto the time integrator #this. This is like
    //! a copy, but a copy constructor is not permitted, because
    //! #TimInt is pure virtual.
    //! Usually this is not wanted when copying, but here it is
    //! highly appreciated. #TimInt contains only pointers (of the
    //! RefCount type) and can thus link -- or merge -- the data
    //! of #tis with #this. Practically, this turns up with time
    //! adaptivity in which #tis is the marching integrator
    //! and #this is the auxiliar method, which shares the marching data.
    void Merge
    (
      TimInt& tis  //!< existing integrator to merge against
    )
    {
      // copy it
      *this = tis;
      // A slice of lemon?
      return;
    }

    //@}

    //! @name Actions
    //@{

    //! Equilibrate the initial state by identifying the consistent
    //! initial accelerations and (if applicable) internal variables
    //! Make damping and mass matrix
    void DetermineMassDampConsistAccel();

    //! Apply Dirichlet boundary conditions on provided state vectors
    void ApplyDirichletBC
    (
      const double time,  //!< at time
      Teuchos::RCP<Epetra_Vector> dis,  //!< displacements
                                        //!< (may be Teuchos::null)
      Teuchos::RCP<Epetra_Vector> vel,  //!< velocities
                                        //!< (may be Teuchos::null)
      Teuchos::RCP<Epetra_Vector> acc,  //!< accelerations
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
    virtual void UpdateStep() = 0;

    //! Update configuration and time after time step
    virtual void UpdateStepAndTime() = 0;

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
    //! \f$D_{n}, V_{n}, A_{n}\f$. So, #UpdateIncrement should be called
    //! upon object prior to writing stuff here.
    //! \author mwgee (originally) \date 03/07
    void OutputStep();

    //! Write restart
    //! \author mwgee (originally) \date 03/07
    void OutputRestart
    (
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );

    //! Output displacements, velocities and accelerations
    //! and more system vectors
    //! \author mwgee (originally) \date 03/07
    void OutputState
    (
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );

    //! Stress & strain output
    //! \author lw (originally)
    void OutputStressStrain
    (
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

    //! @name Forces and Stiffnesses
    //@{

    //! Apply external force
    void ApplyForceExternal
    (
      const double time,  //!< evaluation time
      const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
      const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
      Teuchos::RCP<Epetra_Vector>& fext  //!< external force
    );

    //! Evaluate ordinary internal force, its stiffness at state
    void ApplyForceStiffInternal
    (
      const double time,  //!< evaluation time
      const double dt,  //!< step size
      const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
      const Teuchos::RCP<Epetra_Vector> disi,  //!< residual displacements
      const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
      Teuchos::RCP<Epetra_Vector> fint,  //!< internal force
      Teuchos::RCP<LINALG::SparseMatrix> stiff  //!< stiffness matrix
    );

    //! Evaluate ordinary internal force
    //!
    //! We need incremental displacements, because the internal variables,
    //! chiefly EAS parameters with an algebraic constraint, are treated
    //! as well. They are not treated perfectly, ie they are not iteratively
    //! equilibriated according to their (non-linear) constraint and
    //! the pre-determined displacements -- we talk explicit time integration
    //! here, but they are applied in linearised manner. The linearised
    //! manner means the static condensation is applied once with
    //! residual displacements replaced by the full-step displacement
    //! increment \f$D_{n+1}-D_{n}\f$.
    void ApplyForceInternal
    (
      const double time,  //!< evaluation time
      const double dt,  //!< step size
      const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
      const Teuchos::RCP<Epetra_Vector> disi,  //!< incremental displacements
      const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
      Teuchos::RCP<Epetra_Vector> fint  //!< internal force
    );

    //@}

    //! @name Attributes
    //@{

    //! Provide Name
    virtual enum NameEnum MethodName() = 0;

    //! Provide title
    std::string MethodTitle()
    {
      return MapNameEnumToString(MethodName());
    }

    //! Return true, if time integrator is implicit
    virtual bool MethodImplicit() = 0;

    //! Return true, if time integrator is explicit
    bool MethodExplicit() { return (not MethodImplicit()); }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a \f$m\f$-multistep method returns \f$m\f$
    virtual int MethodSteps() = 0;

    //! Give order of accuracy
    int MethodOrderOfAccuracy()
    {
      return min(MethodOrderOfAccuracyDis(), MethodOrderOfAccuracyVel());
    }

    //! Give local order of accuracy of displacement part
    virtual int MethodOrderOfAccuracyDis() = 0;

    //! Give local order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() = 0;

    //! Return linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis() = 0;

    //! Return linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel() = 0;

    //@}

    //! @name Access methods
    //@{

    //! Access dicretisation
    Teuchos::RCP<DRT::Discretization> Discretization()
    {
      // Here a 'false' must be used. This is due to
      // the fact that TimInt possesses a references
      // on the discretisation #discret_ and not
      // a Teuchos::RefCountPointer. Eventually, TimInt
      // will be destroyed and it will immediately destroy
      // its #discret_ member. However, #discret_ is handed down
      // to the ConstrManager and kept there as a RefCountPointer.
      // The object #discret_ is gone, when ConstrManager tries
      // to kill it. We achieve a nice segmentation fault.
      // The 'false' prevents ConstrManager of trying to kill it.
      //return Teuchos::rcp(&discret_, false);

      // Now, the discretisation is stored as Teuchos::RefCountPointer,
      // thus
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
    void ReadRestart
    (
      const int step //!< restart step
    );

    //! Read and set restart state
    void ReadRestartState();

    //! Read and set restart forces
    virtual void ReadRestartForce() = 0;

    //! Read and set restart values for constraints
    void ReadRestartConstraint();

    //! Read and set restart values for surface stresses
    void ReadRestartSurfstress();

    //! Read and set restart values for multi scale materials
    void ReadRestartMultiScale();

    //! Return displacements \f$D_{n}\f$
    Teuchos::RCP<Epetra_Vector> Dis() { return (*dis_)(0); }

    //! Return displacements \f$D_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> DisNew() { return disn_; }

    //! Return velocities \f$V_{n}\f$
    Teuchos::RCP<Epetra_Vector> Vel() { return (*vel_)(0); }

    //! Return velocities \f$V_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> VelNew() { return veln_; }

    //! Return accelerations \f$A_{n}\f$
    Teuchos::RCP<Epetra_Vector> Acc() { return (*acc_)(0); }

    //! Return accelerations \f$A_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> AccNew() { return accn_; }

    //! Return external force \f$F_{ext,n}\f$
    virtual Teuchos::RCP<Epetra_Vector> Fext() = 0;

    //! Return stiffness,
    //! i.e. force residual differentiated by displacements
    Teuchos::RCP<LINALG::SparseMatrix> Stiff() { return stiff_; }

    //! Return domain map
    const Epetra_Map& GetDomainMap() { return mass_->DomainMap(); }

    //! Return current step number map
    int GetStep() { return step_; }

    //! Return current step number map
    double GetTime() { return (*time_)[0]; }


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
    bool writestate_;  //!< write state on/off
    int writestateevery_;  //!< write state every given step
    int writestrevery_;  //!< write stress/strain every given step
    enum StressEnum writestress_;  //!< stress output type
    enum StrainEnum writestrain_;  //!< strain output type
    int writeenergyevery_;  //!< write system energy every given step
    std::ofstream* energyfile_;  //!< outputfile for energy
    //@}

    //! @name Damping
    //!
    //! Rayleigh damping means \f${C} = c_\text{K} {K} + c_\text{M} {M}\f$
    //@{
    enum DampEnum damping_;  //!< damping type
    double dampk_;  //!< damping factor for stiffness \f$c_\text{K}\f$
    double dampm_;  //!< damping factor for mass \f$c_\text{M}\f$
    //@}

    //! @name Managed stuff
    //@{

    //! whatever constraints
    Teuchos::RCP<UTILS::ConstrManager> conman_;  //!< constraint manager
    Teuchos::RCP<UTILS::ConstraintSolver> consolv_;  //!< constraint solver

    //! surface stresses
    Teuchos::RCP<UTILS::SurfStressManager> surfstressman_;

    //! potential conditions
    Teuchos::RCP<UTILS::PotentialManager> potman_;

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

    //! global displacements \f${D}_{n}, D_{n-1}, ...\f$
    Teuchos::RCP<TimIntMStep<Epetra_Vector> > dis_;
    //! global velocities \f${V}_{n}, V_{n-1}, ...\f$
    Teuchos::RCP<TimIntMStep<Epetra_Vector> > vel_;
    //! global accelerations \f${A}_{n}, A_{n-1}, ...\f$
    Teuchos::RCP<TimIntMStep<Epetra_Vector> > acc_;
    Teuchos::RCP<Epetra_Vector> disn_;  //!< global displacements
                                        //!< \f${D}_{n+1}\f$
                                        //!< at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> veln_;  //!< global velocities
                                        //!< \f${V}_{n+1}\f$
                                        //!< at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> accn_;  //!< global accelerations
                                        //!< \f${A}_{n+1}\f$
                                        //!< at \f$t_{n+1}\f$

    //@}

    //! @name System matrices
    //@{
    //! holds eventually effective stiffness
    Teuchos::RCP<LINALG::SparseMatrix> stiff_;
    //! mass matrix (constant)
    Teuchos::RCP<LINALG::SparseMatrix> mass_;
    //! Rayleigh damping matrix (constant)
    Teuchos::RCP<LINALG::SparseMatrix> damp_;
    //@}

  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_H
#endif  // #ifdef CCADISCRET
