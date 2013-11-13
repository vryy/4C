/*----------------------------------------------------------------------*/
/*!
\file strtimint.H
\brief Time integration for structural dynamics

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef STRTIMINT_H
#define STRTIMINT_H

/*----------------------------------------------------------------------*/
/* headers */
#include <string>
#include <iostream>
#include <fstream>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Time.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <ml_common.h>
#include <ml_include.h>
#include <nlnml_preconditioner.H>
#include <nlnml_linearsystem.H>
#include <nlnml_finelevelnoxinterface.H>

#include "../drt_adapter/ad_str_structure.H"
#include "strtimint_mstep.H"
#include "../drt_inpar/inpar_structure.H"


// forward declarations
namespace DRT
{
  class Discretization;

  namespace UTILS
  {
    class LocsysManager;
  }
}

namespace UTILS
{
  class ConstrManager;
  class ConstraintSolver;
  class SurfStressManager;
//  class GrowthDisplManager;
}

namespace MORTAR
{
  class ManagerBase;
}

namespace CONTACT
{
  class Beam3cmanager;
}

namespace STATMECH
{
  class StatMechManager;
}

namespace UTILS
{
  class PlastSsnManager;
}

namespace POTENTIAL
{
  class PotentialManager;
}

namespace LINALG
{
  class Solver;
  class MapExtractor;
  class SparseMatrix;
  class SparseOperator;
  class BlockSparseMatrixBase;
}

namespace IO
{
  class DiscretizationWriter;
}


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
   * the tangent system matrix #stiff_ by #SystemMatrix(). This selective access
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
  class TimInt : public ::ADAPTER::Structure
  {

    //! Structural time adaptivity is friend
    friend class TimAda;

    //! Joint auxiliar schemes are friends
    template <typename T> friend class TimAdaJoint;

  public:

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
      Teuchos::RCP<LINALG::Solver> contactsolver,  //!< the solver for contact/meshtying
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

    //! Resize #TimIntMStep<T> multi-step quantities, needed for fsi time adaptivity
    virtual void ResizeMStepTimAda();

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

    /// start new time step
    virtual void PrepareTimeStep() = 0;

    //! Do time integration of multiple steps
    int Integrate();

    //! do something in case nonlinear soluttion does not converge for some reason
    int PerformErrorAction(int nonlinsoldiv);

    //! Do time integration of single step
    virtual int IntegrateStep() = 0;


    /// non-linear solve
    ///
    /// Do the nonlinear solve for the time step. All boundary conditions have
    /// been set.
    virtual int Solve() = 0;

    //! Linear structure solve with just an interface load
    virtual Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear() = 0;

    /// update displacement and evaluate elements
    ///
    /// there are two displacement increments possible
    ///
    /// x^n+1_i+1 = x^n+1_i + disiterinc  (sometimes referred to as residual increment), and
    ///
    /// x^n+1_i+1 = x^n     + disstepinc
    ///
    /// with n and i being time and Newton iteration step
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) = 0;

    /// update at time step end
    virtual void Update() = 0;

    /// update at time step end in case of FSI time adaptivity
    virtual void Update(const double endtime) = 0;

    /// Update iteration
    /// Add residual increment to Lagrange multipliers stored in Constraint manager
    virtual void UpdateIterIncrConstr
    (
      Teuchos::RCP<Epetra_Vector> lagrincr  ///< Lagrange multiplier increment
    ) = 0;

    /// output results
    virtual void Output(bool forced_writerestart = false) = 0;

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awating the next
    //! time step.
    virtual void UpdateStepState() = 0;

    //! Update everything on element level after time step and after output
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awating the next
    //! time step.
    virtual void UpdateStepElement() = 0;

    //! Update time and step counter
    void UpdateStepTime();

    //! Update step for contact / meshtying
    void UpdateStepContactMeshtying();

    //! Velocity update method (VUM) for contact
    //!
    //! The VUM is an explicit update method at the end of each time step
    //! which is supposed to assure exact algorithmic conservation of total
    //! energy during contact. Further details can be found in the original
    //! paper by Laursen and Love (IJNME, 2002) and in the more recent and
    //! mortar-related paper by Hartmann et al. (IJNME, 2007).
    //! CAUTION: The VUM is only available for Gen-Alpha and GEMM.
    void UpdateStepContactVUM();

    //! Update step for beam contact
    void UpdateStepBeamContact();

    //! Reset configuration after time step
    //!
    //! Thus the last converged state is copied back on the predictor
    //! for current time step. This applies only to elemet-wise
    //! quantities
    void ResetStep();

    //! Set initial fields in structure (e.g. initial velocities)
    void SetInitialFields();

    //@}

    //! @name Determination of output quantities that depend on
    // the constitutive model
    //@{
    //! Calculate all output quantities depending on the constitutive model
    //  (and, hence, on a potential material history)
    void PrepareOutput();

    //! Calculate stresses, strains on micro-scale
    void PrepareOutputMicro();

    //! Calculate stresses and strains
    void DetermineStressStrain();

    //! Calculate kinetic, internal and external energy
    virtual void DetermineEnergy();

    /// create result test for encapsulated structure algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

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
    void OutputStep(bool forced_writerestart = false);

    //! Write restart
    //! \author mwgee (originally) \date 03/07
    virtual void OutputRestart
    (
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );
    //! Get data that is written during restart
    //! \author biehler \data 06/13
    virtual void GetRestartData
    (
      Teuchos::RCP<int> step,
      Teuchos::RCP<double> time,
      Teuchos::RCP<Epetra_Vector> disn, //!< new displacement state
      Teuchos::RCP<Epetra_Vector> veln, //!< new velocity state
      Teuchos::RCP<Epetra_Vector> accn, //!< new acceleration state
      Teuchos::RCP<std::vector<char> > elementdata //!< internal element/history variables e.g. F_prestress
    );

    //! Output displacements, velocities and accelerations
    //! and more system vectors
    //! \author mwgee (originally) \date 03/07
    virtual void OutputState
    (
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );

    //! Add restart information to OutputState
    void AddRestartToOutputState();

    //! Stress & strain output
    //! \author lw (originally)
    void OutputStressStrain
    (
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );

    //! Energy output
    void OutputEnergy();

    //! Active set, energy and momentum output for contact
    void OutputContact();

    //! Error norm output
    void OutputErrorNorms();

    //! Nodal positions output
    void OutputNodalPositions();

    //! Output on the micro-scale (multi-scale analysis)
    void OutputMicro();

    //! Patient Specific Output
    void OutputPatspec();

    //! Write internal and external forces (if necessary for restart)
    virtual void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output) = 0;

    //! Check whether energy output file is attached
    bool AttachedEnergyFile()
    {
      if (energyfile_) return true;
      else return false;
    }

    //! Attach file handle for energy file #energyfile_
    virtual void AttachEnergyFile(std::string name = "");

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
      const Teuchos::RCP<Epetra_Vector> dis,  //!< old displacement state
      const Teuchos::RCP<Epetra_Vector> disn,  //!< new displacement state
      const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
      Teuchos::RCP<Epetra_Vector>& fext,  //!< external force
      Teuchos::RCP<LINALG::SparseOperator>& fextlin //!<linearization of external force
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
      Teuchos::RCP<LINALG::SparseOperator> stiff,  //!< stiffness matrix
      Teuchos::RCP<LINALG::SparseOperator> damp = Teuchos::null  //!< material damping matrix
    );

    //! Evaluate internal and inertial forces and their linearizations
    void ApplyForceStiffInternalAndInertial
    (
      const double time,  //!< evaluation time
      const double dt,  //!< step size
      const double timintfac_dis,  //!< time integration factor for d(Res) / d dis
      const double timintfac_vel,  //!< time integration factor for d(Res) / d vel
      const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
      const Teuchos::RCP<Epetra_Vector> disi,  //!< residual displacements
      const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
      const Teuchos::RCP<Epetra_Vector> acc,  // acceleration state
      Teuchos::RCP<Epetra_Vector> fint,  //!< internal force
      Teuchos::RCP<Epetra_Vector> finert,  //!< inertial force
      Teuchos::RCP<LINALG::SparseOperator> stiff,  //!< stiffness matrix
      Teuchos::RCP<LINALG::SparseOperator> mass  //!< mass matrix
    );


    //! Return bool indicating if we have nonlinear inertia forces
    bool HaveNonlinearMass();

    //! check wether the initial conditions are fulfilled */
    void NonlinearMassSanityCheck( Teuchos::RCP<Epetra_Vector> fext,
                                   Teuchos::RCP<Epetra_Vector> dis,
                                   Teuchos::RCP<Epetra_Vector> vel,
                                   Teuchos::RCP<Epetra_Vector> acc
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

    //! @name porous media specific methods
    //@{

    //@}

    //! Set forces due to interface with fluid, the force is expected external-force-like
    virtual void SetForceInterface
    (
      Teuchos::RCP<Epetra_MultiVector> iforce  ///< the force on interface
    );

    //! @name Attributes
    //@{

    //! Provide Name
    virtual enum INPAR::STR::DynamicType MethodName() const = 0;

    //! Provide title
    std::string MethodTitle() const
    {
      return INPAR::STR::DynamicTypeString(MethodName());
    }

    //! Return true, if time integrator is implicit
    virtual bool MethodImplicit() = 0;

    //! Return true, if time integrator is explicit
    bool MethodExplicit() { return (not MethodImplicit()); }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a \f$m\f$-multistep method returns \f$m\f$
    virtual int MethodSteps() = 0;

    //! return time integration factor
    virtual const double TimIntParam() const = 0;

    //! return Newmark parameter \f$beta\f$
    virtual const double GetNewmarkBeta() const = 0;

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

    //! Access to rcp pointer to dofrowmap of discretization
    const Epetra_Map* DofRowMapView();

    //! Access solver, one of these have to be removed (see below)
    Teuchos::RCP<LINALG::Solver> Solver()
    {
      return solver_;
    }

    //! Access solver, one of these have to be removed (see above)
    Teuchos::RCP<LINALG::Solver> LinearSolver()
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

    //! Set restart values
    void SetRestart
    (
      int step,   //!< restart step
      double time, //!< restart time
      Teuchos::RCP<Epetra_Vector> disn, //!< restart displacements
      Teuchos::RCP<Epetra_Vector> veln, //!< restart velocities
      Teuchos::RCP<Epetra_Vector> accn, //!< restart accelerations
      Teuchos::RCP<std::vector<char> > elementdata  //!< restart element data
    );

    //! Read and set restart state
    virtual void ReadRestartState();

    //! Sset restart state
    virtual void SetRestartState
    (
      Teuchos::RCP<Epetra_Vector> disn, //!< restart displacements
      Teuchos::RCP<Epetra_Vector> veln, //!< restart velocities
      Teuchos::RCP<Epetra_Vector> accn, //!< restart accelerations
      Teuchos::RCP<std::vector<char> > elementdata  //!< restart element data
    );

    //! Read and set restart forces
    virtual void ReadRestartForce() = 0;

    //! Read and set restart values for constraints
    void ReadRestartConstraint();

    //! Read and set restart values for contact / meshtying
    void ReadRestartContactMeshtying();

    //! Read and set restart values for statmech
    void ReadRestartStatMech();

    //! Read and set restart values for surface stresses
    void ReadRestartSurfstress();

    //! Read and set restart values for multi scale materials
    void ReadRestartMultiScale();

    //! initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() = 0;

    //! right-hand-side of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS() = 0;


    //! @name Access from outside via adapter (needed for coupled problems)
    //@{

    //! unknown displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Dispnp(){ return disn_; }

    //! known displacements at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Dispn(){ return (*dis_)(0); }

    //! unknown velocity at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Velnp() { return veln_; }

    //! unknown velocity at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Veln() { return (*vel_)(0); }

    //! known velocity at \f$t_{n-1}\f$
    Teuchos::RCP<const Epetra_Vector> Velnm() { return (*vel_)(-1); }

    //! unknown accelerations at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Accnp() { return accn_; }

    //! known accelerations at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Accn() { return (*acc_)(0); }

    //@}

    //! Access from inside of the structural time integrator
    //@{

    //! Return material displacements \f$D_{n}\f$
    Teuchos::RCP<Epetra_Vector> Dismat() { return (*dism_)(0); }

    //! Return displacements \f$D_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> DisNew() { return disn_; }

    //! Return displacements \f$D_{n}\f$
    Teuchos::RCP<Epetra_Vector> Dis() { return (*dis_)(0); }

    //! Return velocities \f$V_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> VelNew() { return veln_; }

    //! Return velocities \f$V_{n}\f$
    Teuchos::RCP<Epetra_Vector> Vel() { return (*vel_)(0); }

    //! Return velocities \f$V_{n+1-a_f}\f$
    virtual Teuchos::RCP<Epetra_Vector> Velaf() = 0;

    //! Return accelerations \f$A_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> AccNew() { return accn_; }

    //! Return accelerations \f$A_{n}\f$
    Teuchos::RCP<Epetra_Vector> Acc() { return (*acc_)(0); }

    //@}


    //! Return external force \f$F_{ext,n}\f$
    virtual Teuchos::RCP<Epetra_Vector> Fext() = 0;

    //! Return external force \f$F_{ext,n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> FextNew() = 0;

    //! Return reaction forces
    virtual Teuchos::RCP<Epetra_Vector> Freact() = 0;

    //! Return element data
    //Teuchos::RCP<std::vector<char> > ElementData() {return discret_->PackMyElements();}

    //! dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> DofRowMap();

    //! dof map of vector of unknowns
    // new method for multiple dofsets
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds);

    //! Return stiffness,
    //! i.e. force residual differentiated by displacements
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    //! Return stiffness,
    //! i.e. force residual differentiated by displacements
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

    //! switch structure field to block matrix in fsi simulations
    virtual void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
                                Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps) = 0;

    //! Return sparse mass matrix
    Teuchos::RCP<LINALG::SparseMatrix> MassMatrix();

    //! domain map of system matrix
    virtual const Epetra_Map& DomainMap() const;

    //! are there any algebraic constraints?
    virtual bool HaveConstraint() = 0;

    //! get constraint manager defined in the structure
    virtual const RCP<UTILS::ConstrManager> GetConstraintManager() = 0;

    //! get type of thickness scaling for thin shell structures
    virtual INPAR::STR::STC_Scale GetSTCAlgo() = 0;

    //! Access to scaling matrix for STC
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat() = 0 ;

    //! extract displacements (! not only (*dis_)(0) but all time steps
    virtual Teuchos::RCP<TimIntMStep<Epetra_Vector> > DispMStep() { return dis_; };

    //! times for all the steps
    virtual Teuchos::RCP<TimIntMStep<double> > TimeMStep() { return time_; };


    //@}

    //! @name Time step helpers
    //@{

    //! Return current step number map
    double GetTime() const { return (*time_)[0]; }

    //! Sets the current time \f$t_{n}\f$
    virtual void SetTime(const double time) { (*time_)[0]=time; }

    //! Sets the target time \f$t_{n+1}\f$ of this time step
    virtual void SetTimen(const double time) { timen_ = time; }

    //! Get upper limit of time range of interest
    double GetTimeEnd() const { return timemax_; }

    //! Set upper limit of time range of interest
    void SetTimeEnd(double timemax){timemax_=timemax; }

    //! Get time step size \f$\Delta t_n\f$
    double GetTimeStepSize() const { return (*dt_)[0]; }

    //! Set time step size \f$\Delta t_n\f$
    void SetDt(const double dtnew) { (*dt_)[0] = dtnew; }

    //! Set time step size and time in the first step \f$\Delta t_n\f$
    void SetInitialTimeStepAndTime(const double dt)
    {
    	(*dt_)[0] = dt;
    	timen_ = (*time_)[0] + (*dt_)[0];
    }

    //! Return current step number map
    int GetStep() const { return step_; }

    //! Get current time step \f$n\f$ (sometimes called \f$k\f$ [but not in Baci])
    virtual int GetTimeStep() const { return step_; }

    //! Get number of time steps
    int GetTimeNumStep() const { return stepmax_; }


    //! Return MapExtractor for Dirichlet boundary conditions
    const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return dbcmaps_;
    }

    /*
    \brief Shorten the Dirichlet DOF set

    The method shortens the DOF set (map) which contains the DOFs
    subjected to Dirichlet boundary conditions. For instance, the method is
    called by the monolithic FSI where fluid boundary conditions on the
    FSI interface have to be removed if it is a fluid split.

    */
    void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

    //! Return (rotatory) transformation matrix of local co-ordinate systems
    Teuchos::RCP<const LINALG::SparseMatrix> GetLocSysTrafo() const;

    //! Return locsys manager
    Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() {return locsysman_;}

    //@}

    //! @name Extract methods
    //@{

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp(){ return DisNew(); }

    //! extract displacements
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispn(){ return Dis(); }

    //! extract velocities
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVeln(){ return Vel(); }

    //! extract velocities
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelnp(){ return VelNew(); }

    //@}

    //! @name TSI specific methods
    //@{

    //! specific method for iterative staggered partitioned TSI

    /// Identify residual
    /// This method does not predict the target solution but
    /// evaluates the residual and the stiffness matrix.
    /// In partitioned solution schemes, it is better to keep the current
    /// solution instead of evaluating the initial guess (as the predictor)
    /// does.
    virtual void PreparePartitionStep() = 0;

    /// recalculate structural matrices for tsi simulations
    virtual void TSIMatrix() = 0;

    //@}

    //! @name volume coupling specific methods
    //@{
    //! Set the coupling state (e.g. temperatures for tsi, fluid velocities for poro) onto the structure discretization
    void SetCouplingState();

    //! Apply the coupling state (e.g. temperatures for tsi, fluid velocities for poro) to the structure time integrator
    void ApplyCouplingState
    (
      Teuchos::RCP<const Epetra_Vector> state,  //!< (i) current state
      const std::string& name,                       //!< (i) name of state
      unsigned dofset=1                         //!< (i) number of dofset state is applied to
    );
    //@}

    //! @name Contact and meshtying specific methods
    //@{

    //! use contact solver or not
    bool UseContactSolver();

    //! return bool indicating if contact or meshtying are defined
    bool HaveContactMeshtying() { return (cmtman_ != Teuchos::null); }

    //! return bool indicating if contact or meshtying are defined
    bool HaveSemiSmoothPlasticity() { return (plastman_ != Teuchos::null); }

    //! return contact/meshtying manager
    Teuchos::RCP<MORTAR::ManagerBase> ContactManager() {return cmtman_;}

    //! Check if contact / meshtying is chosen in input file and
    //! create manager object + initialize all relevant stuff if so
    void PrepareContactMeshtying(const Teuchos::ParameterList& sdynparams);

    //! Check if plasticity elements are involved and
    //! create manager object + initialize all relevant stuff if so
    void  PrepareSemiSmoothPlasticity();

    //! Prepare contact at the beginning of each new time step
    //! (call dynamic redistribution of contact interface(s) AND
    //! evaluate reference state for frictional contact at t=0)
    void PrepareStepContact();

    //! Only evaluate reference state for frictional contact at t=0
    void EvaluateReferenceState();

    //@}

    //! @name Beam contact specific methods
    //@{

    //! return bool indicating if beam contact is defined
    bool HaveBeamContact() { return (beamcman_ != Teuchos::null); }

    //! return beam contact manager
    Teuchos::RCP<CONTACT::Beam3cmanager> BeamContactManager() {return beamcman_;}

    //! Check if beam contact is chosen in input file and
    //! create manager object + initialize all relevant stuff if so
    void PrepareBeamContact(const Teuchos::ParameterList& sdynparams);

    //@}

    //@}

    //! @name Statistical Mechanics methods
    //@{

    //! return bool indicating if statmech is defined
    bool HaveStatMech() { return (statmechman_ != Teuchos::null); }

    //! return statmech manager
    Teuchos::RCP<STATMECH::StatMechManager> StatMechManager() {return statmechman_;}

    //! Check if statistical mechanics is chosen in input file and
    //! create manager object + initialize all relevant stuff if so
    void PrepareStatMech();

    //@}


    //! @name Structure with ale specific methods
    //@{

    //! material displacements (structure with ale)
   Teuchos::RCP<Epetra_Vector> DispMat() { return dismatn_; }

    //! apply material displacements to structure field (structure with ale)
    void ApplyDisMat(Teuchos::RCP<Epetra_Vector> dismat);

    //@}

    //! @name Biofilm methods
    //@{
    // reset everything (needed for biofilm simulations)
    void Reset();

    // set structure displacement vector due to biofilm growth
    void SetStrGrDisp(Teuchos::RCP<Epetra_Vector> struct_growth_disp);
    //@}

  protected:

    //! @name General purpose algorithm members
    //@{
    Teuchos::RCP<DRT::Discretization> discret_;  //!< attached discretisation
    int myrank_;  //!< ID of actual processor in parallel
    Teuchos::RCP<LINALG::Solver> solver_;  //!< linear algebraic solver (no contact/meshtying)
    Teuchos::RCP<LINALG::Solver> contactsolver_;  //!< linear algebraic solver (for contact/meshtying)
    bool solveradapttol_;  //!< adapt solver tolerance
    double solveradaptolbetter_;  //!< tolerance to which is adpated ????
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;  //!< map extractor object
                                                  //!< containing non-overlapping
                                                  //!< map of global DOFs on Dirichlet
                                                  //!< boundary conditions
    Teuchos::RCP<Epetra_Vector> dirichtoggle_;    //!< old style dirichtoggle vector
                                                  //!< (only used for contact and supposed
                                                  //!< to go away sonner or later)
    enum INPAR::STR::DivContAct divcontype_ ; // what to do when nonlinear solution fails
    //@}

    //! @name Printing and output
    //@{
    Teuchos::RCP<IO::DiscretizationWriter> output_;  //!< binary output
    bool printlogo_;  //!< true: enjoy your cuppa
    int printscreen_;  //!< print infos to standard out every n steps
    FILE* errfile_;  //!< error file handle
    bool printerrfile_;  //!< print infos to error file #errfile_
    bool printiter_;  //!< print intermediate iterations during solution
    int writerestartevery_;  //!< write restart every given step;
                             //!< if 0, restart is not written
    int writereducedrestart_; // Write reduced Restart Data i.e no mesh for mlmc purposes
    bool writestate_;  //!< write state on/off
    bool writevelacc_;  //!< write velocity and acceleration on/off
    int writeresultsevery_;  //!< write state/stress/strain every given step
    INPAR::STR::StressType writestress_;  //!< stress output type
    INPAR::STR::StressType writecouplstress_;  //!< output type of coupling stress
    INPAR::STR::StrainType writestrain_;  //!< strain output type
    INPAR::STR::StrainType writeplstrain_;  //!< plastic strain output type
    int writeenergyevery_;  //!< write system energy every given step
    bool writesurfactant_;  //!< write surfactant output
    std::ofstream* energyfile_;  //!< outputfile for energy

    Teuchos::RCP<std::vector<char> > stressdata_; //!< container for element GP stresses
    Teuchos::RCP<std::vector<char> > couplstressdata_; //!< container for element GP coupling stresses
    Teuchos::RCP<std::vector<char> > straindata_; //!< container for element GP strains
    Teuchos::RCP<std::vector<char> > plstraindata_; //!< container for element GP plastic strains
    double kinergy_; //!< kinetic energy
    double intergy_; //!< internal energy
    double extergy_; //!< external energy
    //@}

    //! @name Damping
    //!
    //! Rayleigh damping means \f${C} = c_\text{K} {K} + c_\text{M} {M}\f$
    //@{
    enum INPAR::STR::DampKind damping_;  //!< damping type
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
    Teuchos::RCP<POTENTIAL::PotentialManager>   potman_;

    //! contact / meshtying
    Teuchos::RCP<MORTAR::ManagerBase>   cmtman_;

    //! beam contact
    Teuchos::RCP<CONTACT::Beam3cmanager>   beamcman_;

    //! semi-smooth plasticity
    Teuchos::RCP<UTILS::PlastSsnManager>   plastman_;

    //! statistical mechanics
    Teuchos::RCP<STATMECH::StatMechManager>   statmechman_;

    //! Dirichlet BCs with local co-ordinate system
    Teuchos::RCP<DRT::UTILS::LocsysManager> locsysman_;

    //! Map to differentiate pressure and displacement/velocity DOFs
    Teuchos::RCP<LINALG::MapExtractor> pressure_;

    //! elements with micro-materials
    bool havemicromat_;

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
    bool firstoutputofrun_;  //!< flag whether this output step is the first one (restarted or not)
    bool lumpmass_;  //!< flag for lumping the mass matrix, default: false
    bool young_temp_;  //!< flag deciding if young's modulus is temperature-dependent, default: false
    //@}

    //! @name Global vectors
    //@{
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
    //@}

    //! @name Global state vectors
    //@{

    //! global displacements \f${D}_{n}, D_{n-1}, ...\f$
    Teuchos::RCP<TimIntMStep<Epetra_Vector> > dis_;
    //! global material displacements \f${D}_{n}, D_{n-1}, ...\f$
    Teuchos::RCP<TimIntMStep<Epetra_Vector> > dism_;
    //! global velocities \f${V}_{n}, V_{n-1}, ...\f$
    Teuchos::RCP<TimIntMStep<Epetra_Vector> > vel_;
    //! global accelerations \f${A}_{n}, A_{n-1}, ...\f$
    Teuchos::RCP<TimIntMStep<Epetra_Vector> > acc_;
    Teuchos::RCP<Epetra_Vector> disn_;   //!< global displacements
                                         //!< \f${D}_{n+1}\f$
                                         //!< at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> dismatn_;//!< global material displacements
    Teuchos::RCP<Epetra_Vector> veln_;   //!< global velocities
                                         //!< \f${V}_{n+1}\f$
                                         //!< at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> accn_;   //!< global accelerations
                                         //!< \f${A}_{n+1}\f$
                                         //!< at \f$t_{n+1}\f$
    //! additional external forces
    Teuchos::RCP<Epetra_Vector> fifc_;   //!< interface force in e.g. FSI
    //@}


    //! @name System matrices
    //@{
    //! holds eventually effective stiffness
    Teuchos::RCP<LINALG::SparseOperator> stiff_;
    //! mass matrix (constant)
    Teuchos::RCP<LINALG::SparseOperator> mass_;
    //! damping matrix
    Teuchos::RCP<LINALG::SparseOperator> damp_;
    //@}

    //! @name Time measurement
    //@{
    RCP<Epetra_Time>  timer_;           //!< timer for solution technique
    double            dtsolve_;         //!< linear solver time
    double            dtele_;           //!< element evaluation time
    double            dtcmt_;           //!< contact / meshtying evaluation time
    double            inttime_global_;  //!< global integration time for contact evaluation
    //@}

    //! @name volume coupling specific stuff
    //@{
    // global coupling variable
    std::vector<std::map<std::string,Teuchos::RCP<const Epetra_Vector> > > couplstate_; // map of references to coupling states
    //@}

    //! @name Patient specific stuff
    //@{
    Teuchos::RCP<Teuchos::ParameterList> pslist_;
    //@}

    //! @name Biofilm specific stuff
    //@{
    Teuchos::RCP<Epetra_Vector> strgrdisp_;
    //@}

    //! @name porous media specific stuff
    //@{
    Teuchos::RCP<LINALG::MapExtractor>      porositysplitter_;
    //@}
  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_H
