/*----------------------------------------------------------------------*/
/*! \file

\level 1

\brief Time integration for structural dynamics
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef BACI_STRUCTURE_TIMINT_HPP
#define BACI_STRUCTURE_TIMINT_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_config.hpp"

#include "baci_adapter_str_structure.hpp"
#include "baci_inpar_structure.hpp"
#include "baci_lib_elements_paramsinterface.hpp"
#include "baci_timestepping_mstep.hpp"

#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>
#include <ml_common.h>
#include <ml_include.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Time.hpp>

#include <fstream>
#include <iostream>
#include <string>

BACI_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
  class DiscretizationFaces;

  namespace UTILS
  {
    class LocsysManager;
    class PlastSsnManager;
  }  // namespace UTILS
}  // namespace DRT

namespace UTILS
{
  class Cardiovascular0DManager;
}  // namespace UTILS

namespace CONSTRAINTS
{
  class ConstrManager;
  class ConstraintSolver;
  class SpringDashpotManager;
}  // namespace CONSTRAINTS

namespace CONTACT
{
  class Beam3cmanager;
  class MeshtyingContactBridge;
}  // namespace CONTACT

namespace CORE::LINALG
{
  class Solver;
  class MapExtractor;
  class SparseMatrix;
  class SparseOperator;
  class BlockSparseMatrixBase;
}  // namespace CORE::LINALG

namespace IO
{
  class DiscretizationWriter;
}

namespace MOR
{
  class ProperOrthogonalDecomposition;
}

/*----------------------------------------------------------------------*/
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
   * <h3> Responsibilities </h3>
   * Most importantly the base integrator manages the system state vectors and
   * matrices. It also deals with the output to files and offers method to
   * determine forces and stiffnesses (tangents).
   *
   * \author bborn
   * \date 06/08
   */
  class TimInt : public ADAPTER::Structure
  {
    //! Structural time adaptivity is friend
    friend class TimAda;

    //! Joint auxiliary schemes are friends
    template <typename T>
    friend class TimAdaJoint;

   public:
    //! @name Life
    //@{

    //! Print tea time logo
    void Logo();

    //! Constructor
    TimInt(const Teuchos::ParameterList& timeparams,
        const Teuchos::ParameterList& ioparams,            //!< ioflags
        const Teuchos::ParameterList& sdynparams,          //!< input parameters
        const Teuchos::ParameterList& xparams,             //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,          //!< current discretisation
        Teuchos::RCP<CORE::LINALG::Solver> solver,         //!< the solver
        Teuchos::RCP<CORE::LINALG::Solver> contactsolver,  //!< the solver for contact/meshtying
        Teuchos::RCP<IO::DiscretizationWriter> output      //!< the output
    );

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in Init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in Init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref Setup().

    \warning none
    \return bool
    \date 08/16
    \author rauch  */
    virtual void Init(const Teuchos::ParameterList& timeparams,
        const Teuchos::ParameterList& sdynparams, const Teuchos::ParameterList& xparams,
        Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<CORE::LINALG::Solver> solver);

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    void Setup() override;

    //! create fields, based on dofrowmap, whose previous time step values are unimportant
    virtual void CreateFields();

    //! Construct all solution vectors
    void CreateAllSolutionVectors();

    //! Resize \p TimIntMStep<T> multi-step quantities
    virtual void ResizeMStep() = 0;

    //! Resize \p TimIntMStep<T> multi-step quantities, needed for fsi time adaptivity
    void ResizeMStepTimAda() override;

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
    //! and #this is the auxiliary method, which shares the marching data.
    void Merge(TimInt& tis  //!< existing integrator to merge against
    )
    {
      // copy it
      *this = tis;

      return;
    }

    //@}

    //! @name Actions
    //@{

    //! Equilibrate the initial state by identifying the consistent
    //! initial accelerations and (if applicable) internal variables
    //! Make damping and mass matrix
    void DetermineMassDampConsistAccel();

    //! Clear mass matrix and evaluate mass matrix again.
    //! \note not implemented in base class.
    virtual void DetermineMass();

    //! Apply Dirichlet boundary conditions on provided state vectors
    //! (reimplemented in static time integrator)
    virtual void ApplyDirichletBC(const double time,  //!< at time
        Teuchos::RCP<Epetra_Vector> dis,              //!< displacements
                                                      //!< (may be Teuchos::null)
        Teuchos::RCP<Epetra_Vector> vel,              //!< velocities
                                                      //!< (may be Teuchos::null)
        Teuchos::RCP<Epetra_Vector> acc,              //!< accelerations
                                                      //!< (may be Teuchos::null)
        bool recreatemap                              //!< recreate map extractor/toggle-vector
                                                      //!< which stores the DOF IDs subjected
                                                      //!< to Dirichlet BCs
                                                      //!< This needs to be true if the bounded DOFs
                                                      //!< have been changed.
    );

    /// start new time step
    void PrepareTimeStep() override = 0;

    //! Do time integration of multiple steps
    int Integrate() override
    {
      dserror("time loop moved to separate adapter");
      return 0;
    }

    /// tests if there are more time steps to do
    bool NotFinished() const override
    {
      return timen_ <= timemax_ + 1.0e-8 * (*dt_)[0] and stepn_ <= stepmax_;
    }

    //! do something in case nonlinear solution does not converge for some reason
    INPAR::STR::ConvergenceStatus PerformErrorAction(
        INPAR::STR::ConvergenceStatus nonlinsoldiv) override;

    //! Do time integration of single step
    virtual int IntegrateStep() = 0;


    /*! \brief Non-linear solve
     *
     *  Do the nonlinear solve, i.e. (multiple) corrector,
     *  for the time step. All boundary conditions have been set.
     */
    INPAR::STR::ConvergenceStatus Solve() override = 0;

    //! Linear structure solve with just an interface load
    Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear() override = 0;

    /*! \brief Update displacement in case of coupled problems
     *
     *  We always need iterative displacement increments here:
     *
     *  x^n+1_i+1 = x^n+1_i + disiterinc (sometimes referred to as residual increment)
     *
     *  with n and i being time and Newton iteration step
     */
    void UpdateStateIncrementally(Teuchos::RCP<const Epetra_Vector> disiterinc) override = 0;

    /*! \brief Update displacement and evaluate elements in case of coupled problems
     *
     *  We always need iterative displacement increments here:
     *
     *  x^n+1_i+1 = x^n+1_i + disiterinc (sometimes referred to as residual increment)
     *
     *  with n and i being time and Newton iteration step
     */
    void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) override = 0;

    /// don't update displacement but evaluate elements (implicit only)
    void Evaluate() override { dserror("new structural time integration only"); }

    /// update at time step end
    void Update() override = 0;

    /// update at time step end in case of FSI time adaptivity
    void Update(const double endtime) override = 0;

    /// Update iteration
    /// Add residual increment to Lagrange multipliers stored in Constraint manager
    void UpdateIterIncrConstr(
        Teuchos::RCP<Epetra_Vector> lagrincr  ///< Lagrange multiplier increment
        ) override = 0;

    /// output results
    void Output(bool forced_writerestart = false) override = 0;

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awaiting the next
    //! time step.
    virtual void UpdateStepState() = 0;

    //! Update everything on element level after time step and after output
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awaiting the next
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
    //! CAUTION: The VUM is only available for GenAlpha and GEMM.
    void UpdateStepContactVUM();

    //! Update step for beam contact
    void UpdateStepBeamContact();

    //! Reset configuration after time step
    //!
    //! Thus the last converged state is copied back on the predictor
    //! for current time step. This applies only to element-wise
    //! quantities
    void ResetStep() override;

    //! Set initial fields in structure (e.g. initial velocities)
    void SetInitialFields();

    //@}

    //! @name Determination of output quantities that depend on
    // the constitutive model
    //@{
    //! Calculate all output quantities depending on the constitutive model
    //  (and, hence, on a potential material history)
    void PrepareOutput(bool force_prepare_timestep) override;

    //! Calculate stresses, strains on micro-scale
    void PrepareOutputMicro();

    //! Calculate stresses and strains
    void DetermineStressStrain() override;

    //! Calculate kinetic, internal and external energy
    virtual void DetermineEnergy();

    //! Calculate an optional quantity
    void DetermineOptionalQuantity();

    /// create result test for encapsulated structure algorithm
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest() override;

    //@}


    //! @name Output
    //@{

    //! print summary after step
    void PrintStep() override = 0;

    //! Output to file
    //! This routine always prints the last converged state, i.e.
    //! \f$D_{n}, V_{n}, A_{n}\f$. So, #UpdateIncrement should be called
    //! upon object prior to writing stuff here.
    //! \author mwgee (originally) \date 03/07
    void OutputStep(const bool forced_writerestart = false  ///< [in] Force writing of restart data
    );

    bool HasFinalStateBeenWritten() const override;

    //! Write output for every Newton or line search iteration
    //! The step numbers are formated in the following manner:
    //!  n    5               4 2                     1 0
    //!  00..00               000                     00
    //! |__ ___|             |_ _|                   |_ |
    //!    V                   V                       V
    //! digits n to 5       digits 4 to 2            digits 1 to 0
    //! represent the       represent the            represent the
    //! time steps          Newton steps             line search steps
    void OutputEveryIter(bool nw = false, bool ls = false);

    //! write output of step to the Gmsh format
    void WriteGmshStrucOutputStep() override;

    //! Write restart
    //! \author mwgee (originally) \date 03/07
    virtual void OutputRestart(bool& datawritten  //!< (in/out) read and append if
                                                  //!< it was written at this time step
    );
    //! Get data that is written during restart
    //! \author biehler \data 06/13
    void GetRestartData(Teuchos::RCP<int> step, Teuchos::RCP<double> time,
        Teuchos::RCP<Epetra_Vector> disn,  //!< new displacement state
        Teuchos::RCP<Epetra_Vector> veln,  //!< new velocity state
        Teuchos::RCP<Epetra_Vector> accn,  //!< new acceleration state
        Teuchos::RCP<std::vector<char>>
            elementdata,  //!< internal element/history variables e.g. F_prestress
        Teuchos::RCP<std::vector<char>> nodedata  //
        ) override;

    //! Output displacements, velocities and accelerations
    //! and more system vectors
    //! \author mwgee (originally) \date 03/07
    virtual void OutputState(bool& datawritten  //!< (in/out) read and append if
                                                //!< it was written at this time step
    );

    //! Add restart information to OutputState
    void AddRestartToOutputState();

    //! Stress & strain output
    //! \author lw (originally)
    void OutputStressStrain(bool& datawritten  //!< (in/out) read and append if
                                               //!< it was written at this time step
    );

    //! Energy output
    void OutputEnergy();

    //! Optional quantity output
    void OutputOptQuantity(bool& datawritten  //!< (in/out) read and append if
                                              //!< it was written at this time step
    );

    //! Active set, energy and momentum output for contact
    void OutputContact();

    //! Error norm output
    void OutputErrorNorms();

    //! Nodal positions output
    void OutputNodalPositions();

    //! Error norm output
    void OutputVolumeMass();

    //! Output on the micro-scale (multi-scale analysis)
    void OutputMicro();

    //! Write internal and external forces (if necessary for restart)
    virtual void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output) = 0;

    //! Check whether energy output file is attached
    bool AttachedEnergyFile()
    {
      if (not energyfile_.is_null())
        return true;
      else
        return false;
    }

    //! Attach file handle for energy file #energyfile_
    virtual void AttachEnergyFile();

    //@}

    /*! @name Forces
     *
     *  Apply all sets of forces (external, internal, damping, inertia, ...)
     *  based on the current solution state.
     *
     *  On the level of STR::TimInt, we only deal with forces. There are no
     *  stiffnesses since thay are not needed in a general time integration
     *  scheme, but only in an implicit one.
     *
     *  For the application of forces AND
     *  stiffnesses, see STR::TimIntImpl.
     *
     *  \sa STR::TimIntImpl
     */
    //@{

    //! Apply external force
    void ApplyForceExternal(const double time,   //!< evaluation time
        const Teuchos::RCP<Epetra_Vector> dis,   //!< old displacement state
        const Teuchos::RCP<Epetra_Vector> disn,  //!< new displacement state
        const Teuchos::RCP<Epetra_Vector> vel,   // velocity state
        Teuchos::RCP<Epetra_Vector>& fext        //!< external force
    );

    /*! \brief Evaluate ordinary internal force
     *
     *  We need incremental displacements, because the internal variables,
     *  chiefly EAS parameters with an algebraic constraint, are treated
     *  as well. They are not treated perfectly, i.e. they are not iteratively
     *  equilibrated according to their (non-linear) constraint and
     *  the pre-determined displacements -- we talk explicit time integration
     *  here, but they are applied in linearised manner. The linearised
     *  manner means the static condensation is applied once with
     *  residual displacements replaced by the full-step displacement
     *  increment \f$D_{n+1}-D_{n}\f$.
     */
    void ApplyForceInternal(const double time,   //!< evaluation time
        const double dt,                         //!< step size
        Teuchos::RCP<const Epetra_Vector> dis,   //!< displacement state
        Teuchos::RCP<const Epetra_Vector> disi,  //!< incremental displacements
        Teuchos::RCP<const Epetra_Vector> vel,   // velocity state
        Teuchos::RCP<Epetra_Vector> fint         //!< internal force
    );

    //@}

    //! @name Nonlinear mass
    //@{

    //! Return bool indicating if we have nonlinear inertia forces
    int HaveNonlinearMass() const;

    //! check whether the initial conditions are fulfilled */
    virtual void NonlinearMassSanityCheck(
        Teuchos::RCP<const Epetra_Vector> fext,             ///< external forces
        Teuchos::RCP<const Epetra_Vector> dis,              ///< displacements
        Teuchos::RCP<const Epetra_Vector> vel,              ///< velocities
        Teuchos::RCP<const Epetra_Vector> acc,              ///< accelerations
        const Teuchos::ParameterList* sdynparams = nullptr  ///< structural dynamics parameter list
    ) const;

    //@}

    //! Set forces due to interface with fluid, the force is expected external-force-like
    void SetForceInterface(Teuchos::RCP<Epetra_MultiVector> iforce  ///< the force on interface
        ) override;

    //! @name Attributes
    //@{

    //! Provide Name
    virtual enum INPAR::STR::DynamicType MethodName() const = 0;

    //! Provide title
    std::string MethodTitle() const { return INPAR::STR::DynamicTypeString(MethodName()); }

    //! Return true, if time integrator is implicit
    virtual bool MethodImplicit() = 0;

    //! Return true, if time integrator is explicit
    bool MethodExplicit() { return (not MethodImplicit()); }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a \f$m\f$-multistep method returns \f$m\f$
    virtual int MethodSteps() const = 0;

    //! return time integration factor
    double TimIntParam() const override = 0;

    //! Give order of accuracy
    int MethodOrderOfAccuracy() const
    {
      return std::min(MethodOrderOfAccuracyDis(), MethodOrderOfAccuracyVel());
    }

    //! Give local order of accuracy of displacement part
    virtual int MethodOrderOfAccuracyDis() const = 0;

    //! Give local order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() const = 0;

    //! Return linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis() const = 0;

    //! Return linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel() const = 0;

    //@}

    //! @name Access methods
    //@{

    //! Access discretisation
    Teuchos::RCP<DRT::Discretization> Discretization() override
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
      // return Teuchos::rcp(&discret_, false);

      // Now, the discretisation is stored as Teuchos::RefCountPointer,
      // thus
      return discret_;
    }

    //! Access to dofrowmap of discretization via const raw pointer
    const Epetra_Map* DofRowMapView() override;

    //! Access solver, one of these have to be removed (see below)
    Teuchos::RCP<CORE::LINALG::Solver> Solver() { return solver_; }

    //! Access solver, one of these have to be removed (see above)
    Teuchos::RCP<CORE::LINALG::Solver> LinearSolver() override { return solver_; }

    //! Access solver for contact/meshtying problems
    Teuchos::RCP<CORE::LINALG::Solver> ContactSolver() { return contactsolver_; }

    //! Access output object
    Teuchos::RCP<IO::DiscretizationWriter> DiscWriter() override { return output_; }

    //! Read restart values
    void ReadRestart(const int step  //!< restart step
        ) override;

    //! Set restart values
    void SetRestart(int step,                         //!< restart step
        double time,                                  //!< restart time
        Teuchos::RCP<Epetra_Vector> disn,             //!< restart displacements
        Teuchos::RCP<Epetra_Vector> veln,             //!< restart velocities
        Teuchos::RCP<Epetra_Vector> accn,             //!< restart accelerations
        Teuchos::RCP<std::vector<char>> elementdata,  //!< restart element data
        Teuchos::RCP<std::vector<char>> nodedata      //!< restart element data
        ) override;

    //! Set the state of the nox group and the global state data container (implicit only)
    void SetState(const Teuchos::RCP<Epetra_Vector>& x) override
    {
      dserror("new structural time integration only...");
    }

    //! Read and set restart state
    virtual void ReadRestartState();

    //! Set restart state
    virtual void SetRestartState(Teuchos::RCP<Epetra_Vector> disn,  //!< restart displacements
        Teuchos::RCP<Epetra_Vector> veln,                           //!< restart velocities
        Teuchos::RCP<Epetra_Vector> accn,                           //!< restart accelerations
        Teuchos::RCP<std::vector<char>> elementdata,                //!< restart element data
        Teuchos::RCP<std::vector<char>> nodedata                    //!< restart element data
    );

    //! Read and set restart forces
    virtual void ReadRestartForce() = 0;

    //! Read and set restart values for constraints
    void ReadRestartConstraint();

    //! Read and set restart values for Cardiovascular0D
    void ReadRestartCardiovascular0D();

    //! Read and set restart values for Spring Dashpot
    void ReadRestartSpringDashpot();

    //! Read and set restart values for contact / meshtying
    void ReadRestartContactMeshtying();

    //! Read and set restart values for beam contact
    void ReadRestartBeamContact();

    //! Read and set restart values for multi scale materials
    void ReadRestartMultiScale();

    //! initial guess of Newton's method
    Teuchos::RCP<const Epetra_Vector> InitialGuess() override = 0;

    //! right-hand-side of Newton's method
    Teuchos::RCP<const Epetra_Vector> RHS() override = 0;

    /// set evaluation action
    void SetActionType(const DRT::ELEMENTS::ActionType& action) override
    {
      dserror("new structural time integration only...");
    }

    //! @name Access from outside via adapter (needed for coupled problems)
    //@{

    //! unknown displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Dispnp() const override { return disn_; }

    //! known displacements at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Dispn() const override { return (*dis_)(0); }

    //! unknown velocity at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Velnp() const override { return veln_; }

    //! unknown velocity at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Veln() const override { return (*vel_)(0); }

    //! known velocity at \f$t_{n-1}\f$
    Teuchos::RCP<const Epetra_Vector> Velnm() const override { return (*vel_)(-1); }

    //! unknown accelerations at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Accnp() const override { return accn_; }

    //! known accelerations at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Accn() const override { return (*acc_)(0); }

    //@}


    //! Access from inside of the structural time integrator
    //@{

    //! Return material displacements \f$D_{n}\f$
    Teuchos::RCP<Epetra_Vector> Dismat() { return (*dismat_)(0); }

    //! Return displacements \f$D_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> DisNew() { return disn_; }

    //! Return displacements \f$D_{n}\f$
    Teuchos::RCP<Epetra_Vector> Dis() { return (*dis_)(0); }

    //! Return velocities \f$V_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> VelNew() { return veln_; }

    //! Return velocities \f$V_{n}\f$
    Teuchos::RCP<Epetra_Vector> Vel() { return (*vel_)(0); }

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
    Teuchos::RCP<Epetra_Vector> Freact() override = 0;

    //! Return element data
    // Teuchos::RCP<std::vector<char> > ElementData() {return discret_->PackMyElements();}

    //! dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> DofRowMap() override;

    //! dof map of vector of unknowns
    // new method for multiple dofsets
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) override;

    //! Return stiffness,
    //! i.e. force residual differentiated by displacements
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override;

    //! Return stiffness,
    //! i.e. force residual differentiated by displacements
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override;

    //! switch structure field to block matrix in fsi simulations
    void UseBlockMatrix(Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> domainmaps,
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> rangemaps) override = 0;

    //! Return sparse mass matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> MassMatrix();

    //! domain map of system matrix
    const Epetra_Map& DomainMap() const override;

    //! are there any algebraic constraints?
    bool HaveConstraint() override = 0;

    //! are there any spring dashpot BCs?
    bool HaveSpringDashpot() override = 0;

    //! get constraint manager defined in the structure
    Teuchos::RCP<CONSTRAINTS::ConstrManager> GetConstraintManager() override = 0;

    //! get spring dashpot manager defined in the structure
    Teuchos::RCP<CONSTRAINTS::SpringDashpotManager> GetSpringDashpotManager() override = 0;

    //! get type of thickness scaling for thin shell structures
    INPAR::STR::STC_Scale GetSTCAlgo() override = 0;

    //! Access to scaling matrix for STC
    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetSTCMat() override = 0;


    //@}

    //! @name Time step helpers
    //@{

    //! Return current time \f$t_{n}\f$
    double TimeOld() const override { return (*time_)[0]; }

    //! Return target time \f$t_{n+1}\f$
    double Time() const override { return timen_; }

    //! Sets the current time \f$t_{n}\f$
    void SetTime(const double time) override { (*time_)[0] = time; }

    //! Sets the target time \f$t_{n+1}\f$ of this time step
    void SetTimen(const double time) override { timen_ = time; }

    //! Sets the current step \f$n+1\f$
    void SetStep(int step) override { step_ = step; }

    //! Sets the current step \f$n+1\f$
    void SetStepn(int step) override { stepn_ = step; }

    //! Get upper limit of time range of interest
    double GetTimeEnd() const override { return timemax_; }

    //! Set upper limit of time range of interest
    void SetTimeEnd(double timemax) override { timemax_ = timemax; }

    //! Get time step size \f$\Delta t_n\f$
    double Dt() const override { return (*dt_)[0]; }

    //! Set time step size \f$\Delta t_n\f$
    void SetDt(const double dtnew) override { (*dt_)[0] = dtnew; }

    //! Return current step number $n$
    int StepOld() const override { return step_; }

    //! Return current step number $n+1$
    int Step() const override { return stepn_; }

    //! Get number of time steps
    int NumStep() const override { return stepmax_; }


    //! Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor() override
    {
      return dbcmaps_;
    }

    //! Return (rotatory) transformation matrix of local co-ordinate systems
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> GetLocSysTrafo() const;

    //! Return locsys manager
    Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() override { return locsysman_; }

    //@}

    //! @name Write access to field solution variables at \f$t^{n+1}\f$
    //@{

    /// write access to displacements at \f$t^{n+1}\f$
    Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() override { return DisNew(); }

    //! write access to velocities at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> WriteAccessVelnp() override { return VelNew(); }

    /// write access to displacements at \f$t^{n+1}\f$
    Teuchos::RCP<Epetra_Vector> WriteAccessDispn() override { return Dis(); }

    //! write access to velocities at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> WriteAccessVeln() override { return Vel(); }

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
    void PreparePartitionStep() override = 0;

    //@}

    //! @name Contact and meshtying specific methods
    //@{

    //! return bool indicating if contact or meshtying are defined
    bool HaveContactMeshtying() { return (cmtbridge_ != Teuchos::null); }

    //! return contact/meshtying manager
    Teuchos::RCP<CONTACT::MeshtyingContactBridge> MeshtyingContactBridge() override
    {
      return cmtbridge_;
    }

    /// do we have this model
    bool HaveModel(INPAR::STR::ModelType model) override
    {
      dserror("new structural time integration only");
      return false;
    }

    STR::MODELEVALUATOR::Generic& ModelEvaluator(INPAR::STR::ModelType mtype) override
    {
      dserror("new time integration only");
      exit(EXIT_FAILURE);
    }

    /*!
    \brief Prepare time integration for contact/meshtying

    Check if contact / meshtying is chosen in input file. If yes, create manager object and
    initialize all relevant stuff.

    @param[in] sdynparams Structural dynamics input parameter list
    */
    void PrepareContactMeshtying(const Teuchos::ParameterList& sdynparams);

    /*!
    \brief Apply results of mesh initialization to the underlying problem discretization

    \note This is only necessary in case of a mortar method.

    \warning This routine modifies the reference coordinates of slave nodes at the meshtying
    interface.

    @param[in] Xslavemod Vector with modified nodal positions
    */
    void ApplyMeshInitialization(Teuchos::RCP<const Epetra_Vector> Xslavemod);

    //! Prepare contact at the beginning of each new time step
    //! (call dynamic redistribution of contact interface(s) AND
    //! evaluate reference state for frictional contact at t=0)
    void PrepareStepContact();

    /// wrapper for things that should be done before PrepareTimeStep is called
    void PrePredict() final{};

    /// wrapper for things that should be done before solving the nonlinear iterations
    void PreSolve() final{};

    /// wrapper for things that should be done before updating
    void PreUpdate() final{};

    /// wrapper for things that should be done after solving the update
    void PostUpdate() final{};

    /// wrapper for things that should be done after convergence of Newton scheme
    void PostOutput() final{};

    /// wrapper for things that should be done after the actual time loop is finished
    void PostTimeLoop() final;

    //@}

    //! @name Beam contact specific methods
    //@{

    //! return bool indicating if beam contact is defined
    bool HaveBeamContact() { return (beamcman_ != Teuchos::null); }

    //! return beam contact manager
    Teuchos::RCP<CONTACT::Beam3cmanager> BeamContactManager() { return beamcman_; }

    //! Check if beam contact is chosen in input file and
    //! create manager object + initialize all relevant stuff if so
    void PrepareBeamContact(const Teuchos::ParameterList& sdynparams);

    //@}

    //! @name Structure with ale specific methods
    //@{

    //! material displacements (structure with ale)
    Teuchos::RCP<Epetra_Vector> DispMat() override { return dismatn_; }

    //! apply material displacements to structure field (structure with ale)
    void ApplyDisMat(Teuchos::RCP<Epetra_Vector> dismat) override;

    //@}

    //! @name Biofilm methods
    //@{

    // reset everything (needed for biofilm simulations)
    void Reset() override;

    // set structure displacement vector due to biofilm growth
    void SetStrGrDisp(Teuchos::RCP<Epetra_Vector> struct_growth_disp) override;

    virtual bool HaveBiofilmGrowth() const { return (not strgrdisp_.is_null()); }

    //@}

    //! @name Micro material methods
    //@{

    //! bool indicating if micro material is used
    bool HaveMicroMat() override { return havemicromat_; }

    //}

   protected:
    /// Expand the dbc map by dofs provided in Epetra_Map maptoadd
    void AddDirichDofs(const Teuchos::RCP<const Epetra_Map> maptoadd) override;

    /// Contract the dbc map by dofs provided in Epetra_Map maptoremove
    void RemoveDirichDofs(const Teuchos::RCP<const Epetra_Map> maptoremove) override;

    //! @name General purpose algorithm members
    //@{
    Teuchos::RCP<DRT::Discretization> discret_;  //!< attached discretisation

    // face discretization (only initialized for edge-based stabilization)
    Teuchos::RCP<DRT::DiscretizationFaces> facediscret_;

    int myrank_;                                 //!< ID of actual processor in parallel
    Teuchos::RCP<CORE::LINALG::Solver> solver_;  //!< linear algebraic solver (no contact/meshtying)
    Teuchos::RCP<CORE::LINALG::Solver>
        contactsolver_;           //!< linear algebraic solver (for contact/meshtying)
    bool solveradapttol_;         //!< adapt solver tolerance
    double solveradaptolbetter_;  //!< tolerance to which is adapted ????
    Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps_;  //!< map extractor object
                                                        //!< containing non-overlapping
                                                        //!< map of global DOFs on Dirichlet
                                                        //!< boundary conditions

    enum INPAR::STR::DivContAct divcontype_;  //!< what to do when nonlinear solution fails
    int divconrefinementlevel_;  //!< number of refinement level in case of divercontype_ ==
                                 //!< adapt_step
    int divconnumfinestep_;      //!< number of converged time steps on current refinement level
                                 //!< in case of divercontype_ == adapt_step

    //! structural dynamic parameter list
    Teuchos::ParameterList sdynparams_;

    //@}

    //! @name Printing and output
    //@{
    Teuchos::RCP<IO::DiscretizationWriter> output_;  //!< binary output
    int printscreen_;                                //!< print infos to standard out every n steps
    bool printlogo_;                                 //!< print the logo (or not)?
    bool printiter_;                           //!< print intermediate iterations during solution
    bool outputeveryiter_;                     //!< switch
    int oei_filecounter_;                      //!< filename counter
    int outputcounter_;                        //!< output counter for OutputEveryIter
    int writerestartevery_;                    //!< write restart every given step;
                                               //!< if 0, restart is not written
    bool writeele_;                            //!< write elements on/off
    bool writestate_;                          //!< write state on/off
    bool writevelacc_;                         //!< write velocity and acceleration on/off
    int writeresultsevery_;                    //!< write state/stress/strain every given step
    INPAR::STR::StressType writestress_;       //!< stress output type
    INPAR::STR::StressType writecouplstress_;  //!< output type of coupling stress
    INPAR::STR::StrainType writestrain_;       //!< strain output type
    INPAR::STR::StrainType writeplstrain_;     //!< plastic strain output type
    INPAR::STR::OptQuantityType writeoptquantity_;  //!< stress output type
    int writeenergyevery_;                          //!< write system energy every given step
    bool writesurfactant_;                          //!< write surfactant output
    bool writerotation_;                            //!< write strutural rotation tensor output
    Teuchos::RCP<std::ofstream> energyfile_;        //!< outputfile for energy

    Teuchos::RCP<std::vector<char>> stressdata_;  //!< container for element GP stresses
    Teuchos::RCP<std::vector<char>>
        couplstressdata_;                           //!< container for element GP coupling stresses
    Teuchos::RCP<std::vector<char>> straindata_;    //!< container for element GP strains
    Teuchos::RCP<std::vector<char>> plstraindata_;  //!< container for element GP plastic strains
    Teuchos::RCP<std::vector<char>> rotdata_;       //!< container for element rotation tensor
    Teuchos::RCP<std::vector<char>>
        optquantitydata_;  //!< container for element GP optional quantities
    double kinergy_;       //!< kinetic energy
    double intergy_;       //!< internal energy
    double extergy_;       //!< external energy
    //@}

    //! @name Damping
    //!
    //! Rayleigh damping means \f${C} = c_\text{K} {K} + c_\text{M} {M}\f$
    //@{
    enum INPAR::STR::DampKind damping_;  //!< damping type
    double dampk_;                       //!< damping factor for stiffness \f$c_\text{K}\f$
    double dampm_;                       //!< damping factor for mass \f$c_\text{M}\f$
    //@}

    //! @name Managed stuff
    //@{

    //! whatever constraints
    Teuchos::RCP<CONSTRAINTS::ConstrManager> conman_;      //!< constraint manager
    Teuchos::RCP<CONSTRAINTS::ConstraintSolver> consolv_;  //!< constraint solver

    // for 0D cardiovascular models
    Teuchos::RCP<UTILS::Cardiovascular0DManager> cardvasc0dman_;  //!< Cardiovascular0D manager

    // spring dashpot
    Teuchos::RCP<CONSTRAINTS::SpringDashpotManager> springman_;

    //! bridge for meshtying and contact
    Teuchos::RCP<CONTACT::MeshtyingContactBridge> cmtbridge_;

    //! beam contact
    Teuchos::RCP<CONTACT::Beam3cmanager> beamcman_;

    //! Dirichlet BCs with local co-ordinate system
    Teuchos::RCP<DRT::UTILS::LocsysManager> locsysman_;

    //! Map to differentiate pressure and displacement/velocity DOFs
    Teuchos::RCP<CORE::LINALG::MapExtractor> pressure_;

    //! elements with micro-materials
    bool havemicromat_;

    //! Is GMSH output of displacements required?
    bool gmsh_out_;

    //@}

    //! @name General control parameters
    //@{
    Teuchos::RCP<TIMESTEPPING::TimIntMStep<double>>
        time_;      //!< time \f$t_{n}\f$ of last converged step
    double timen_;  //!< target time \f$t_{n+1}\f$
    Teuchos::RCP<TIMESTEPPING::TimIntMStep<double>> dt_;  //!< time step size \f$\Delta t\f$
    double timemax_;                                      //!< final time \f$t_\text{fin}\f$
    int stepmax_;                                         //!< final step \f$N\f$
    int step_;                                            //!< time step index \f$n\f$
    int stepn_;                                           //!< time step index \f$n+1\f$
    double rand_tsfac_;      //!< random factor for modifying time-step size in case this way of
                             //!< continuing non-linear iteration was chosen
    bool firstoutputofrun_;  //!< flag whether this output step is the first one (restarted or not)
    bool lumpmass_;          //!< flag for lumping the mass matrix, default: false
    //@}

    //! @name Global vectors
    //@{
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
    //@}

    //! @name Global state vectors
    //@{

    //! global displacements \f${D}_{n}, D_{n-1}, ...\f$
    Teuchos::RCP<TIMESTEPPING::TimIntMStep<Epetra_Vector>> dis_;

    //! global material displacements \f${D}_{n}, D_{n-1}, ...\f$
    Teuchos::RCP<TIMESTEPPING::TimIntMStep<Epetra_Vector>> dismat_;

    //! global velocities \f${V}_{n}, V_{n-1}, ...\f$
    Teuchos::RCP<TIMESTEPPING::TimIntMStep<Epetra_Vector>> vel_;

    //! global accelerations \f${A}_{n}, A_{n-1}, ...\f$
    Teuchos::RCP<TIMESTEPPING::TimIntMStep<Epetra_Vector>> acc_;

    //!< global displacements \f${D}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> disn_;

    //!< global material displacements
    Teuchos::RCP<Epetra_Vector> dismatn_;  //!< global material displacements

    //!< global velocities \f${V}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> veln_;

    //!< global accelerations \f${A}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> accn_;

    //!< global internal force
    Teuchos::RCP<Epetra_Vector> fint_;

    //! additional external forces (e.g. interface force in FSI)
    Teuchos::RCP<Epetra_Vector> fifc_;

    //!< pure structural global internal force, i.e. no condensation of EAS, plasticity,...
    Teuchos::RCP<Epetra_Vector> fresn_str_;

    //!< pure structural global internal force at \f$t_n\f$, i.e. no condensation of EAS,
    //!< plasticity,...
    Teuchos::RCP<Epetra_Vector> fintn_str_;

    //@}


    //! @name System matrices
    //@{
    //! holds eventually effective stiffness
    Teuchos::RCP<CORE::LINALG::SparseOperator> stiff_;

    //! mass matrix (constant)
    Teuchos::RCP<CORE::LINALG::SparseOperator> mass_;

    //! damping matrix
    Teuchos::RCP<CORE::LINALG::SparseOperator> damp_;
    //@}

    //! @name Time measurement
    //@{
    Teuchos::RCP<Teuchos::Time> timer_;  //!< timer for solution technique
    double dtsolve_;                     //!< linear solver time
    double dtele_;                       //!< element evaluation time
    double dtcmt_;                       //!< contact / meshtying evaluation time
    double inttime_global_;              //!< global integration time for contact evaluation
    //@}

    //! @name Biofilm specific stuff
    //@{
    Teuchos::RCP<Epetra_Vector> strgrdisp_;
    //@}

    //! @name porous media specific stuff
    //@{
    Teuchos::RCP<CORE::LINALG::MapExtractor> porositysplitter_;
    //@}

    Teuchos::RCP<MOR::ProperOrthogonalDecomposition> mor_;  //!< model order reduction

   private:
    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

    //! load/time step of the last written results
    int lastwrittenresultsstep_;

   protected:
    //! returns true if Setup() was called and is still valid
    bool IsSetup() { return issetup_; };

    //! returns true if Init(..) was called and is still valid
    bool IsInit() const { return isinit_; };

    //! check if \ref Setup() was called
    void CheckIsSetup()
    {
      if (not IsSetup()) dserror("Setup() was not called.");
    };

    //! check if \ref Init() was called
    void CheckIsInit() const
    {
      if (not IsInit()) dserror("Init(...) was not called.");
    };

   public:
    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };

  };  // class TimInt

}  // namespace STR

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // STRUCTURE_TIMINT_H
