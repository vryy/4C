/*----------------------------------------------------------------------*/
/*!
\file strutimint.H

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
#ifndef STRUTIMINT_H
#define STRUTIMINT_H

/*----------------------------------------------------------------------*/
/* headers */
#include <iostream>
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include "../drt_io/io.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_constraint/uzawasolver.H"
#include "../drt_surfstress/drt_surfstress_manager.H"
#include "../drt_surfstress/drt_potential_manager.H"
#include "../drt_contact/drt_contact_manager.H"

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

#include "strutimint_vector.H"
#include "strutimint_state.H"
//#include "strutimada.H"

/*======================================================================*/
/*!
\brief Front-end for structural dynamics by structure time integration. 
\author bborn
\date 06/08
*/
namespace STR  // belongs to structure namespace
{

  class StruTimInt
  {

    //! Structural time adaptivity is friend
    friend class StruTimAda;

  public:

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
      name_genalpha,  //!< generalised-alpha
      name_onesteptheta,  //!< one-step-theta
      name_ab2  //!< Adams-Bashforth 2nd order
    };

    //! Map name enum term to string
    static std::string MapNameEnumToString
    (
      const enum NameEnum term  //!< the enum
    )
    {
      switch (term)
      {
      case name_genalpha :
        return "GenAlpha";
        break;
      case name_onesteptheta :
        return "OneStepTheta";
        break;
      case name_ab2 :
        return "AdamsBashforth2";
        break;
      default :
        dserror("Cannot cope with name enum %d", term);
        return "";
        break;
      }
    }   

    //! Constructor
    StruTimInt
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      DRT::Discretization& actdis,  //!< current discretisation
      LINALG::Solver& solver,  //!< the solver
      IO::DiscretizationWriter& output  //!< the output
    );

    //! Destructor
    virtual ~StruTimInt() { ; }

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
      Teuchos::RCP<Epetra_Vector> acc  //!< accelerations 
                                       //!< (may be Teuchos::null)
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

    //! print summary after step
    virtual void PrintStep() = 0;

    //! Output to file
    virtual void OutputStep() = 0;

    //@}

    //! @name Forces and Stiffnesses
    //@{

    //! Apply external force
    void ApplyForceExternal
    (
      const double time,  //!< evaluation time
      const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
      Teuchos::RCP<Epetra_Vector>& fext  //!< external force
    );

    //! Evaluate ordinary internal force, its stiffness at state
    void ApplyForceStiffInternal
    (
      const double time,  //!< evaluation time
      const double dt,  //!< step size
      const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
      const Teuchos::RCP<Epetra_Vector> disi,  //!< residual displacements
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

    //! Give local order of accuracy of displacement part
    virtual int MethodOrderOfAccuracyDis() = 0;

    //! Give local order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() = 0;

    //! Access dicretisation
    Teuchos::RCP<DRT::Discretization> Discretization()
    { 
      // Here a 'false' must be used. This is due to
      // the fact that StruTimInt possesses a references
      // on the discretisation #discret_ and not
      // a Teuchos::RefCountPointer. Eventually, StruTimInt
      // will be destroyed and it will immediately destroy
      // its #discret_ member. However, #discret_ is handed down
      // to the ConstrManager and kept there as a RefCountPointer.
      // The object #discret_ is gone, when ConstrManager tries
      // to kill it. We achieve a nice segmentation fault.
      // The 'false' prevents ConstrManager of trying to kill it.
      return Teuchos::rcp(&discret_, false);
    }

    //! Access solver
    Teuchos::RCP<LINALG::Solver> GetSolver()
    {
      return Teuchos::rcp(&solver_, false);
    }

    //! Access output object
    Teuchos::RCP<IO::DiscretizationWriter> GetDiscretizationWriter()
    {
      return Teuchos::rcp(&output_, false);
    }

    //@}

    //! @name Phase states
    //@{

    //! Return displacements \f$D_{n}\f$
    Teuchos::RCP<Epetra_Vector> Disp()
    {
      return dis_();
    }

    //! Return velocities \f$V_{n}\f$
    Teuchos::RCP<Epetra_Vector> Vel()
    {
      return vel_();
    }

    //! Return accelerations \f$A_{n}\f$
    Teuchos::RCP<Epetra_Vector> Acc()
    {
      return acc_();
    }

    //@}

  protected:

    //! @name General purpose algorithm members
    //@{
    DRT::Discretization& discret_;  //!< attached discretisation
    int myrank_;  //!< ID of actual processor in parallel
    const Epetra_Map* dofrowmap_;  //!< vector layout from the discretization 
    //!< to construct matching
    //!< vectors and matrices
    LINALG::Solver& solver_;  //!< linear algebraic solver
    bool solveradapttol_;  //!< adapt solver tolerance
    double solveradaptolbetter_;  //!< tolerance to which is adapated ????
    //@}

    //! @name Printing and output
    //@{
    IO::DiscretizationWriter& output_;  //!< binary output
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
    //@}

    //! @name Rayleigh damping
    //!       \f${C} = c_\text{K} {K} + c_\text{M} {M}\f$
    //@{
    bool damping_;  //!< damping active/inactive
    double dampk_;  //!< damping factor for stiffness \f$c_\text{K}\f$
    double dampm_;  //!< damping factor for mass \f$c_\text{M}\f$
    //@}

    //! @name Managed stuff
    //@{

    //! whatever constraints
    Teuchos::RCP<ConstrManager> conman_;  //!< constraint manager
    Teuchos::RCP<UzawaSolver> uzawasolv_;  //!< Uzawa constraint solver

    //! surface stresses
    Teuchos::RCP<DRT::SurfStressManager> surfstressman_;

    //! potential conditions
    Teuchos::RCP<DRT::PotentialManager> potman_;

    //@}

    //! @name General control parameters
    //@{
    double time_;  //!< time \f$t_{n}\f$ of last converged step
    double timen_;  //!< target time \f$t_{n+1}\f$
    double dt_;  //!< time step size \f$\Delta t\f$
    double timemax_;  //!< final time \f$t_\text{fin}\f$
    int stepmax_;  //!< final step \f$N\f$
    int step_;  //!< time step index \f$n\f$
    int stepn_;  //!< time step index \f$n+1\f$
    //@}

    //! @name Global vectors
    //@{
    Teuchos::RCP<Epetra_Vector> dirichtoggle_;  //!< At a free DOF this 
                                                //!< vector holds a 0; 
                                                //!< at a supported DOF
                                                //!< it is a 1.
    Teuchos::RCP<Epetra_Vector> invtoggle_;  //!< Opposite of 
                                             //!< #dirichtoggle_ vector:
                                             //!< i.e.\ for each component 
                                             //!< a 0 indicates a supported 
                                             //!< corresponding DOF;
                                             //!< a 1 indicates
                                             //!< a free/unsupported DOF.
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
    //@}


    //! @name Global state vectors
    //@{

    StruTimIntState dis_;  //!< global displacements
                           //!< \f${D}_{n}, D_{n-1}, ...\f$
    StruTimIntState vel_;  //!< global velocities
                           //!< \f${V}_{n}, V_{n-1}, ...\f$
    StruTimIntState acc_;  //!< global accelerations
                           //!< \f${A}_{n}, A_{n-1}, ...\f$
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
#endif  // #ifndef STRUTIMINT_H
#endif  // #ifdef CCADISCRET
