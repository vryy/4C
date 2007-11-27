/*----------------------------------------------------------------------*/
/*!
\file mfsi_algorithm.H

\brief Solve FSI problems in a monolithic solver

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef MFSI_ALGORITHM_H
#define MFSI_ALGORITHM_H

#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Thyra_LinearOpWithSolveFactoryBase.hpp>
#include <Thyra_PreconditionerFactoryBase.hpp>
#include <Thyra_StateFuncModelEvaluatorBase.hpp>
#include <Thyra_VectorSpaceBase.hpp>

#include <Thyra_DefaultProductVectorSpace.hpp>
#include <Thyra_DefaultProductVector.hpp>
#include <Thyra_DefaultBlockedLinearOp.hpp>
#include <Thyra_DefaultLinearOpSource.hpp>

#include <NOX_Thyra_Group.H>
#include <NOX_Thyra_Vector.H>

#include "mfsi_structure.H"
#include "mfsi_fluid.H"
#include "mfsi_coupling.H"
#include "mfsi_debug.H"

#include "../drt_fsi/fsi_ale.H"


/// monolithic FSI
/*!

A highly experimental FSI solution method that solves the whole
nonlinear system in one huge Newton method. The implementation is
based on NOX, the nonlinear solver, and uses thyra to glue things
together. The Algorithm class is the main class of the monolithic
solver and contains further explanation.

\author u.kue
\date 11/07
 */
namespace MFSI {

  /// monolithic FSI algorithm class
  /*!

    Here we have it. Monolithic FSI based on independent field
    solvers. This class makes heavy use of Thyra to glue the linear
    systems of the three participating fields together. The result is
    a global linear block system at the core of a (NOX based) Newton
    scheme. The composite system really is meant to be solved by a
    Newton method. Other methods might follow later.

    Communication with the field algorithms is done via adapter
    classes that constrain the interfaces the algorithms present.

    The central point of the monolithic FSI algorithm is the
    Timeloop() method. Here the global NOX solver is constructed in
    each time step. We use the thyra interface to NOX and present
    ourself as thyra model evaluator to NOX.

    The point of being a model evaluator is the ability to calculate
    r(x) and J(x), where J(x)=r'(x). In addition to that the model
    evaluator can tell what argument it expects and what results it
    produces. We use a very basic version here. We expect just a
    vector x from NOX and return the rhs r or the system matrix
    J. There is much more to the thyra model evaluator, but right now
    that is all we need. (And it is all NOX can handle, btw.)

    We have to cheat a little, however, as NOX will ask us of the
    right hand side and for the system matrix in separate calls,
    however our element calls always give us both at the same time and
    we do not want to call the elements twice. That is why there is a
    NOX group subclass NOX_Thyra_Group here, that implicitly updates
    both variables. A detail only.

    So NOX does not know anything about a block system but applies the
    Newton method to a nonlinear system given with thyra objects that
    can give its linearization.

    Underneath, however, the thyra objects are block objects composed
    of three field equations and two coupling conditions applied in
    Lagrangian manner. The coupling conditions are formulated in terms
    of the structures interface map. Right now only matching meshes
    are supported. But it is simple to insert any mortar coupling
    matrix.

    One point with the combined monolithic system is the linear
    solver. We are bound to iterative solvers here. And we need
    suitable block preconditioners to solve the system. There is the
    preconditioner factory that gives us a very nice preconditioner
    based on linear solves of the field equations. However, the
    coupling conditions introduce a large disturbance. A perfect
    preconditioner is still an open issue.

    On a finer level there are a lot of details to get right. All
    fields have to be glued together in the right way. So some
    knowledge of what the field solver actually does goes into the
    coupling as well.

    Here, both fluid and structural solver solve for increments. After
    each Newton step they want to get the next delta. NOX, on the
    other hand, accumulates all increments over a nonlinar solve and
    hands in the summed value.

    Another issue is the convergance test. We need to check the
    convergance of structural displacement equations, fluid velocity
    equations, fluid pressure equations and interface equations, where
    structural displacements and fluid velocities meet. The last ones
    should be at the same point in time, probably at the end of the
    time step.

    \warning Right now the coupling is wrong unless backward Euler is
    used in both fluid and structural field.

    \note This class binds two levels of detail together. On the lower
    level, the field and coupling level, epetra objects are used
    exclusively. There is no numerical work going on that is not
    handled by epetra. On the higher level, the nonlinear solver
    (NOX), its linear solver and preconditioner, only thyra objects
    are used. Thyra is needed, because it provides the ability to
    build block systems very easily. But these thyra objects are just
    wrappers around epetra objects. Changing the values of a thyra
    object will always affect its epetra object. You have to
    understand the relation between thyra and epetra to work with this
    implementation.

    \author u.kue
    \date 11/07
  */
  class Algorithm : public Thyra::StateFuncModelEvaluatorBase<double>
  {
  public:
    Algorithm(Epetra_Comm& comm);

    /// The time loop. Start this. Everything is in here.
    void Timeloop();

    //! @name  StateFuncModelEvaluatorBase
    /// methods needed to implement the thyra model evaluator

    /// tell the map of our vector of unknowns
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > get_x_space() const;

    /// tell the map of our right hand side
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > get_f_space() const;

    /// give a description of the arguments this model excepts (just x)
    Thyra::ModelEvaluatorBase::InArgs<double> createInArgs() const;

    /// give a description of what this model can compute (r and J)
    Thyra::ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const;

    /// actually call the model and calculate whatever is requested
    void evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
                       const Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs) const;

    /// return the linear solver to solve this model's linear problem
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > create_W() const;

    /// return the linear problem (system matrix) of this model
    Teuchos::RCP<Thyra::LinearOpBase<double> > create_W_op() const;

    //@}

    /// Evaluate all fields at position x
    void Evaluate(Teuchos::RCP<const Thyra::DefaultProductVector<double> > x) const;

    /// Extract current unknown from fields
    //void CurrentX(Thyra::DefaultProductVector<double>& x);

    /// Extract initial guess from fields
    void InitialGuess(Thyra::DefaultProductVector<double>& ig);

    // Apply current field state to system

    /// setup composed right hand side from field solvers
    void SetupRHS(Thyra::DefaultProductVector<double> &f) const;

    /// setup composed system matrix from field solvers
    void SetupSysMat(Thyra::DefaultBlockedLinearOp<double>& mat) const;

    /// the composed system matrix
    Teuchos::RCP<Thyra::DefaultBlockedLinearOp<double> > SysMat() const { return mat_; }

  private:

    /// setup structure algorithm
    void SetupStructure();

    /// setup fluid algorithm
    void SetupFluid();

    /// setup ale algorithm
    void SetupAle();

    /// start a new time step
    void PrepareTimeStep();

    /// take current results for converged and save for next time step
    void Update();

    /// write output
    void Output();

    /// setup of NOX convergance tests
    Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Thyra::Group> grp);

    /// communication (mainly for screen output)
    Epetra_Comm& comm_;

    /// structural field solver
    Teuchos::RCP<Structure> structure_;

    /// fluid field solver
    Teuchos::RCP<Fluid>     fluid_;

    /// ale field solver
    Teuchos::RCP<FSI::Ale>  ale_;

    /// coupling of structure and fluid at the interface
    Coupling               coupsf_;

    /// coupling of structure and ale at the interface
    Coupling               coupsa_;

    /// coupling of fluid and ale
    Coupling               coupfa_;

    /// the dof row map of the coupled system
    Teuchos::RCP<Thyra::DefaultProductVectorSpace<double> > dofrowmap_;

    //! @name Wrapped Thyra maps
    /// These are the maps that constitude the composed system

    /// structure field map
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > smap_;

    /// fluid field map
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > fmap_;

    /// ale field map
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > amap_;

    /// structure interface map, included twice in composed system
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > sfmap_;

    //@}

    //! @name global system

    /// identity matrix at structural interface
    Teuchos::RCP<const Thyra::LinearOpBase<double> > sfidentity_;

    /// Lagrangian multiplier for fluid-structure coupling
    Teuchos::RCP< Epetra_Vector > ifstruct_;

    /// Lagrangian multiplier for ale-structure coupling
    Teuchos::RCP< Epetra_Vector > iastruct_;

    /// composed global system matrix
    Teuchos::RCP<Thyra::DefaultBlockedLinearOp<double> > mat_;

    //@}

    //! @name global system solvers

    /// preconditioner factory for global system solver
    Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > preconditionerfactory_;

    /// solver factory for global system solver
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solverfactory_;

    /// solver factory for structure system in preconditioner
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > structsolverfactory_;

    /// solver factory for fluid system in preconditioner
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > fluidsolverfactory_;

    /// solver factory for ale system in preconditioner
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > alesolverfactory_;

    //@}

    /// output utilities
    Teuchos::RefCountPtr<NOX::Utils> utils_;

    /// very simple one processor debug output facility
    Debug debug;

    //! @name Time step valiables
    /// guess what
    int step_;
    int nstep_;
    double time_;
    double maxtime_;
    double dt_;
    //@}
  };
}

#endif
#endif
