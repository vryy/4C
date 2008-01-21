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
namespace MFSI
{

  /// base class for all monolithic FSI algorithms
  /*!

  A common base class for all monolithic FSI algorithms. The idea is
  to keep all algorithms as close together as possible. This way it is
  easier to change and understand them. Hopefully.

  \author u.kue
  \date 12/07
  */
  class Algorithm : public Thyra::StateFuncModelEvaluatorBase<double>
  {
  public:
    Algorithm(Epetra_Comm& comm);

    /// The time loop. Start this. Everything is in here.
    virtual void Timeloop();

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

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Thyra::DefaultProductVector<double> &f) const = 0;

    /// setup composed system matrix from field solvers
    virtual void SetupSysMat(Thyra::DefaultBlockedLinearOp<double>& mat) const = 0;

    //@}

    /// the composed system matrix
    Teuchos::RCP<Thyra::DefaultBlockedLinearOp<double> > SysMat() const { return mat_; }

    /// Evaluate all fields at position x
    virtual void Evaluate(Teuchos::RCP<const Thyra::DefaultProductVector<double> > x) const;

    /// Extract initial guess from fields
    virtual void InitialGuess(Thyra::DefaultProductVector<double>& ig) = 0;

  protected:

    /// translate baci solver variables to Thyra
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
    CreateSolverFactory(struct _SOLVAR* actsolv);

    /// setup of NOX convergance tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Thyra::Group> grp) = 0;

    /// extract the three field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param sx (o) structural displacements
      \param fx (o) fluid velocities and pressure
      \param ax (o) ale displacements
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Thyra::DefaultProductVector<double> > x,
                                     Teuchos::RCP<const Epetra_Vector>& sx,
                                     Teuchos::RCP<const Epetra_Vector>& fx,
                                     Teuchos::RCP<const Epetra_Vector>& ax) const = 0;

    //! @name Access methods for subclasses

    Epetra_Comm& Comm() { return comm_; }
    double Dt() const { return dt_; }
    int Step() const { return step_; }

    Teuchos::RefCountPtr<NOX::Utils> Utils() const { return utils_; }

    Teuchos::RCP<Structure> StructureField() const { return structure_; }
    Teuchos::RCP<Fluid>     FluidField()     const { return fluid_; }
    Teuchos::RCP<FSI::Ale>  AleField()       const { return ale_; }

    Teuchos::RCP<Thyra::DefaultProductVectorSpace<double> > DofRowMap() const { return dofrowmap_; }
    void SetDofRowMap(Teuchos::RCP<Thyra::DefaultProductVectorSpace<double> > drm) { dofrowmap_ = drm; }

    Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > PreconditionerFactory() const { return preconditionerfactory_; }
    void SetPreconditionerFactory(Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > pf) { preconditionerfactory_ = pf; }

    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > SolverFactory() const { return solverfactory_; }
    void SetSolverFactory(Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > sf) { solverfactory_ = sf; }

    //@}

    //! @name Time loop building blocks

    /// tests if there are more time steps to do
    bool NotFinished() { return step_ < nstep_ and time_ <= maxtime_; }

    /// start a new time step
    void PrepareTimeStep();

    /// take current results for converged and save for next time step
    void Update();

    /// write output
    void Output();

    //@}

    //! @name Transfer helpers

    Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<Epetra_Vector> iv) const;
    Teuchos::RCP<Epetra_Vector> AleToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
    Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const;
    Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
    Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const;

    Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const;
    Teuchos::RCP<Epetra_Vector> AleToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;
    Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;
    Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;
    Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

    Coupling& StructureFluidCoupling() { return coupsf_; }
    Coupling& StructureAleCoupling()   { return coupsa_; }
    Coupling& FluidAleCoupling()       { return coupfa_; }

    const Coupling& StructureFluidCoupling() const { return coupsf_; }
    const Coupling& StructureAleCoupling()   const { return coupsa_; }
    const Coupling& FluidAleCoupling()       const { return coupfa_; }

  private:

    /// setup structure algorithm
    void SetupStructure();

    /// setup fluid algorithm
    void SetupFluid();

    /// setup ale algorithm
    void SetupAle();

    /// communication (mainly for screen output)
    Epetra_Comm& comm_;

    /// coupling of structure and fluid at the interface
    Coupling               coupsf_;

    /// coupling of structure and ale at the interface
    Coupling               coupsa_;

    /// coupling of fluid and ale
    Coupling               coupfa_;

    //! @name global system

    /// the dof row map of the coupled system
    Teuchos::RCP<Thyra::DefaultProductVectorSpace<double> > dofrowmap_;

    /// composed global system matrix
    Teuchos::RCP<Thyra::DefaultBlockedLinearOp<double> > mat_;

    //@}

    //! @name global system solvers

    /// preconditioner factory for global system solver
    Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > preconditionerfactory_;

    /// solver factory for global system solver
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solverfactory_;

    //@}

    /// structural field solver
    Teuchos::RCP<Structure> structure_;

    /// fluid field solver
    Teuchos::RCP<Fluid>     fluid_;

    /// ale field solver
    Teuchos::RCP<FSI::Ale>  ale_;

    /// output utilities
    Teuchos::RefCountPtr<NOX::Utils> utils_;

    //! @name Time step valiables
    /// guess what
    int step_;
    int nstep_;
    double time_;
    double maxtime_;
    double dt_;
    int nlnstep_;
    //@}

    //! @name time measurement

    Teuchos::RCP<Teuchos::Time> evaluatetimer_;

    //@}

  protected:

    /// very simple one processor debug output facility
    Debug debug_;
  };
}

#endif
#endif
