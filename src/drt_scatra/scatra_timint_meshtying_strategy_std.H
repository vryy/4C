/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_std.H

\brief Standard solution strategy for standard scalar transport problems (without meshtying)

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_STD_H
#define SCATRA_TIMINT_MESHTYING_STRATEGY_STD_H

#include "scatra_timint_meshtying_strategy_base.H"

namespace SCATRA
{

/*!
\brief Standard solution strategy for standard scalar transport problems (without meshtying)

To keep the scalar transport time integrator class and derived classes as plain as possible,
several algorithmic parts have been encapsulated within separate meshtying strategy classes.
These algorithmic parts include initializing the system matrix and other relevant objects,
computing meshtying residual terms and their linearizations, and solving the resulting
linear system of equations. By introducing a hierarchy of strategies for these algorithmic
parts, a bunch of unhandy if-else selections within the time integrator classes themselves
can be circumvented. This class contains the standard solution strategy for standard scalar
transport problems without meshtying.

\author fang (fang@lnm.mw.tum.de)
*/

class MeshtyingStrategyStd : public MeshtyingStrategyBase
{
  public:

    //! constructor
    explicit MeshtyingStrategyStd(
        SCATRA::ScaTraTimIntImpl* scatratimint
        );

    //! destructor
    virtual ~MeshtyingStrategyStd(){return;};

    //! compute meshtying residual terms and their linearizations
    void EvaluateMeshtying();

    //! initialize meshtying objects
    void InitMeshtying();

    //! initialize system matrix
    virtual Teuchos::RCP<LINALG::SparseOperator> InitSystemMatrix() const;

    //! solve resulting linear system of equations
    void Solve(
        const Teuchos::RCP<LINALG::Solver>&            solver,         //!< solver
        const Teuchos::RCP<LINALG::SparseOperator>&    systemmatrix,   //!< system matrix
        const Teuchos::RCP<Epetra_Vector>&             increment,      //!< increment vector
        const Teuchos::RCP<Epetra_Vector>&             residual,       //!< residual vector
        const Teuchos::RCP<Epetra_Vector>&             phinp,          //!< state vector at time n+1
        const int&                                     iteration,      //!< number of current Newton-Raphson iteration
        const Teuchos::RCP<LINALG::KrylovProjector>&   projector       //!< Krylov projector
        ) const;

  protected:

    //! instantiate strategy for Newton-Raphson convergence check
    virtual void InitConvCheckStrategy();

  private:

    //! copy constructor
    MeshtyingStrategyStd(const MeshtyingStrategyStd& old);
}; // class MeshtyingStrategyStd
} // namespace SCATRA
#endif // #ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_STD_H
