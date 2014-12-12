/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_s2i.H

\brief Scatra-scatra interface coupling strategy for standard scalar transport problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_H
#define SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_H

#include "scatra_timint_meshtying_strategy_base.H"

// forward declarations
namespace ADAPTER
{
  class Coupling;
}

namespace FSI
{
  namespace UTILS
  {
    class MatrixColTransform;
  }
}

namespace LINALG
{
  class MultiMapExtractor;
  class MapExtractor;
  class BlockSparseMatrixBase;
}

namespace SCATRA
{

/*!
\brief Scatra-scatra interface coupling strategy for standard scalar transport problems

To keep the scalar transport time integrator class and derived classes as plain as possible,
several algorithmic parts have been encapsulated within separate meshtying strategy classes.
These algorithmic parts include initializing the system matrix and other relevant objects,
computing meshtying residual terms and their linearizations, and solving the resulting
linear system of equations. By introducing a hierarchy of strategies for these algorithmic
parts, a bunch of unhandy if-else selections within the time integrator classes themselves
can be circumvented. This class contains the scatra-scatra interface coupling strategy for
standard scalar transport problems.

\author fang (fang@lnm.mw.tum.de)
*/

class MeshtyingStrategyS2I : public MeshtyingStrategyBase
{
  public:

    //! constructor
    MeshtyingStrategyS2I(
        SCATRA::ScaTraTimIntImpl* scatratimint
        );

    //! destructor
    virtual ~MeshtyingStrategyS2I() {};

    //! compute meshtying residual terms and their linearizations
    void EvaluateMeshtying() const;

    //! initialize meshtying objects
    void InitMeshtying();

    //! initialize system matrix
    Teuchos::RCP<LINALG::SparseOperator> InitSystemMatrix() const;

    //! solve resulting linear system of equations
    void Solve(
        const Teuchos::RCP<LINALG::Solver>&            solver,         //! solver
        const Teuchos::RCP<LINALG::SparseOperator>&    systemmatrix,   //! system matrix
        const Teuchos::RCP<Epetra_Vector>&             increment,      //! increment vector
        const Teuchos::RCP<Epetra_Vector>&             residual,       //! residual vector
        const Teuchos::RCP<Epetra_Vector>&             phinp,          //! state vector at time n+1
        const int&                                     iteration,      //! number of current Newton-Raphson iteration
        const Teuchos::RCP<LINALG::KrylovProjector>&   projector       //! Krylov projector
        ) const;

  private:

    //! copy constructor
    MeshtyingStrategyS2I(const MeshtyingStrategyS2I& old);

    //! global map extractor (0: other, 1: slave, 2: master)
    Teuchos::RCP<LINALG::MultiMapExtractor> maps_;

    //! interface map extractor (0: master, 1: slave)
    Teuchos::RCP<LINALG::MapExtractor> imaps_;

    //! interface coupling adapter
    Teuchos::RCP<ADAPTER::Coupling> icoup_;

    //! auxiliary system matrix for coupling terms
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> auxmat_;

    //! transformed state variables at scatra-scatra interface
    Teuchos::RCP<Epetra_Vector> iphinp_;

    //! transformation operators for auxiliary system matrix
    Teuchos::RCP<FSI::UTILS::MatrixColTransform> imastertoslavetransform_;
    Teuchos::RCP<FSI::UTILS::MatrixColTransform> islavetomastertransform_;

}; // class MeshtyingStrategyS2I
} // namespace SCATRA
#endif // #ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_H
