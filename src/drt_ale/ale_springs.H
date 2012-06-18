/*----------------------------------------------------------------------*/
/*!
\file ale_springs.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ALE_SPRINGS_H
#define ALE_SPRINGS_H

#include "ale.H"

/// forward declarations
namespace LINALG
{
  class Solver;
  class Preconditioner;
  class SparseMatrix;
  class BlockSparseMatrixBase;
  class SparseOperator;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace ALE {


/// springs ALE algorithm
/*!
  Incremental springs based implementation.
 */
class AleSprings : public Ale
{

public:

  AleSprings(Teuchos::RCP<DRT::Discretization> actdis,
             Teuchos::RCP<LINALG::Solver> solver,
             Teuchos::RCP<Teuchos::ParameterList> params,
             Teuchos::RCP<IO::DiscretizationWriter> output,
             bool dirichletcond=true);

  /// build system matrix
  virtual void BuildSystemMatrix(bool full=true);

  //! @name Algorithm core methods
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp);
  virtual void Solve();
  virtual void Update();
  virtual void Output();
  //@}

  Teuchos::RCP<LINALG::Preconditioner> ConstPreconditioner() { return Teuchos::null; }

  virtual Teuchos::RCP<Epetra_Vector> ExtractDisplacement() const;

private:

  /// element call
  void EvaluateElements();

  Teuchos::RCP<Epetra_Vector>    incr_;

};

}

#endif
