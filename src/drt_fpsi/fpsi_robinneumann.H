/*----------------------------------------------------------------------*/
/*!
\file fpsi_robinneumann.H

\brief Solve FSI problems using a Robin-Neumann partitioning approach

<pre>
Maintainer: Andreas Rauch
            rauch.andreas@mytum.de

</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FPSI_ROBINNEUMANN_H
#define FPSI_ROBINNEUMANN_H

#include "fpsi_partitioned.H"

namespace FPSI {


class RobinNeumann : public Partitioned
{
public:

  explicit RobinNeumann(const Epetra_Comm& comm,
                        const Teuchos::ParameterList& fpsidynparams);

protected:

  /// composed FSI operator
  virtual void FSIOp(const Epetra_Vector &x, Epetra_Vector &F);

  /// interface fluid operator
  virtual Teuchos::RCP<Epetra_Vector> FluidOp(Teuchos::RCP<Epetra_Vector> idisp);

  /// interface structural operator
  virtual void PoroOp(Teuchos::RCP<Epetra_Vector> iforce);

  virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

private:

  bool displacementcoupling_;

};

}

#endif // FPSI_ROBINNEUMANN_H
