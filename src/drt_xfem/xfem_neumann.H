/*----------------------------------------------------------------------*/
/*!
\file xfem_neumann.H

\brief base XFEM Neumann boundary conditions

\level 2

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>

\warning think about removing these routines!!!

*/
/*----------------------------------------------------------------------*/


#ifndef XFEM_NEUMANN_H_
#define XFEM_NEUMANN_H_


#include <Teuchos_RCP.hpp>


class Epetra_Vector;

namespace DRT
{
  class Condition;
  class Discretization;
  class Element;
}  // namespace DRT

namespace LINALG
{
  class SparseOperator;
}
namespace Teuchos
{
  class ParameterList;
}


namespace XFEM
{
  /// evaluate Neumann boundary conditions
  void EvaluateNeumann(Teuchos::ParameterList& params, Teuchos::RCP<DRT::Discretization> discret,
      Teuchos::RCP<Epetra_Vector> systemvector,
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix = Teuchos::null);

  /// evaluate Neumann boundary conditions
  void EvaluateNeumann(Teuchos::ParameterList& params, Teuchos::RCP<DRT::Discretization> discret,
      Epetra_Vector& systemvector, LINALG::SparseOperator* systemmatrix = NULL);

  /// evaluate standard Neumann boundary conditions
  void EvaluateNeumannStandard(std::multimap<std::string, DRT::Condition*>& condition,
      const double time, bool assemblemat, Teuchos::ParameterList& params,
      Teuchos::RCP<DRT::Discretization> discret, Epetra_Vector& systemvector,
      LINALG::SparseOperator* systemmatrix);


}  // namespace XFEM

#endif /* XFEM_NEUMANN_H_ */
