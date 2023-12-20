/*----------------------------------------------------------------------*/
/*! \file

\brief base XFEM Neumann boundary conditions

\level 2


\warning think about removing these routines!!!

*/
/*----------------------------------------------------------------------*/


#ifndef BACI_XFEM_NEUMANN_H
#define BACI_XFEM_NEUMANN_H


#include "baci_config.H"

#include <Teuchos_RCP.hpp>

class Epetra_Vector;
namespace Teuchos
{
  class ParameterList;
}

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Condition;
  class Discretization;
  class Element;
}  // namespace DRT

namespace CORE::LINALG
{
  class SparseOperator;
}

namespace XFEM
{
  /// evaluate Neumann boundary conditions
  void EvaluateNeumann(Teuchos::ParameterList& params, Teuchos::RCP<DRT::Discretization> discret,
      Teuchos::RCP<Epetra_Vector> systemvector,
      Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix = Teuchos::null);

  /// evaluate Neumann boundary conditions
  void EvaluateNeumann(Teuchos::ParameterList& params, Teuchos::RCP<DRT::Discretization> discret,
      Epetra_Vector& systemvector, CORE::LINALG::SparseOperator* systemmatrix = nullptr);

  /// evaluate standard Neumann boundary conditions
  void EvaluateNeumannStandard(std::multimap<std::string, DRT::Condition*>& condition,
      const double time, bool assemblemat, Teuchos::ParameterList& params,
      Teuchos::RCP<DRT::Discretization> discret, Epetra_Vector& systemvector,
      CORE::LINALG::SparseOperator* systemmatrix);


}  // namespace XFEM

BACI_NAMESPACE_CLOSE

#endif  // XFEM_NEUMANN_H
