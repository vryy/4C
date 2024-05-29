/*----------------------------------------------------------------------*/
/*! \file

\brief base XFEM Neumann boundary conditions

\level 2


\warning think about removing these routines!!!

*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_XFEM_NEUMANN_HPP
#define FOUR_C_XFEM_NEUMANN_HPP


#include "4C_config.hpp"

#include "4C_discretization_condition.hpp"

#include <Teuchos_RCP.hpp>

class Epetra_Vector;
namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}  // namespace DRT

namespace CORE::Elements
{
  class Element;
}

namespace CORE::LINALG
{
  class SparseOperator;
}

namespace XFEM
{
  /// evaluate Neumann boundary conditions
  void evaluate_neumann(Teuchos::ParameterList& params, Teuchos::RCP<DRT::Discretization> discret,
      Teuchos::RCP<Epetra_Vector> systemvector,
      Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix = Teuchos::null);

  /// evaluate Neumann boundary conditions
  void evaluate_neumann(Teuchos::ParameterList& params, Teuchos::RCP<DRT::Discretization> discret,
      Epetra_Vector& systemvector, CORE::LINALG::SparseOperator* systemmatrix = nullptr);

  /// evaluate standard Neumann boundary conditions
  void EvaluateNeumannStandard(std::multimap<std::string, CORE::Conditions::Condition*>& condition,
      const double time, bool assemblemat, Teuchos::ParameterList& params,
      Teuchos::RCP<DRT::Discretization> discret, Epetra_Vector& systemvector,
      CORE::LINALG::SparseOperator* systemmatrix);


}  // namespace XFEM

FOUR_C_NAMESPACE_CLOSE

#endif
