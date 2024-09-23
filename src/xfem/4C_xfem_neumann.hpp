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

#include "4C_fem_condition.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_RCP.hpp>

class Epetra_Vector;

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace Core::LinAlg
{
  class SparseOperator;
}

namespace XFEM
{
  /// evaluate Neumann boundary conditions
  void evaluate_neumann(Teuchos::ParameterList& params,
      Teuchos::RCP<Core::FE::Discretization> discret, Teuchos::RCP<Epetra_Vector> systemvector,
      Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix = Teuchos::null);

  /// evaluate Neumann boundary conditions
  void evaluate_neumann(Teuchos::ParameterList& params,
      Teuchos::RCP<Core::FE::Discretization> discret, Epetra_Vector& systemvector,
      Core::LinAlg::SparseOperator* systemmatrix = nullptr);

  /// evaluate standard Neumann boundary conditions
  void evaluate_neumann_standard(
      std::multimap<std::string, Core::Conditions::Condition*>& condition, const double time,
      bool assemblemat, Teuchos::ParameterList& params,
      Teuchos::RCP<Core::FE::Discretization> discret, Epetra_Vector& systemvector,
      Core::LinAlg::SparseOperator* systemmatrix);


}  // namespace XFEM

FOUR_C_NAMESPACE_CLOSE

#endif
