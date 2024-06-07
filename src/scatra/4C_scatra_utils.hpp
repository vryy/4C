/*----------------------------------------------------------------------*/
/*! \file
\brief Performs ScaTra specifc functions not yet generalized for other fields.

\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_UTILS_HPP
#define FOUR_C_SCATRA_UTILS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"

#include <Epetra_MultiVector.h>

FOUR_C_NAMESPACE_OPEN

namespace ScaTra::ScaTraUtils
{
  //! check, if s2i condition definition is consistent
  void CheckConsistencyOfS2IConditions(Teuchos::RCP<Discret::Discretization> discretization);

  //! check, if nodes of input conditions equal s2i kinetics condition
  void CheckConsistencyWithS2IKineticsCondition(const std::string& condition_to_be_tested,
      Teuchos::RCP<Discret::Discretization> discretization);

  //! Calculate the reconstructed nodal gradient at a node by means of mean value averaging
  template <const int dim>
  Teuchos::RCP<Epetra_MultiVector> ComputeGradientAtNodesMeanAverage(
      Teuchos::RCP<Discret::Discretization> discret, const Teuchos::RCP<const Epetra_Vector> state,
      const int scatra_dofid);

  //! Calculate the reconstructed nodal gradient at a node by means of mean value averaging
  template <const int dim, Core::FE::CellType DISTYPE>
  Core::LinAlg::Matrix<dim, 1> DoMeanValueAveragingOfElementGradientNode(
      Teuchos::RCP<Discret::Discretization> discret,
      std::vector<const Core::Elements::Element*> elements, Teuchos::RCP<Epetra_Vector> phinp_node,
      const int nodegid, const int scatra_dofid);

}  // namespace ScaTra::ScaTraUtils
FOUR_C_NAMESPACE_CLOSE

#endif
