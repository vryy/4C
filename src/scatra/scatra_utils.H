/*----------------------------------------------------------------------*/
/*! \file
\brief Performs ScaTra specifc functions not yet generalized for other fields.

\level 2


*----------------------------------------------------------------------*/

#ifndef SCATRA_UTILS_H
#define SCATRA_UTILS_H

#include <Epetra_MultiVector.h>

#include "lib_element.H"

namespace SCATRA
{
  namespace SCATRAUTILS
  {
    //! check, if nodes of input conditions equal s2i kinetics condition
    void CheckConsistencyWithS2IKineticsCondition(const std::string& condition_to_be_tested,
        Teuchos::RCP<DRT::Discretization> discretization);

    //! Calculate the reconstructed nodal gradient at a node by means of mean value averaging
    template <const int dim>
    Teuchos::RCP<Epetra_MultiVector> ComputeGradientAtNodesMeanAverage(
        Teuchos::RCP<DRT::Discretization> discret, const Teuchos::RCP<const Epetra_Vector> state,
        const int scatra_dofid);

    //! Calculate the reconstructed nodal gradient at a node by means of mean value averaging
    template <const int dim, DRT::Element::DiscretizationType DISTYPE>
    LINALG::Matrix<dim, 1> DoMeanValueAveragingOfElementGradientNode(
        Teuchos::RCP<DRT::Discretization> discret, std::vector<const DRT::Element*> elements,
        Teuchos::RCP<Epetra_Vector> phinp_node, const int nodegid, const int scatra_dofid);

  }  // namespace SCATRAUTILS
}  // namespace SCATRA
#endif  // SCATRA_UTILS_H
