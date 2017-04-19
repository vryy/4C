/*!----------------------------------------------------------------------
\file scatra_utils.H
\brief Performs ScaTra specifc functions not yet generalized for other fields.

\level 2

\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15236

*----------------------------------------------------------------------*/

#ifndef SCATRA_UTILS_H
#define SCATRA_UTILS_H

#include <Epetra_MultiVector.h>

#include "../drt_lib/drt_element.H"

/*==========================================================================*/
//! Algorithm and time integration classes for scalar transport problems
/*==========================================================================*/
namespace SCATRA
{

namespace SCATRAUTILS
{

//! Calculate the reconstructed nodal gradient at a node by means of mean value averaging
template <const int dim>
Teuchos::RCP<Epetra_MultiVector> ComputeGradientAtNodesMeanAverage(
    Teuchos::RCP<DRT::Discretization>  discret,
    const Teuchos::RCP<const Epetra_Vector> state,
    const int scatra_dofid
    );

//! Calculate the reconstructed nodal gradient at a node by means of mean value averaging
template <const int dim, DRT::Element::DiscretizationType DISTYPE>
LINALG::Matrix<dim,1> DoMeanValueAveragingOfElementGradientNode(
    Teuchos::RCP<DRT::Discretization>    discret,
    std::vector<const DRT::Element*> elements,
    Teuchos::RCP<Epetra_Vector> phinp_node,
    const int nodegid,
    const int scatra_dofid);

}
} // namespace SCATRA
#endif  // #ifndef SCATRA_UTILS_H
