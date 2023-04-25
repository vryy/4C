/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of free service functions for elements

\level 3
*/
/*----------------------------------------------------------------------*/

#include "so3_element_service.H"
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseSolver.h>
#include <Teuchos_RCPDecl.hpp>
#include <type_traits>
#include "so3_hex8.H"
#include "so3_tet10.H"
#include "linalg_serialdensevector.H"
#include "fem_general_utils_integration.H"
#include "fem_general_utils_gausspoints.H"

#include "fem_general_utils_fem_shapefunctions.H"

#include "so3_hex8.H"
#include "so3_tet10.H"

#include "fem_general_utils_fem_shapefunctions.H"

void DRT::ELEMENTS::AssembleGaussPointValues(
    std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data,
    const Epetra_SerialDenseMatrix& gp_data, const DRT::Element& ele)
{
  for (int gp = 0; gp < gp_data.M(); ++gp)
  {
    const Epetra_BlockMap& elemap = global_data[gp]->Map();
    int lid = elemap.LID(ele.Id());
    if (lid != -1)
    {
      for (int i = 0; i < gp_data.N(); ++i)
      {
        (*((*global_data[gp])(i)))[lid] += gp_data(gp, i);
      }
    }
  }
}

void DRT::ELEMENTS::AssembleNodalElementCount(
    Epetra_IntVector& global_count, const DRT::Element& ele)
{
  for (int n = 0; n < ele.NumNode(); ++n)
  {
    const int lid = global_count.Map().LID(ele.NodeIds()[n]);

    if (lid != -1)
    {
      global_count[lid] += 1;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
std::vector<double> DRT::ELEMENTS::ProjectNodalQuantityToXi(
    const LINALG::Matrix<3, 1>& xi, const std::vector<double>& nodal_quantity)
{
  const int numNodesPerElement = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  LINALG::Matrix<numNodesPerElement, 1> shapefunct(true);
  DRT::UTILS::shape_function<distype>(xi, shapefunct);

  const int num_dof_per_node = static_cast<int>(nodal_quantity.size()) / numNodesPerElement;
  std::vector<double> projected_quantities(num_dof_per_node, 0.0);

  for (int dof = 0; dof < num_dof_per_node; ++dof)
  {
    for (int i = 0; i < numNodesPerElement; ++i)
    {
      projected_quantities[dof] += nodal_quantity[i * num_dof_per_node + dof] * shapefunct(i);
    }
  }

  return projected_quantities;
}

template std::vector<double> DRT::ELEMENTS::ProjectNodalQuantityToXi<DRT::Element::hex8>(
    const LINALG::Matrix<3, 1>&, const std::vector<double>&);
template std::vector<double> DRT::ELEMENTS::ProjectNodalQuantityToXi<DRT::Element::hex27>(
    const LINALG::Matrix<3, 1>&, const std::vector<double>&);
template std::vector<double> DRT::ELEMENTS::ProjectNodalQuantityToXi<DRT::Element::tet4>(
    const LINALG::Matrix<3, 1>&, const std::vector<double>&);
template std::vector<double> DRT::ELEMENTS::ProjectNodalQuantityToXi<DRT::Element::tet10>(
    const LINALG::Matrix<3, 1>&, const std::vector<double>&);
template std::vector<double> DRT::ELEMENTS::ProjectNodalQuantityToXi<DRT::Element::wedge6>(
    const LINALG::Matrix<3, 1>&, const std::vector<double>&);
