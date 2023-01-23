/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of free service functions for elements

\level 3
*/
/*----------------------------------------------------------------------*/

#include "so3_element_service.H"
#include <Epetra_SerialDenseMatrix.h>
#include "so3_hex8.H"
#include "so3_tet10.H"

#include "fem_general_utils_fem_shapefunctions.H"

template <class T>
void DRT::ELEMENTS::AssembleExtrapolatedNodalValues(Epetra_MultiVector& global_data,
    const T& nodal_data, const DRT::ELEMENTS::So_base* ele, bool nodal_average)
{
  for (decltype(nodal_data.M()) i = 0; i < nodal_data.M(); ++i)
  {
    const int lid = global_data.Map().LID(ele->NodeIds()[i]);
    if (lid >= 0)  // rownode
    {
      const double invmyadjele = (nodal_average) ? 1.0 / ele->Nodes()[i]->NumElement() : 1.0;
      for (decltype(nodal_data.N()) j = 0; j < nodal_data.N(); ++j)
      {
        (*(global_data(j)))[lid] += nodal_data(i, j) * invmyadjele;
      }
    }
  }
}

template <class T>
void DRT::ELEMENTS::AssembleAveragedElementValues(
    Epetra_MultiVector& global_data, const T& gp_data, const DRT::Element& ele)
{
  const Epetra_BlockMap& elemap = global_data.Map();
  int lid = elemap.LID(ele.Id());
  if (lid != -1)
  {
    for (decltype(gp_data.N()) i = 0; i < gp_data.N(); ++i)
    {
      double& s = (*(global_data(i)))[lid];  // resolve pointer for faster access
      s = 0.;
      for (decltype(gp_data.M()) j = 0; j < gp_data.M(); ++j)
      {
        s += gp_data(j, i);
      }
      s *= 1.0 / gp_data.M();
    }
  }
}

void DRT::ELEMENTS::AssembleGaussPointValues(
    std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data,
    const Epetra_SerialDenseMatrix& gp_data, const DRT::ELEMENTS::So_base* ele)
{
  for (int gp = 0; gp < gp_data.M(); ++gp)
  {
    const Epetra_BlockMap& elemap = global_data[gp]->Map();
    int lid = elemap.LID(ele->Id());
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
    Epetra_IntVector& global_count, const DRT::ELEMENTS::So_base* ele)
{
  for (int n = 0; n < ele->NumNode(); ++n)
  {
    const int lid = global_count.Map().LID(ele->NodeIds()[n]);

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


// explicit template instantiations
template void DRT::ELEMENTS::AssembleExtrapolatedNodalValues(Epetra_MultiVector&,
    const LINALG::Matrix<NUMNOD_SOH8, MAT::NUM_STRESS_3D>&, const DRT::ELEMENTS::So_base*, bool);
template void DRT::ELEMENTS::AssembleExtrapolatedNodalValues(Epetra_MultiVector&,
    const LINALG::Matrix<NUMNOD_SOTET10, MAT::NUM_STRESS_3D>&, const DRT::ELEMENTS::So_base*, bool);
template void DRT::ELEMENTS::AssembleExtrapolatedNodalValues(
    Epetra_MultiVector&, const LINALG::SerialDenseMatrix&, const DRT::ELEMENTS::So_base*, bool);

template void DRT::ELEMENTS::AssembleAveragedElementValues(Epetra_MultiVector&,
    const LINALG::Matrix<NUMNOD_SOH8, MAT::NUM_STRESS_3D>&, const DRT::Element&);
template void DRT::ELEMENTS::AssembleAveragedElementValues(
    Epetra_MultiVector&, const LINALG::SerialDenseMatrix&, const DRT::Element&);
template void DRT::ELEMENTS::AssembleAveragedElementValues(
    Epetra_MultiVector&, const Epetra_SerialDenseMatrix&, const DRT::Element&);

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