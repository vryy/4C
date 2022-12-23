/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of free service functions for elements

\level 3
*/
/*----------------------------------------------------------------------*/

#include "so_element_service.H"
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseSolver.h>
#include <Teuchos_RCPDecl.hpp>
#include <type_traits>
#include "so_hex8.H"
#include "so_tet10.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

#include "utils_fem_shapefunctions.H"

namespace
{
  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DRT::is_tet_v<distype>, bool> = true>
  inline DRT::Element::DiscretizationType GetGaussPointExtrapolationBaseDistype(unsigned numgp)
  {
    if (numgp < 4) return DRT::Element::DiscretizationType::point1;
    if (numgp < 10) return DRT::Element::DiscretizationType::tet4;
    return DRT::Element::DiscretizationType::tet10;
  }

  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DRT::is_hex_v<distype>, bool> = true>
  inline DRT::Element::DiscretizationType GetGaussPointExtrapolationBaseDistype(unsigned numgp)
  {
    if (numgp < 8) return DRT::Element::DiscretizationType::point1;
    if (numgp < 20) return DRT::Element::DiscretizationType::hex8;
    if (numgp < 27) return DRT::Element::DiscretizationType::hex20;
    return DRT::Element::DiscretizationType::hex27;
  }

  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DRT::is_quad_v<distype>, bool> = true>
  inline DRT::Element::DiscretizationType GetGaussPointExtrapolationBaseDistype(unsigned numgp)
  {
    if (numgp < 4) return DRT::Element::DiscretizationType::point1;
    if (numgp < 8) return DRT::Element::DiscretizationType::quad4;
    if (numgp < 9) return DRT::Element::DiscretizationType::quad8;
    return DRT::Element::DiscretizationType::quad9;
  }

  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DRT::is_tri_v<distype>, bool> = true>
  inline DRT::Element::DiscretizationType GetGaussPointExtrapolationBaseDistype(unsigned numgp)
  {
    if (numgp < 3) return DRT::Element::DiscretizationType::point1;
    if (numgp < 6) return DRT::Element::DiscretizationType::tri3;
    return DRT::Element::DiscretizationType::tri6;
  }

  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DRT::is_wedge_v<distype>, bool> = true>
  inline DRT::Element::DiscretizationType GetGaussPointExtrapolationBaseDistype(unsigned numgp)
  {
    if (numgp < 6) return DRT::Element::DiscretizationType::point1;
    if (numgp < 15) return DRT::Element::DiscretizationType::wedge6;
    return DRT::Element::DiscretizationType::wedge15;
  }

  template <DRT::Element::DiscretizationType distype,
      std::enable_if_t<DRT::is_pyramid_v<distype>, bool> = true>
  inline DRT::Element::DiscretizationType GetGaussPointExtrapolationBaseDistype(unsigned numgp)
  {
    if (numgp < 5) return DRT::Element::DiscretizationType::point1;
    return DRT::Element::DiscretizationType::pyramid5;
  }

  template <DRT::Element::DiscretizationType distype>
  LINALG::SerialDenseMatrix EvaluateBaseShapeFunctionsAtGaussPoints(
      const DRT::Element::DiscretizationType base_distype,
      const DRT::UTILS::GaussIntegration& intpoints)
  {
    constexpr int nsd = DRT::UTILS::DisTypeToDim<distype>::dim;
    int base_numnod = DRT::UTILS::getNumberOfElementNodes(base_distype);
    LINALG::SerialDenseMatrix mat(intpoints.NumPoints(), base_numnod);


    for (int gp = 0; gp < intpoints.NumPoints(); ++gp)
    {
      LINALG::Matrix<nsd, 1> xi(intpoints.Point(gp), true);

      LINALG::SerialDenseVector shape_functions(base_numnod);
      DRT::UTILS::shape_function_dim<LINALG::Matrix<nsd, 1>, LINALG::SerialDenseVector, nsd>(
          xi, shape_functions, base_distype);

      for (int inode = 0; inode < base_numnod; ++inode)
      {
        mat(gp, inode) = shape_functions(inode);
      }
    }
    return mat;
  }

  LINALG::SerialDenseMatrix EvaluateProjectionGaussPointsToBaseDistype(
      const LINALG::SerialDenseMatrix& shapefcns_at_gps)
  {
    LINALG::SerialDenseMatrix shapefunctions_at_gps_copy(shapefcns_at_gps);
    if (shapefcns_at_gps.M() == shapefcns_at_gps.N())
    {
      // Todo: simply invert the matrix
      Epetra_SerialDenseSolver matrixInverter;
      matrixInverter.SetMatrix(shapefunctions_at_gps_copy);
      matrixInverter.Invert();

      if (!matrixInverter.Inverted())
      {
        dserror(
            "Failed to invert the matrix of the shapefunctions evaluated at the Gauss points. It "
            "looks like this element does not support the default way to extrapolate quantities "
            "from Gauss points to nodes.");
      }

      dsassert(shapefunctions_at_gps_copy.A() == matrixInverter.AF(),
          "Inverse of the matrix was not computed in place, but we expect that. Unfortunately, the "
          "Trilinos documentation is ambiguous here.");

      return shapefunctions_at_gps_copy;
    }

    // solve least square algorithm
    LINALG::SerialDenseMatrix matTmat(shapefcns_at_gps.N(), shapefcns_at_gps.N());
    matTmat.Multiply('T', 'N', 1.0, shapefcns_at_gps, shapefcns_at_gps, 0.0);

    {
      Epetra_SerialDenseSolver matrixInverter;
      matrixInverter.SetMatrix(matTmat);
      matrixInverter.Invert();

      if (!matrixInverter.Inverted())
      {
        dserror(
            "Failed to invert the matrix of the shapefunctions evaluated at the Gauss points. It "
            "looks like this element does not support the default way to extrapolate quantities "
            "from Gauss points to nodes.");
      }

      dsassert(matTmat.A() == matrixInverter.AF(),
          "Inverse of the matrix was not computed in place, but we expect that. Unfortunately, the "
          "Trilinos documentation is ambiguous here.");
    }


    LINALG::SerialDenseMatrix matrix_gp_to_base(shapefcns_at_gps.M(), shapefcns_at_gps.N());
    matrix_gp_to_base.Multiply('N', 'T', 1.0, matTmat, shapefunctions_at_gps_copy, 0.0);

    return matrix_gp_to_base;
  }

  template <DRT::Element::DiscretizationType distype>
  LINALG::SerialDenseMatrix EvaluateProjectionGaussPointsToDistype(
      const LINALG::SerialDenseMatrix& matrix_gp_to_base,
      DRT::Element::DiscretizationType base_distype)
  {
    if (base_distype == distype)
    {
      return matrix_gp_to_base;
    }
    constexpr int nsd = DRT::UTILS::DisTypeToDim<distype>::dim;
    int base_numnod = DRT::UTILS::getNumberOfElementNodes(base_distype);

    LINALG::SerialDenseMatrix matrix_base_to_dis(
        DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, base_numnod);

    for (int dis_inode = 0;
         dis_inode < DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement; ++dis_inode)
    {
      LINALG::SerialDenseMatrix reference_nodes =
          DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

      if (nsd == 3)
      {
        LINALG::SerialDenseVector shape_functions(base_numnod);
        DRT::UTILS::shape_function_3D(shape_functions, reference_nodes(0, dis_inode),
            reference_nodes(1, dis_inode), reference_nodes(2, dis_inode), base_distype);

        for (int basedis_inode = 0; basedis_inode < base_numnod; ++basedis_inode)
        {
          matrix_base_to_dis(dis_inode, basedis_inode) = shape_functions(basedis_inode);
        }
      }
    }

    LINALG::SerialDenseMatrix matrix_gp_to_nodes(
        DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, matrix_gp_to_base.N());

    // extend matrix from base_distype to distype
    matrix_gp_to_nodes.Multiply('N', 'N', 1.0, matrix_base_to_dis, matrix_gp_to_base, 0.0);

    return matrix_gp_to_nodes;
  }

}  // namespace

template <class T>
void DRT::ELEMENTS::AssembleExtrapolatedNodalValues(Epetra_MultiVector& global_data,
    const T& nodal_data, const DRT::Element& ele, bool nodal_average)
{
  for (decltype(nodal_data.M()) i = 0; i < nodal_data.M(); ++i)
  {
    const int lid = global_data.Map().LID(ele.NodeIds()[i]);
    if (lid >= 0)  // rownode
    {
      const double invmyadjele = (nodal_average) ? 1.0 / ele.Nodes()[i]->NumElement() : 1.0;
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

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ExtrapolateGPQuantityToNodesAndAssemble(const DRT::Element& ele,
    const LINALG::SerialDenseMatrix& gp_data, Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::GaussIntegration& integration)
{
  LINALG::SerialDenseMatrix nodal_quantity(
      DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, gp_data.N());
  nodal_quantity.Multiply('N', 'N', 1.0,
      EvaluateGaussPointsToNodesExtrapolationMatrix<distype>(integration), gp_data, 0.0);

  AssembleExtrapolatedNodalValues(global_data, nodal_quantity, ele, nodal_average);
}

template <DRT::Element::DiscretizationType distype>
LINALG::SerialDenseMatrix DRT::ELEMENTS::EvaluateGaussPointsToNodesExtrapolationMatrix(
    const DRT::UTILS::GaussIntegration& intpoints)
{
  // TODO: This has to be done only once per simulation (or theoretically during compile time for
  // every distype)
  DRT::Element::DiscretizationType base_distype =
      GetGaussPointExtrapolationBaseDistype<distype>(intpoints.NumPoints());

  dsassert(DRT::UTILS::getNumberOfElementNodes(base_distype) <= intpoints.NumPoints(),
      "The base discretization has more nodes than Gauß points. The extrapolation is not unique! "
      "This should not happen. The evaluation of the base extrapolation type for the number of "
      "gauss points is not correct.");

  LINALG::SerialDenseMatrix shapefcns_at_gps =
      EvaluateBaseShapeFunctionsAtGaussPoints<distype>(base_distype, intpoints);

  LINALG::SerialDenseMatrix matrix_gp_to_base =
      EvaluateProjectionGaussPointsToBaseDistype(shapefcns_at_gps);

  return EvaluateProjectionGaussPointsToDistype<distype>(matrix_gp_to_base, base_distype);
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
    const LINALG::Matrix<NUMNOD_SOH8, MAT::NUM_STRESS_3D>&, const DRT::Element&, bool);
template void DRT::ELEMENTS::AssembleExtrapolatedNodalValues(Epetra_MultiVector&,
    const LINALG::Matrix<NUMNOD_SOTET10, MAT::NUM_STRESS_3D>&, const DRT::Element&, bool);
template void DRT::ELEMENTS::AssembleExtrapolatedNodalValues(
    Epetra_MultiVector&, const LINALG::SerialDenseMatrix&, const DRT::Element&, bool);

template void DRT::ELEMENTS::AssembleAveragedElementValues(Epetra_MultiVector&,
    const LINALG::Matrix<NUMNOD_SOH8, MAT::NUM_STRESS_3D>&, const DRT::Element&);
template void DRT::ELEMENTS::AssembleAveragedElementValues(
    Epetra_MultiVector&, const LINALG::SerialDenseMatrix&, const DRT::Element&);

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

// template DRT::Element::DiscretizationType
// GetGaussPointExtrapolationBaseDistype<DRT::Element::hex8>(
//     unsigned);

template void DRT::ELEMENTS::ExtrapolateGPQuantityToNodesAndAssemble<DRT::Element::hex8>(
    const DRT::Element&, const LINALG::SerialDenseMatrix&, Epetra_MultiVector&, bool,
    const DRT::UTILS::GaussIntegration&);
template void DRT::ELEMENTS::ExtrapolateGPQuantityToNodesAndAssemble<DRT::Element::hex18>(
    const DRT::Element&, const LINALG::SerialDenseMatrix&, Epetra_MultiVector&, bool,
    const DRT::UTILS::GaussIntegration&);
template void DRT::ELEMENTS::ExtrapolateGPQuantityToNodesAndAssemble<DRT::Element::hex20>(
    const DRT::Element&, const LINALG::SerialDenseMatrix&, Epetra_MultiVector&, bool,
    const DRT::UTILS::GaussIntegration&);
template void DRT::ELEMENTS::ExtrapolateGPQuantityToNodesAndAssemble<DRT::Element::hex27>(
    const DRT::Element&, const LINALG::SerialDenseMatrix&, Epetra_MultiVector&, bool,
    const DRT::UTILS::GaussIntegration&);
template void DRT::ELEMENTS::ExtrapolateGPQuantityToNodesAndAssemble<DRT::Element::tet4>(
    const DRT::Element&, const LINALG::SerialDenseMatrix&, Epetra_MultiVector&, bool,
    const DRT::UTILS::GaussIntegration&);
template void DRT::ELEMENTS::ExtrapolateGPQuantityToNodesAndAssemble<DRT::Element::tet10>(
    const DRT::Element&, const LINALG::SerialDenseMatrix&, Epetra_MultiVector&, bool,
    const DRT::UTILS::GaussIntegration&);
template void DRT::ELEMENTS::ExtrapolateGPQuantityToNodesAndAssemble<DRT::Element::wedge6>(
    const DRT::Element&, const LINALG::SerialDenseMatrix&, Epetra_MultiVector&, bool,
    const DRT::UTILS::GaussIntegration&);
template void DRT::ELEMENTS::ExtrapolateGPQuantityToNodesAndAssemble<DRT::Element::pyramid5>(
    const DRT::Element&, const LINALG::SerialDenseMatrix&, Epetra_MultiVector&, bool,
    const DRT::UTILS::GaussIntegration&);
