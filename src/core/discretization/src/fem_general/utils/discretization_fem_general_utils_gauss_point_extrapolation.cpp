/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of utility functions for Gauss point data extrapolation

\level 2

*----------------------------------------------------------------------*/
#include "discretization_fem_general_utils_gauss_point_extrapolation.H"
#include <Epetra_SerialDenseSolver.h>
#include "linalg_serialdensevector.H"
#include "discretization_fem_general_utils_gausspoints.H"
#include "discretization_fem_general_utils_integration.H"
#include "discretization_fem_general_utils_local_connectivity_matrices.H"
#include "discretization_fem_general_utils_fem_shapefunctions.H"
#include "lib_node.H"

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

  template <DRT::Element::DiscretizationType distype, class GaussIntegration>
  LINALG::SerialDenseMatrix EvaluateBaseShapeFunctionsAtGaussPoints(
      const DRT::Element::DiscretizationType base_distype, const GaussIntegration& intpoints)
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
      Epetra_SerialDenseSolver matrixInverter;
      matrixInverter.SetMatrix(shapefunctions_at_gps_copy);
      int error_code = matrixInverter.Invert();

      if (error_code != 0)
      {
        dserror(
            "Failed to invert the matrix of the shapefunctions evaluated at the Gauss points. It "
            "looks like this element does not support the default way to extrapolate quantities "
            "from Gauss points to nodes. Error code: %d",
            error_code);
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
      int error_code = matrixInverter.Invert();

      if (error_code != 0)
      {
        dserror(
            "Failed to invert the matrix of the shapefunctions evaluated at the Gauss points. It "
            "looks like this element does not support the default way to extrapolate quantities "
            "from Gauss points to nodes. Error code %d",
            error_code);
      }

      dsassert(matTmat.A() == matrixInverter.AF(),
          "Inverse of the matrix was not computed in place, but we expect that. Unfortunately, the "
          "Trilinos documentation is ambiguous here.");
    }

    LINALG::SerialDenseMatrix matrix_gp_to_base(shapefcns_at_gps.N(), shapefcns_at_gps.M());
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

      LINALG::SerialDenseVector shape_functions(base_numnod);
      switch (nsd)
      {
        case 3:
        {
          DRT::UTILS::shape_function_3D(shape_functions, reference_nodes(0, dis_inode),
              reference_nodes(1, dis_inode), reference_nodes(2, dis_inode), base_distype);
        }
        break;
        case 2:
        {
          DRT::UTILS::shape_function_2D(shape_functions, reference_nodes(0, dis_inode),
              reference_nodes(1, dis_inode), base_distype);
        }
        break;
        case 1:
        {
          DRT::UTILS::shape_function_1D(
              shape_functions, reference_nodes(0, dis_inode), base_distype);
        }
        break;
        default:
          dserror("This function is not implemented for space dimension %d.", nsd);
      }
      for (int basedis_inode = 0; basedis_inode < base_numnod; ++basedis_inode)
      {
        matrix_base_to_dis(dis_inode, basedis_inode) = shape_functions(basedis_inode);
      }
    }

    LINALG::SerialDenseMatrix matrix_gp_to_nodes(
        DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, matrix_gp_to_base.N());

    // extend matrix from base_distype to distype
    matrix_gp_to_nodes.Multiply('N', 'N', 1.0, matrix_base_to_dis, matrix_gp_to_base, 0.0);

    return matrix_gp_to_nodes;
  }

  template <class T>
  void AssembleExtrapolatedNodalValues(Epetra_MultiVector& global_data, const T& nodal_data,
      const DRT::Element& ele, bool nodal_average)
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
}  // namespace


template <DRT::Element::DiscretizationType distype, class GaussIntegration>
LINALG::SerialDenseMatrix DRT::UTILS::EvaluateGaussPointsToNodesExtrapolationMatrix(
    const GaussIntegration& intpoints)
{
  static std::unordered_map<unsigned, LINALG::SerialDenseMatrix> extrapolation_matrix_cache{};

  if (extrapolation_matrix_cache.find(intpoints.NumPoints()) == extrapolation_matrix_cache.end())
  {
    DRT::Element::DiscretizationType base_distype =
        GetGaussPointExtrapolationBaseDistype<distype>(intpoints.NumPoints());

    dsassert(DRT::UTILS::getNumberOfElementNodes(base_distype) <= intpoints.NumPoints(),
        "The base discretization has more nodes than Gauss points. The extrapolation is not "
        "unique! "
        "This should not happen. The evaluation of the base extrapolation type for the number of "
        "gauss points is not correct.");

    LINALG::SerialDenseMatrix shapefcns_at_gps =
        EvaluateBaseShapeFunctionsAtGaussPoints<distype>(base_distype, intpoints);

    LINALG::SerialDenseMatrix matrix_gp_to_base =
        EvaluateProjectionGaussPointsToBaseDistype(shapefcns_at_gps);

    extrapolation_matrix_cache[intpoints.NumPoints()] =
        EvaluateProjectionGaussPointsToDistype<distype>(matrix_gp_to_base, base_distype);
  }

  return extrapolation_matrix_cache[intpoints.NumPoints()];
}

// template specialization for pyramid 5 elements
// The default procedure of extrapolation by using the shape functions results in different results
// than with our previous method. The 8 Gauss points create a HEX-element inside the pyramid. The
// extrapolation matrix holds the shapefunction-values of the HEX-element, evaluated at the
// pyramid-nodes.
template <>
LINALG::SerialDenseMatrix DRT::UTILS::EvaluateGaussPointsToNodesExtrapolationMatrix<
    DRT::Element::DiscretizationType::pyramid5>(const DRT::UTILS::IntegrationPoints3D& intpoints)
{
  if (intpoints.NumPoints() != 8)
  {
    dserror(
        "Gauss point extrapolation is not yet implemented for Pyramid5 elements with %d Gauss "
        "points. Currently, only 8 are supported",
        intpoints.NumPoints());
  }

  static LINALG::SerialDenseMatrix extrapolation_matrix = std::invoke(
      []()
      {
        LINALG::SerialDenseMatrix expol(5, 8);
        expol(0, 0) = 2.408235313815748;
        expol(0, 1) = -0.6452847075210328;
        expol(0, 2) = 0.1729035162684118;
        expol(0, 3) = -0.6452847075210328;
        expol(0, 4) = -0.542209910031327;
        expol(0, 5) = 0.1452847075210439;
        expol(0, 6) = -0.03892892005285509;
        expol(0, 7) = 0.1452847075210439;
        expol(1, 0) = -0.6452847075210328;
        expol(1, 1) = 2.408235313815748;
        expol(1, 2) = -0.6452847075210328;
        expol(1, 3) = 0.1729035162684118;
        expol(1, 4) = 0.1452847075210439;
        expol(1, 5) = -0.542209910031327;
        expol(1, 6) = 0.1452847075210439;
        expol(1, 7) = -0.03892892005285509;
        expol(2, 0) = 0.1729035162684118;
        expol(2, 1) = -0.6452847075210328;
        expol(2, 2) = 2.408235313815748;
        expol(2, 3) = -0.6452847075210328;
        expol(2, 4) = -0.03892892005285509;
        expol(2, 5) = 0.1452847075210439;
        expol(2, 6) = -0.542209910031327;
        expol(2, 7) = 0.1452847075210439;
        expol(3, 0) = -0.6452847075210328;
        expol(3, 1) = 0.1729035162684118;
        expol(3, 2) = -0.6452847075210328;
        expol(3, 3) = 2.408235313815748;
        expol(3, 4) = 0.1452847075210439;
        expol(3, 5) = -0.03892892005285509;
        expol(3, 6) = 0.1452847075210439;
        expol(3, 7) = -0.542209910031327;
        expol(4, 0) = -0.2702847075210531;
        expol(4, 1) = -0.2702847075210531;
        expol(4, 2) = -0.2702847075210531;
        expol(4, 3) = -0.2702847075210531;
        expol(4, 4) = 0.520284707521053;
        expol(4, 5) = 0.520284707521053;
        expol(4, 6) = 0.520284707521053;
        expol(4, 7) = 0.520284707521053;
        return expol;
      });

  return extrapolation_matrix;
}

template <DRT::Element::DiscretizationType distype, class GaussIntegration>
void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble(const DRT::Element& ele,
    const LINALG::SerialDenseMatrix& gp_data, Epetra_MultiVector& global_data, bool nodal_average,
    const GaussIntegration& integration)
{
  LINALG::SerialDenseMatrix nodal_quantity(
      DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, gp_data.N());
  nodal_quantity.Multiply('N', 'N', 1.0,
      EvaluateGaussPointsToNodesExtrapolationMatrix<distype>(integration), gp_data, 0.0);

  AssembleExtrapolatedNodalValues(global_data, nodal_quantity, ele, nodal_average);
}

template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::hex8, DRT::UTILS::IntegrationPoints3D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints3D& integration);
template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::hex27, DRT::UTILS::IntegrationPoints3D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints3D& integration);
template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::tet10, DRT::UTILS::IntegrationPoints3D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints3D& integration);
template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::hex20, DRT::UTILS::IntegrationPoints3D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints3D& integration);
template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::tet4, DRT::UTILS::IntegrationPoints3D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints3D& integration);
template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::wedge6, DRT::UTILS::IntegrationPoints3D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints3D& integration);
template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::wedge15, DRT::UTILS::IntegrationPoints3D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints3D& integration);
template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::pyramid5, DRT::UTILS::IntegrationPoints3D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints3D& integration);
template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::quad4, DRT::UTILS::IntegrationPoints2D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints2D& integration);
template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::quad8, DRT::UTILS::IntegrationPoints2D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints2D& integration);
template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::quad9, DRT::UTILS::IntegrationPoints2D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints2D& integration);
template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::tri3, DRT::UTILS::IntegrationPoints2D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints2D& integration);
template void DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
    DRT::Element::DiscretizationType::tri6, DRT::UTILS::IntegrationPoints2D>(
    const DRT::Element& ele, const LINALG::SerialDenseMatrix& gp_data,
    Epetra_MultiVector& global_data, bool nodal_average,
    const DRT::UTILS::IntegrationPoints2D& integration);