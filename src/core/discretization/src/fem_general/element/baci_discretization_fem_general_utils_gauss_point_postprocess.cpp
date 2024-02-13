/*----------------------------------------------------------------------*/
/*! \file
\brief Extrapolation of Gauss point quantities for to

\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_discretization_fem_general_utils_gauss_point_postprocess.hpp"

#include "baci_discretization_fem_general_utils_gauss_point_extrapolation.hpp"
#include "baci_discretization_fem_general_utils_gausspoints.hpp"
#include "baci_discretization_fem_general_utils_integration.hpp"
#include "baci_lib_element.hpp"

BACI_NAMESPACE_OPEN

namespace
{
  template <CORE::FE::CellType distype>
  const auto& GetGaussIntegration(unsigned numgp)
  {
    static std::unordered_map<unsigned, CORE::FE::IntegrationPoints<CORE::FE::dim<distype>>>
        gaussIntegrations = {};

    if (gaussIntegrations.find(numgp) == gaussIntegrations.end())
    {
      auto rule = CORE::FE::NumGaussPointsToGaussRule<distype>(numgp);
      gaussIntegrations.emplace(numgp, rule);
    }

    return gaussIntegrations.at(numgp);
  }
}  // namespace

void CORE::FE::ExtrapolateGaussPointQuantityToNodes(DRT::Element& ele,
    const CORE::LINALG::SerialDenseMatrix& data, const DRT::Discretization& dis,
    Epetra_MultiVector& nodal_data)
{
  switch (ele.Shape())
  {
    case CORE::FE::CellType::hex8:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::hex8>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::hex8>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::hex27:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::hex27>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::hex27>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::nurbs27:
    {
      CORE::FE::ExtrapolateGPQuantityToNURBSKnotsAndAssemble<CORE::FE::CellType::nurbs27>(dis, ele,
          data, nodal_data, true, GetGaussIntegration<CORE::FE::CellType::nurbs27>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::tet10:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::tet10>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::tet10>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::tet4:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::tet4>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::tet4>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::hex20:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::hex20>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::hex20>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::wedge6:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::wedge6>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::wedge6>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::wedge15:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::wedge15>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::wedge15>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::pyramid5:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::pyramid5>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::pyramid5>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::quad4:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::quad4>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::quad4>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::quad8:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::quad8>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::quad8>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::quad9:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::quad9>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::quad9>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::tri3:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::tri3>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::tri3>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::tri6:
    {
      CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::tri6>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::tri6>(data.numRows()));
    }
    break;
    default:
      dserror("Your discretization type (%s) is not yet in the list!",
          CORE::FE::CellTypeToString(ele.Shape()).c_str());
  }
}

void CORE::FE::EvaluateGaussPointQuantityAtElementCenter(DRT::Element& ele,
    const CORE::LINALG::SerialDenseMatrix& data, Epetra_MultiVector& element_data)
{
  AssembleAveragedElementValues(element_data, data, ele);
}

BACI_NAMESPACE_CLOSE
