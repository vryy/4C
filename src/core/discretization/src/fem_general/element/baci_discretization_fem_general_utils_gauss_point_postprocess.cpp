/*----------------------------------------------------------------------*/
/*! \file
\brief Extrapolation of Gauss point quantities for to

\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_discretization_fem_general_utils_gauss_point_postprocess.H"

#include "baci_discretization_fem_general_utils_gauss_point_extrapolation.H"
#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_lib_element.H"

namespace
{
  template <CORE::FE::CellType distype>
  const auto& GetGaussIntegration(unsigned numgp)
  {
    static std::unordered_map<unsigned,
        CORE::DRT::UTILS::IntegrationPoints<CORE::DRT::UTILS::DisTypeToDim<distype>::dim>>
        gaussIntegrations = {};

    if (gaussIntegrations.find(numgp) == gaussIntegrations.end())
    {
      auto rule = CORE::DRT::UTILS::NumGaussPointsToGaussRule<distype>(numgp);
      gaussIntegrations.emplace(numgp, rule);
    }

    return gaussIntegrations.at(numgp);
  }
}  // namespace

void CORE::DRT::ELEMENTS::ExtrapolateGaussPointQuantityToNodes(::DRT::Element& ele,
    const CORE::LINALG::SerialDenseMatrix& data, Epetra_MultiVector& nodal_data)
{
  switch (ele.Shape())
  {
    case CORE::FE::CellType::hex8:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::hex8>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::hex8>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::hex27:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::hex27>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::hex27>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::tet10:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::tet10>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::tet10>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::tet4:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::tet4>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::tet4>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::hex20:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::hex20>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::hex20>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::wedge6:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::wedge6>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::wedge6>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::wedge15:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::wedge15>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::wedge15>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::pyramid5:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::pyramid5>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::pyramid5>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::quad4:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::quad4>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::quad4>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::quad8:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::quad8>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::quad8>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::quad9:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::quad9>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::quad9>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::tri3:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::tri3>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::tri3>(data.numRows()));
    }
    break;
    case CORE::FE::CellType::tri6:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::tri6>(ele, data,
          nodal_data, true, GetGaussIntegration<CORE::FE::CellType::tri6>(data.numRows()));
    }
    break;
    default:
      dserror("Your discretization type (%s) is not yet in the list!",
          ::DRT::DistypeToString(ele.Shape()).c_str());
  }
}

void CORE::DRT::ELEMENTS::EvaluateGaussPointQuantityAtElementCenter(::DRT::Element& ele,
    const CORE::LINALG::SerialDenseMatrix& data, Epetra_MultiVector& element_data)
{
  AssembleAveragedElementValues(element_data, data, ele);
}
