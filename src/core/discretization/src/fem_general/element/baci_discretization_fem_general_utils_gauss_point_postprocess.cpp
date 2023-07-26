/*----------------------------------------------------------------------*/
/*! \file
\brief Extrapolation of Gauss point quantities for to

\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_discretization_fem_general_utils_gauss_point_postprocess.H"
#include "baci_lib_element.H"
#include "baci_discretization_fem_general_utils_gauss_point_extrapolation.H"
#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_discretization_fem_general_utils_integration.H"

namespace
{
  template <DRT::Element::DiscretizationType distype>
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
    case ::DRT::Element::DiscretizationType::hex8:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<::DRT::Element::DiscretizationType::hex8>(
          ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::hex8>(data.M()));
    }
    break;
    case ::DRT::Element::DiscretizationType::hex27:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
          ::DRT::Element::DiscretizationType::hex27>(ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::hex27>(data.M()));
    }
    break;
    case ::DRT::Element::DiscretizationType::tet10:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
          ::DRT::Element::DiscretizationType::tet10>(ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::tet10>(data.M()));
    }
    break;
    case ::DRT::Element::DiscretizationType::tet4:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<::DRT::Element::DiscretizationType::tet4>(
          ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::tet4>(data.M()));
    }
    break;
    case ::DRT::Element::DiscretizationType::hex20:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
          ::DRT::Element::DiscretizationType::hex20>(ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::hex20>(data.M()));
    }
    break;
    case ::DRT::Element::DiscretizationType::wedge6:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
          ::DRT::Element::DiscretizationType::wedge6>(ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::wedge6>(data.M()));
    }
    break;
    case ::DRT::Element::DiscretizationType::wedge15:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
          ::DRT::Element::DiscretizationType::wedge15>(ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::wedge15>(data.M()));
    }
    break;
    case ::DRT::Element::DiscretizationType::pyramid5:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
          ::DRT::Element::DiscretizationType::pyramid5>(ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::pyramid5>(data.M()));
    }
    break;
    case ::DRT::Element::DiscretizationType::quad4:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
          ::DRT::Element::DiscretizationType::quad4>(ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::quad4>(data.M()));
    }
    break;
    case ::DRT::Element::DiscretizationType::quad8:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
          ::DRT::Element::DiscretizationType::quad8>(ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::quad8>(data.M()));
    }
    break;
    case ::DRT::Element::DiscretizationType::quad9:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
          ::DRT::Element::DiscretizationType::quad9>(ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::quad9>(data.M()));
    }
    break;
    case ::DRT::Element::DiscretizationType::tri3:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<::DRT::Element::DiscretizationType::tri3>(
          ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::tri3>(data.M()));
    }
    break;
    case ::DRT::Element::DiscretizationType::tri6:
    {
      DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<::DRT::Element::DiscretizationType::tri6>(
          ele, data, nodal_data, true,
          GetGaussIntegration<::DRT::Element::DiscretizationType::tri6>(data.M()));
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
