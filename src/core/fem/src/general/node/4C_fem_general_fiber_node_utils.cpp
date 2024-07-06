/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of utility functions for fiber interpolation

\level 3
*/
/*----------------------------------------------------------------------*/


#include "4C_fem_general_fiber_node_utils.hpp"

#include "4C_fem_general_fiber_node.hpp"
#include "4C_fem_general_fiber_node_holder.hpp"

#include <cmath>

FOUR_C_NAMESPACE_OPEN

template <Core::FE::CellType distype>
void Core::Nodes::ProjectFibersToGaussPoints(const Core::Nodes::Node* const* nodes,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, 1>>& shapefcts,
    NodalFiberHolder& gpFiberHolder)
{
  // number of nodes per element
  constexpr int nen = Core::FE::num_nodes<distype>;
  std::array<const Core::Nodes::FiberNode*, nen> fiberNodes;
  std::vector<std::array<std::array<double, 3>, nen>> fibers;
  std::map<Core::Nodes::CoordinateSystemDirection, std::array<std::array<double, 3>, nen>>
      coordinateSystemDirections;
  std::map<Core::Nodes::AngleType, std::array<double, nen>> angles;

  for (int inode = 0; inode < nen; ++inode)
  {
    fiberNodes[inode] = dynamic_cast<const Core::Nodes::FiberNode*>(nodes[inode]);

    if (fiberNodes[inode] == nullptr)
    {
      FOUR_C_THROW("At least one node of the element does not provide fibers.");
    }

    for (const auto& pair : fiberNodes[inode]->coordinate_system_directions())
    {
      coordinateSystemDirections[pair.first][inode] = pair.second;
    }

    for (std::size_t fiberid = 0; fiberid < fiberNodes[inode]->fibers().size(); ++fiberid)
    {
      if (inode == 0)
      {
        // The first node has to create the array of fibers for each node
        fibers.emplace_back(std::array<std::array<double, 3>, nen>());
      }
      fibers[fiberid][inode] = fiberNodes[inode]->fibers()[fiberid];
    }

    for (const auto& pair : fiberNodes[inode]->angles())
    {
      angles[pair.first][inode] = pair.second;
    }
  }


  // project fibers
  for (const auto& fiber : fibers)
  {
    std::vector<Core::LinAlg::Matrix<3, 1>> gpQuantity;
    ProjectQuantityWithShapeFunctions<distype, 3>(fiber, shapefcts, gpQuantity);

    // normalize fiber vectors
    for (Core::LinAlg::Matrix<3, 1>& quantity : gpQuantity)
    {
      quantity.scale(1.0 / quantity.norm2());
    }

    // add quantity to the fiber holder
    gpFiberHolder.add_fiber(gpQuantity);
  }


  // project coordinate system directions
  for (const auto& pair : coordinateSystemDirections)
  {
    const Core::Nodes::CoordinateSystemDirection type = pair.first;
    std::vector<Core::LinAlg::Matrix<3, 1>> gpQuantity;
    ProjectQuantityWithShapeFunctions<distype, 3>(pair.second, shapefcts, gpQuantity);

    // normalize fiber vectors
    for (Core::LinAlg::Matrix<3, 1>& quantity : gpQuantity)
    {
      quantity.scale(1.0 / quantity.norm2());
    }

    // add quantity to the fiber holder
    gpFiberHolder.set_coordinate_system_direction(type, gpQuantity);
  }

  // project angles
  for (const auto& pair : angles)
  {
    const Core::Nodes::AngleType type = pair.first;
    std::vector<double> gpAngle;
    ProjectQuantityWithShapeFunctions<distype>(pair.second, shapefcts, gpAngle);

    // add quantity to the fiber holder
    gpFiberHolder.set_angle(type, gpAngle);
  }


  // orthogonalize coordinate system
  if (gpFiberHolder.contains_coordinate_system_direction(CoordinateSystemDirection::Circular) &&
      gpFiberHolder.contains_coordinate_system_direction(CoordinateSystemDirection::Tangential))
  {
    const std::vector<Core::LinAlg::Matrix<3, 1>>& cir =
        gpFiberHolder.get_coordinate_system_direction(CoordinateSystemDirection::Circular);
    std::vector<Core::LinAlg::Matrix<3, 1>>& tan =
        gpFiberHolder.get_coordinate_system_direction_mutual(CoordinateSystemDirection::Tangential);

    // orthogonalize tangential vector, preserve circular direction
    for (std::size_t gp = 0; gp < tan.size(); ++gp)
    {
      double tancir = tan[gp].dot(cir[gp]);
      tan[gp].update(-tancir, cir[gp], 1.0);
      tan[gp].scale(1.0 / tan[gp].norm2());
    }

    // orthogonalize radial vector, preserve circular and tangential direction
    if (gpFiberHolder.contains_coordinate_system_direction(CoordinateSystemDirection::Radial))
    {
      std::vector<Core::LinAlg::Matrix<3, 1>>& rad =
          gpFiberHolder.get_coordinate_system_direction_mutual(CoordinateSystemDirection::Radial);
      for (std::size_t gp = 0; gp < tan.size(); ++gp)
      {
        double radcir = rad[gp].dot(cir[gp]);
        double radtan = rad[gp].dot(tan[gp]);
        // double
        rad[gp].update(-radcir, cir[gp], -radtan, tan[gp], 1.0);
        rad[gp].scale(1.0 / rad[gp].norm2());
      }
    }
  }
}

template <Core::FE::CellType distype, std::size_t dim>
void Core::Nodes::ProjectQuantityWithShapeFunctions(
    const std::array<std::array<double, dim>, Core::FE::num_nodes<distype>> quantity,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, 1>>& shapefcts,
    std::vector<Core::LinAlg::Matrix<dim, 1>>& quantityProjected)
{
  quantityProjected.resize(shapefcts.size());

  // number of nodes per element
  constexpr std::size_t nen = Core::FE::num_nodes<distype>;

  // spatial dimension
  constexpr std::size_t nsd = Core::FE::dim<distype>;

  // number of gauss points
  const std::size_t ngp = shapefcts.size();

  for (Core::LinAlg::Matrix<dim, 1>& quantity : quantityProjected)
  {
    quantity.clear();
  }

  for (std::size_t i = 0; i < nsd; ++i)
  {
    for (std::size_t q = 0; q < ngp; ++q)
    {
      for (std::size_t j = 0; j < nen; ++j)
      {
        quantityProjected[q](i) += shapefcts[q](j) * quantity[j][i];
      }
    }
  }
}

template <Core::FE::CellType distype>
void Core::Nodes::ProjectQuantityWithShapeFunctions(
    const std::array<double, Core::FE::num_nodes<distype>> quantity,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, 1>>& shapefcts,
    std::vector<double>& quantityProjected)
{
  quantityProjected.resize(shapefcts.size());

  // number of nodes per element
  constexpr int nen = Core::FE::num_nodes<distype>;

  // number of gauss points
  const int ngp = shapefcts.size();

  for (double& quantity : quantityProjected)
  {
    quantity = 0.0;
  }

  for (int q = 0; q < ngp; ++q)
  {
    for (int j = 0; j < nen; ++j)
    {
      quantityProjected[q] += shapefcts[q](j) * quantity[j];
    }
  }
}

template <Core::FE::CellType distype>
bool Core::Nodes::HaveNodalFibers(const Core::Nodes::Node* const* nodes)
{
  constexpr int numberOfNodes = Core::FE::num_nodes<distype>;

  // check whether node can be casted to FiberNode
  auto nodeHasFibers = [](const Core::Nodes::Node* n)
  { return dynamic_cast<const Core::Nodes::FiberNode*>(n) != nullptr; };
  return std::all_of(nodes, &nodes[numberOfNodes], nodeHasFibers);
}

template void Core::Nodes::ProjectQuantityWithShapeFunctions<Core::FE::CellType::tet4>(
    const std::array<double, Core::FE::num_nodes<Core::FE::CellType::tet4>> quantity,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::tet4>, 1>>&
        shapefcts,
    std::vector<double>& quantityProjected);
template void Core::Nodes::ProjectQuantityWithShapeFunctions<Core::FE::CellType::tet10>(
    const std::array<double, Core::FE::num_nodes<Core::FE::CellType::tet10>> quantity,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::tet10>, 1>>&
        shapefcts,
    std::vector<double>& quantityProjected);
template void Core::Nodes::ProjectQuantityWithShapeFunctions<Core::FE::CellType::hex8>(
    const std::array<double, Core::FE::num_nodes<Core::FE::CellType::hex8>> quantity,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::hex8>, 1>>&
        shapefcts,
    std::vector<double>& quantityProjected);
template void Core::Nodes::ProjectQuantityWithShapeFunctions<Core::FE::CellType::tri3>(
    const std::array<double, Core::FE::num_nodes<Core::FE::CellType::tri3>> quantity,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::tri3>, 1>>&
        shapefcts,
    std::vector<double>& quantityProjected);

template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::tet4>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::tet4>, 1>>&,
    NodalFiberHolder&);
template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::tet10>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::tet10>, 1>>&,
    NodalFiberHolder&);
template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::hex8>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::hex8>, 1>>&,
    NodalFiberHolder&);
template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::hex18>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::hex18>, 1>>&,
    NodalFiberHolder&);
template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::hex20>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::hex20>, 1>>&,
    NodalFiberHolder&);
template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::hex27>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::hex27>, 1>>&,
    NodalFiberHolder&);
template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::nurbs27>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::nurbs27>, 1>>&,
    NodalFiberHolder&);
template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::tri3>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::tri3>, 1>>&,
    NodalFiberHolder&);
template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::quad4>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::quad4>, 1>>&,
    NodalFiberHolder&);
template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::quad9>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::quad9>, 1>>&,
    NodalFiberHolder&);
template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::pyramid5>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::pyramid5>, 1>>&,
    NodalFiberHolder&);
template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::wedge6>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::wedge6>, 1>>&,
    NodalFiberHolder&);
template void Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::nurbs9>(
    const Core::Nodes::Node* const*,
    const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::nurbs9>, 1>>&,
    NodalFiberHolder&);

template bool Core::Nodes::HaveNodalFibers<Core::FE::CellType::tet4>(
    const Core::Nodes::Node* const* nodes);
template bool Core::Nodes::HaveNodalFibers<Core::FE::CellType::tet10>(
    const Core::Nodes::Node* const* nodes);
template bool Core::Nodes::HaveNodalFibers<Core::FE::CellType::hex8>(
    const Core::Nodes::Node* const* nodes);
template bool Core::Nodes::HaveNodalFibers<Core::FE::CellType::hex18>(
    const Core::Nodes::Node* const* nodes);
template bool Core::Nodes::HaveNodalFibers<Core::FE::CellType::hex20>(
    const Core::Nodes::Node* const* nodes);
template bool Core::Nodes::HaveNodalFibers<Core::FE::CellType::hex27>(
    const Core::Nodes::Node* const* nodes);
template bool Core::Nodes::HaveNodalFibers<Core::FE::CellType::nurbs27>(
    const Core::Nodes::Node* const* nodes);
template bool Core::Nodes::HaveNodalFibers<Core::FE::CellType::pyramid5>(
    const Core::Nodes::Node* const* nodes);
template bool Core::Nodes::HaveNodalFibers<Core::FE::CellType::wedge6>(
    const Core::Nodes::Node* const* nodes);

FOUR_C_NAMESPACE_CLOSE
