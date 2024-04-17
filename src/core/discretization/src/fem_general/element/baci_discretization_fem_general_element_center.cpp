/*----------------------------------------------------------------------*/
/*! \file
\brief Center coordinates of an element
\level 0
*/
/*----------------------------------------------------------------------*/

#include "baci_discretization_fem_general_element_center.hpp"

#include "baci_lib_element.hpp"
#include "baci_lib_node.hpp"

FOUR_C_NAMESPACE_OPEN

std::vector<double> CORE::FE::ElementCenterRefeCoords(const DRT::Element& ele)
{
  // get nodes of element
  const DRT::Node* const* nodes = ele.Nodes();
  const int numnodes = ele.NumNode();
  const double invnumnodes = 1.0 / numnodes;

  // calculate mean of node coordinates
  std::vector<double> centercoords(3, 0.0);
  for (int i = 0; i < 3; ++i)
  {
    double var = 0.0;
    for (int j = 0; j < numnodes; ++j)
    {
      const auto& x = nodes[j]->X();
      var += x[i];
    }
    centercoords[i] = var * invnumnodes;
  }

  return centercoords;
}

FOUR_C_NAMESPACE_CLOSE
