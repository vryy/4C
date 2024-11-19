// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_nodereader.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element_definition.hpp"
#include "4C_fem_general_fiber_node.hpp"
#include "4C_fem_general_immersed_node.hpp"
#include "4C_fem_nurbs_discretization_control_point.hpp"
#include "4C_io_input_file.hpp"

#include <istream>

FOUR_C_NAMESPACE_OPEN

namespace
{
  std::vector<std::shared_ptr<Core::FE::Discretization>> find_dis_node(
      const std::vector<Core::IO::ElementReader>& element_readers, int global_node_id)
  {
    std::vector<std::shared_ptr<Core::FE::Discretization>> list_of_discretizations;
    for (const auto& element_reader : element_readers)
      if (element_reader.has_node(global_node_id))
        list_of_discretizations.emplace_back(element_reader.get_dis());

    return list_of_discretizations;
  }

}  // namespace


void Core::IO::read_nodes(Core::IO::InputFile& input, const std::string& node_section_name,
    std::vector<ElementReader>& element_readers, int& max_node_id)
{
  const int myrank = Core::Communication::my_mpi_rank(input.get_comm());
  if (myrank > 0) return;

  std::string tmp;
  std::string tmp2;

  int line_count = 0;
  for (const auto& node_line : input.lines_in_section(node_section_name))
  {
    std::istringstream linestream{std::string{node_line}};
    linestream >> tmp;

    if (tmp == "NODE")
    {
      std::vector<double> coords(3, 0.0);
      int nodeid;
      // read in the node coordinates
      linestream >> nodeid >> tmp >> coords[0] >> coords[1] >> coords[2];

      nodeid--;
      max_node_id = std::max(max_node_id, nodeid) + 1;
      std::vector<std::shared_ptr<Core::FE::Discretization>> dis =
          find_dis_node(element_readers, nodeid);

      for (const auto& di : dis)
      {
        // create node and add to discretization
        std::shared_ptr<Core::Nodes::Node> node =
            std::make_shared<Core::Nodes::Node>(nodeid, coords, myrank);
        di->add_node(node);
      }
    }
    // this is a specialized node for immersed problems
    else if (tmp == "INODE")
    {
      std::vector<double> coords(3, 0.0);
      int nodeid;
      // read in the node coordinates
      linestream >> nodeid >> tmp >> coords[0] >> coords[1] >> coords[2];

      nodeid--;
      max_node_id = std::max(max_node_id, nodeid) + 1;
      std::vector<std::shared_ptr<Core::FE::Discretization>> diss =
          find_dis_node(element_readers, nodeid);

      for (const auto& dis : diss)
      {
        // create node and add to discretization
        std::shared_ptr<Core::Nodes::Node> node =
            std::make_shared<Core::Nodes::ImmersedNode>(nodeid, coords, myrank);
        dis->add_node(node);
      }
    }
    // this node is a Nurbs control point
    else if (tmp == "CP")
    {
      // read control points for isogeometric analysis (Nurbs)
      std::vector<double> coords(3, 0.0);
      double weight;

      int cpid;
      linestream >> cpid >> tmp >> coords[0] >> coords[1] >> coords[2] >> weight;
      cpid--;
      max_node_id = std::max(max_node_id, cpid) + 1;
      if (cpid != line_count)
        FOUR_C_THROW(
            "Reading of control points %d failed: They must be numbered consecutive!!", cpid);
      if (tmp != "COORD") FOUR_C_THROW("failed to read control point %d", cpid);
      std::vector<std::shared_ptr<Core::FE::Discretization>> diss =
          find_dis_node(element_readers, cpid);

      for (auto& dis : diss)
      {
        // create node/control point and add to discretization
        std::shared_ptr<Core::FE::Nurbs::ControlPoint> node =
            std::make_shared<Core::FE::Nurbs::ControlPoint>(cpid, coords, weight, myrank);
        dis->add_node(node);
      }
    }
    // this is a special node with additional fiber information
    else if (tmp == "FNODE")
    {
      enum class FiberType
      {
        Unknown,
        Angle,
        Fiber,
        CosyDirection
      };

      // read fiber node
      std::vector<double> coords(3, 0.0);
      std::map<Core::Nodes::CoordinateSystemDirection, std::array<double, 3>> cosyDirections;
      std::vector<std::array<double, 3>> fibers;
      std::map<Core::Nodes::AngleType, double> angles;

      int nodeid;
      // read in the node coordinates and fiber direction
      linestream >> nodeid >> tmp >> coords[0] >> coords[1] >> coords[2];
      nodeid--;
      max_node_id = std::max(max_node_id, nodeid) + 1;

      while (true)
      {
        // store current position of linestream reader
        auto length = linestream.tellg();
        // try to read new fiber direction or coordinate system
        linestream >> tmp2;

        // Nothing more to read.
        if (linestream.fail()) break;

        Core::Nodes::CoordinateSystemDirection coordinateSystemDirection;
        Core::Nodes::AngleType angleType;
        FiberType type = FiberType::Unknown;

        if (tmp2 == "FIBER" + std::to_string(1 + fibers.size()))
        {
          type = FiberType::Fiber;
        }
        else if (tmp2 == "CIR")
        {
          coordinateSystemDirection = Core::Nodes::CoordinateSystemDirection::Circular;
          type = FiberType::CosyDirection;
        }
        else if (tmp2 == "TAN")
        {
          coordinateSystemDirection = Core::Nodes::CoordinateSystemDirection::Tangential;
          type = FiberType::CosyDirection;
        }
        else if (tmp2 == "RAD")
        {
          coordinateSystemDirection = Core::Nodes::CoordinateSystemDirection::Radial;
          type = FiberType::CosyDirection;
        }
        else if (tmp2 == "HELIX")
        {
          angleType = Core::Nodes::AngleType::Helix;
          type = FiberType::Angle;
        }
        else if (tmp2 == "TRANS")
        {
          angleType = Core::Nodes::AngleType::Transverse;
          type = FiberType::Angle;
        }
        else
        {
          // No more fiber information. Jump to last position.
          linestream.seekg(length);
          break;
        }

        // add fiber / angle to the map
        switch (type)
        {
          case FiberType::Unknown:
          {
            FOUR_C_THROW(
                "Unknown fiber node attribute. Numbered fibers must be in order, i.e. "
                "FIBER1, FIBER2, ...");
          }
          case FiberType::Angle:
          {
            linestream >> angles[angleType];
            break;
          }
          case FiberType::Fiber:
          {
            std::array<double, 3> fiber_components;
            linestream >> fiber_components[0] >> fiber_components[1] >> fiber_components[2];
            fibers.emplace_back(fiber_components);
            break;
          }
          case FiberType::CosyDirection:
          {
            linestream >> cosyDirections[coordinateSystemDirection][0] >>
                cosyDirections[coordinateSystemDirection][1] >>
                cosyDirections[coordinateSystemDirection][2];
            break;
          }
          default:
            FOUR_C_THROW("Unknown number of components");
        }
      }

      // add fiber information to node
      std::vector<std::shared_ptr<Core::FE::Discretization>> discretizations =
          find_dis_node(element_readers, nodeid);
      for (auto& dis : discretizations)
      {
        auto node = std::make_shared<Core::Nodes::FiberNode>(
            nodeid, coords, cosyDirections, fibers, angles, myrank);
        dis->add_node(node);
      }
    }
    else
      FOUR_C_THROW("unexpected word '%s'", tmp.c_str());

    ++line_count;
  }
}

FOUR_C_NAMESPACE_CLOSE
