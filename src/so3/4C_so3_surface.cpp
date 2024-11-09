// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_so3_surface.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_so3_line.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::Elements::StructuralSurfaceType Discret::Elements::StructuralSurfaceType::instance_;

Discret::Elements::StructuralSurfaceType& Discret::Elements::StructuralSurfaceType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::StructuralSurfaceType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::Elements::StructuralSurface(-1, -1);
  object->unpack(buffer);
  return object;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::StructuralSurfaceType::create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new StructuralSurface( id, owner ) );
  return nullptr;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                              gee 04/08|
 *----------------------------------------------------------------------*/
Discret::Elements::StructuralSurface::StructuralSurface(int id, int owner, int nnode,
    const int* nodeids, Core::Nodes::Node** nodes, Core::Elements::Element* parent,
    const int lsurface)
    : Core::Elements::FaceElement(id, owner),
      distype_(Core::FE::CellType::dis_none),
      numdofpernode_(-1),
      gaussrule_(Core::FE::GaussRule2D::undefined)
{
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  set_parent_master_element(parent, lsurface);

  numdofpernode_ = parent_master_element()->num_dof_per_node(*StructuralSurface::nodes()[0]);
  // Safety check if all nodes have the same number of dofs!
  for (int nlid = 1; nlid < num_node(); ++nlid)
  {
    if (numdofpernode_ !=
        parent_master_element()->num_dof_per_node(*StructuralSurface::nodes()[nlid]))
      FOUR_C_THROW(
          "You need different NumDofPerNode for each node on this structural surface? (%d != %d)",
          numdofpernode_,
          parent_master_element()->num_dof_per_node(*StructuralSurface::nodes()[nlid]));
  }

  set_distype();
  set_gaussrule();
  return;
}

/*------------------------------------------------------------------------*
 |  ctor (private) - used by StructuralSurfaceType              ager 12/16|
 *-----------------------------------------------------------------------*/
Discret::Elements::StructuralSurface::StructuralSurface(int id, int owner)
    : Core::Elements::FaceElement(id, owner),
      distype_(Core::FE::CellType::dis_none),
      numdofpernode_(-1),
      gaussrule_(Core::FE::GaussRule2D::undefined)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 04/08|
 *----------------------------------------------------------------------*/
Discret::Elements::StructuralSurface::StructuralSurface(
    const Discret::Elements::StructuralSurface& old)
    : Core::Elements::FaceElement(old),
      distype_(old.distype_),
      numdofpernode_(old.numdofpernode_),
      gaussrule_(old.gaussrule_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 04/08|
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::StructuralSurface::clone() const
{
  auto* newelement = new Discret::Elements::StructuralSurface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::Elements::StructuralSurface::shape() const { return distype_; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
void Discret::Elements::StructuralSurface::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Core::Elements::FaceElement
  Core::Elements::FaceElement::pack(data);
  // add distype_
  add_to_pack(data, distype_);
  // add numdofpernode_
  add_to_pack(data, numdofpernode_);
  // add gaussrule_
  add_to_pack(data, gaussrule_);
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
void Discret::Elements::StructuralSurface::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Core::Elements::FaceElement
  Core::Elements::FaceElement::unpack(buffer);

  // distype_
  extract_from_pack(buffer, distype_);
  // numdofpernode_
  extract_from_pack(buffer, numdofpernode_);
  // gaussrule_
  extract_from_pack(buffer, gaussrule_);
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                gee 04/08|
 *----------------------------------------------------------------------*/
void Discret::Elements::StructuralSurface::print(std::ostream& os) const
{
  os << "StructuralSurface ";
  Element::print(os);
  return;
}

std::vector<std::shared_ptr<Core::Elements::Element>> Discret::Elements::StructuralSurface::lines()
{
  return Core::Communication::element_boundary_factory<Discret::Elements::StructuralLine,
      Discret::Elements::StructuralSurface>(Core::Communication::buildLines, *this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Discret::Elements::StructuralSurface::num_line() const
{
  return Core::FE::get_number_of_element_lines(shape());
}

/*------------------------------------------------------------------------*
 |  Set discretization Type of the Surface Element              ager 12/16|
 *-----------------------------------------------------------------------*/
void Discret::Elements::StructuralSurface::set_distype()
{
  // if NURBS elements:
  if (parent_master_element()->shape() == Core::FE::CellType::nurbs8)
    distype_ = Core::FE::CellType::nurbs4;
  else if (parent_master_element()->shape() == Core::FE::CellType::nurbs27)
    distype_ = Core::FE::CellType::nurbs9;
  // Lagrange elements:
  else
  {
    switch (num_node())
    {
      case 3:
        distype_ = Core::FE::CellType::tri3;
        break;
      case 6:
      {
        if (parent_master_element()->shape() == Core::FE::CellType::tet10)
          distype_ = Core::FE::CellType::tri6;
        else if (parent_master_element()->shape() == Core::FE::CellType::hex18)
          distype_ = Core::FE::CellType::quad6;
        else
        {
          FOUR_C_THROW("what other surface element has 6 nodes???");
          distype_ = Core::FE::CellType::dis_none;
        }
        break;
      }
      case 4:
        distype_ = Core::FE::CellType::quad4;
        break;
      case 8:
        distype_ = Core::FE::CellType::quad8;
        break;
      case 9:
        distype_ = Core::FE::CellType::quad9;
        break;
      default:
        FOUR_C_THROW("Unknown shape of surface element (unknown number of nodes)");
        break;
    }
  }
}

/*------------------------------------------------------------------------*
 |  Set Gaussrule dependent on shape of the structural surface  ager 12/16|
 *-----------------------------------------------------------------------*/
void Discret::Elements::StructuralSurface::set_gaussrule()
{
  // type of gaussian integration
  switch (shape())
  {
    case Core::FE::CellType::tri3:
      gaussrule_ = Core::FE::GaussRule2D::tri_3point;
      break;
    case Core::FE::CellType::tri6:
      gaussrule_ = Core::FE::GaussRule2D::tri_6point;
      break;
    case Core::FE::CellType::quad4:
      gaussrule_ = Core::FE::GaussRule2D::quad_4point;
      break;
    case Core::FE::CellType::quad8:
      gaussrule_ = Core::FE::GaussRule2D::quad_9point;
      break;
    case Core::FE::CellType::quad9:
      gaussrule_ = Core::FE::GaussRule2D::quad_9point;
      break;
    case Core::FE::CellType::nurbs9:
      gaussrule_ = Core::FE::GaussRule2D::quad_9point;
      break;
    case Core::FE::CellType::quad6:
      gaussrule_ = Core::FE::GaussRule2D::quad_6point;
      break;
    default:
      FOUR_C_THROW("shape type unknown!\n");
  }
}

FOUR_C_NAMESPACE_CLOSE
