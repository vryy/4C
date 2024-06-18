/*----------------------------------------------------------------------*/
/*! \file

\brief class for evaluation of equations on the structural surface
\level 1


*----------------------------------------------------------------------*/

#include "4C_so3_surface.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_so3_line.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::StructuralSurfaceType Discret::ELEMENTS::StructuralSurfaceType::instance_;

Discret::ELEMENTS::StructuralSurfaceType& Discret::ELEMENTS::StructuralSurfaceType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::StructuralSurfaceType::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::StructuralSurface(-1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::StructuralSurfaceType::Create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new StructuralSurface( id, owner ) );
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                              gee 04/08|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::StructuralSurface::StructuralSurface(int id, int owner, int nnode,
    const int* nodeids, Core::Nodes::Node** nodes, Core::Elements::Element* parent,
    const int lsurface)
    : Core::Elements::FaceElement(id, owner),
      distype_(Core::FE::CellType::dis_none),
      numdofpernode_(-1),
      gaussrule_(Core::FE::GaussRule2D::undefined)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lsurface);

  numdofpernode_ = ParentMasterElement()->NumDofPerNode(*Nodes()[0]);
  // Safety check if all nodes have the same number of dofs!
  for (int nlid = 1; nlid < num_node(); ++nlid)
  {
    if (numdofpernode_ != ParentMasterElement()->NumDofPerNode(*Nodes()[nlid]))
      FOUR_C_THROW(
          "You need different NumDofPerNode for each node on this structural surface? (%d != %d)",
          numdofpernode_, ParentMasterElement()->NumDofPerNode(*Nodes()[nlid]));
  }

  set_distype();
  set_gaussrule();
  return;
}

/*------------------------------------------------------------------------*
 |  ctor (private) - used by StructuralSurfaceType              ager 12/16|
 *-----------------------------------------------------------------------*/
Discret::ELEMENTS::StructuralSurface::StructuralSurface(int id, int owner)
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
Discret::ELEMENTS::StructuralSurface::StructuralSurface(
    const Discret::ELEMENTS::StructuralSurface& old)
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
Core::Elements::Element* Discret::ELEMENTS::StructuralSurface::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::StructuralSurface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::StructuralSurface::Shape() const { return distype_; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::StructuralSurface::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Core::Elements::FaceElement
  Core::Elements::FaceElement::pack(data);
  // add distype_
  add_to_pack(data, (int)distype_);
  // add numdofpernode_
  add_to_pack(data, numdofpernode_);
  // add gaussrule_
  add_to_pack(data, (int)gaussrule_);
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::StructuralSurface::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Core::Elements::FaceElement
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Core::Elements::FaceElement::unpack(basedata);

  // distype_
  distype_ = static_cast<Core::FE::CellType>(extract_int(position, data));
  // numdofpernode_
  numdofpernode_ = extract_int(position, data);
  // gaussrule_
  gaussrule_ = static_cast<Core::FE::GaussRule2D>(extract_int(position, data));

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                gee 04/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::StructuralSurface::Print(std::ostream& os) const
{
  os << "StructuralSurface ";
  Element::Print(os);
  return;
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::StructuralSurface::Lines()
{
  return Core::Communication::ElementBoundaryFactory<Discret::ELEMENTS::StructuralLine,
      Discret::ELEMENTS::StructuralSurface>(Core::Communication::buildLines, *this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::StructuralSurface::NumLine() const
{
  return Core::FE::getNumberOfElementLines(Shape());
}

/*------------------------------------------------------------------------*
 |  Set discretization Type of the Surface Element              ager 12/16|
 *-----------------------------------------------------------------------*/
void Discret::ELEMENTS::StructuralSurface::set_distype()
{
  // if NURBS elements:
  if (ParentMasterElement()->Shape() == Core::FE::CellType::nurbs8)
    distype_ = Core::FE::CellType::nurbs4;
  else if (ParentMasterElement()->Shape() == Core::FE::CellType::nurbs27)
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
        if (ParentMasterElement()->Shape() == Core::FE::CellType::tet10)
          distype_ = Core::FE::CellType::tri6;
        else if (ParentMasterElement()->Shape() == Core::FE::CellType::hex18)
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
void Discret::ELEMENTS::StructuralSurface::set_gaussrule()
{
  // type of gaussian integration
  switch (Shape())
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
