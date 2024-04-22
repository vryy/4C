/*----------------------------------------------------------------------*/
/*! \file

\brief class for evaluation of equations on the structural surface
\level 1


*----------------------------------------------------------------------*/

#include "baci_so3_surface.hpp"

#include "baci_comm_utils_factory.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_so3_line.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::StructuralSurfaceType DRT::ELEMENTS::StructuralSurfaceType::instance_;

DRT::ELEMENTS::StructuralSurfaceType& DRT::ELEMENTS::StructuralSurfaceType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::StructuralSurfaceType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::StructuralSurface(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::StructuralSurfaceType::Create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new StructuralSurface( id, owner ) );
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                              gee 04/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::StructuralSurface::StructuralSurface(int id, int owner, int nnode,
    const int* nodeids, DRT::Node** nodes, DRT::Element* parent, const int lsurface)
    : DRT::FaceElement(id, owner),
      distype_(CORE::FE::CellType::dis_none),
      numdofpernode_(-1),
      gaussrule_(CORE::FE::GaussRule2D::undefined)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lsurface);

  numdofpernode_ = ParentMasterElement()->NumDofPerNode(*Nodes()[0]);
  // Safety check if all nodes have the same number of dofs!
  for (int nlid = 1; nlid < NumNode(); ++nlid)
  {
    if (numdofpernode_ != ParentMasterElement()->NumDofPerNode(*Nodes()[nlid]))
      FOUR_C_THROW(
          "You need different NumDofPerNode for each node on this structural surface? (%d != %d)",
          numdofpernode_, ParentMasterElement()->NumDofPerNode(*Nodes()[nlid]));
  }

  SetDistype();
  SetGaussrule();
  return;
}

/*------------------------------------------------------------------------*
 |  ctor (private) - used by StructuralSurfaceType              ager 12/16|
 *-----------------------------------------------------------------------*/
DRT::ELEMENTS::StructuralSurface::StructuralSurface(int id, int owner)
    : DRT::FaceElement(id, owner),
      distype_(CORE::FE::CellType::dis_none),
      numdofpernode_(-1),
      gaussrule_(CORE::FE::GaussRule2D::undefined)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 04/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::StructuralSurface::StructuralSurface(const DRT::ELEMENTS::StructuralSurface& old)
    : DRT::FaceElement(old),
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
DRT::Element* DRT::ELEMENTS::StructuralSurface::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::StructuralSurface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::StructuralSurface::Shape() const { return distype_; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class DRT::FaceElement
  DRT::FaceElement::Pack(data);
  // add distype_
  AddtoPack(data, (int)distype_);
  // add numdofpernode_
  AddtoPack(data, numdofpernode_);
  // add gaussrule_
  AddtoPack(data, (int)gaussrule_);
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class DRT::FaceElement
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::FaceElement::Unpack(basedata);

  // distype_
  distype_ = static_cast<CORE::FE::CellType>(ExtractInt(position, data));
  // numdofpernode_
  numdofpernode_ = ExtractInt(position, data);
  // gaussrule_
  gaussrule_ = static_cast<CORE::FE::GaussRule2D>(ExtractInt(position, data));

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::Print(std::ostream& os) const
{
  os << "StructuralSurface ";
  Element::Print(os);
  return;
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::StructuralSurface::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<DRT::ELEMENTS::StructuralLine,
      DRT::ELEMENTS::StructuralSurface>(CORE::COMM::buildLines, *this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::StructuralSurface::NumLine() const
{
  return CORE::FE::getNumberOfElementLines(Shape());
}

/*------------------------------------------------------------------------*
 |  Set Discretization Type of the Surface Element              ager 12/16|
 *-----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::SetDistype()
{
  // if NURBS elements:
  if (ParentMasterElement()->Shape() == CORE::FE::CellType::nurbs8)
    distype_ = CORE::FE::CellType::nurbs4;
  else if (ParentMasterElement()->Shape() == CORE::FE::CellType::nurbs27)
    distype_ = CORE::FE::CellType::nurbs9;
  // Lagrange elements:
  else
  {
    switch (NumNode())
    {
      case 3:
        distype_ = CORE::FE::CellType::tri3;
        break;
      case 6:
      {
        if (ParentMasterElement()->Shape() == CORE::FE::CellType::tet10)
          distype_ = CORE::FE::CellType::tri6;
        else if (ParentMasterElement()->Shape() == CORE::FE::CellType::hex18)
          distype_ = CORE::FE::CellType::quad6;
        else
        {
          FOUR_C_THROW("what other surface element has 6 nodes???");
          distype_ = CORE::FE::CellType::dis_none;
        }
        break;
      }
      case 4:
        distype_ = CORE::FE::CellType::quad4;
        break;
      case 8:
        distype_ = CORE::FE::CellType::quad8;
        break;
      case 9:
        distype_ = CORE::FE::CellType::quad9;
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
void DRT::ELEMENTS::StructuralSurface::SetGaussrule()
{
  // type of gaussian integration
  switch (Shape())
  {
    case CORE::FE::CellType::tri3:
      gaussrule_ = CORE::FE::GaussRule2D::tri_3point;
      break;
    case CORE::FE::CellType::tri6:
      gaussrule_ = CORE::FE::GaussRule2D::tri_6point;
      break;
    case CORE::FE::CellType::quad4:
      gaussrule_ = CORE::FE::GaussRule2D::quad_4point;
      break;
    case CORE::FE::CellType::quad8:
      gaussrule_ = CORE::FE::GaussRule2D::quad_9point;
      break;
    case CORE::FE::CellType::quad9:
      gaussrule_ = CORE::FE::GaussRule2D::quad_9point;
      break;
    case CORE::FE::CellType::nurbs9:
      gaussrule_ = CORE::FE::GaussRule2D::quad_9point;
      break;
    case CORE::FE::CellType::quad6:
      gaussrule_ = CORE::FE::GaussRule2D::quad_6point;
      break;
    default:
      FOUR_C_THROW("shape type unknown!\n");
  }
}

FOUR_C_NAMESPACE_CLOSE
