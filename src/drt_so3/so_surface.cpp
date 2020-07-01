/*----------------------------------------------------------------------*/
/*! \file

\brief class for evaluation of equations on the structural surface
\level 1


*----------------------------------------------------------------------*/

#include "so_surface.H"
#include "so_line.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"

DRT::ELEMENTS::StructuralSurfaceType DRT::ELEMENTS::StructuralSurfaceType::instance_;

DRT::ELEMENTS::StructuralSurfaceType& DRT::ELEMENTS::StructuralSurfaceType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::StructuralSurfaceType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::StructuralSurface* object = new DRT::ELEMENTS::StructuralSurface(-1, -1);
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
      distype_(DRT::Element::dis_none),
      numdofpernode_(-1),
      gaussrule_(DRT::UTILS::intrule2D_undefined)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lsurface);

  numdofpernode_ = ParentMasterElement()->NumDofPerNode(*Nodes()[0]);
  // Safety check if all nodes have the same number of dofs!
  for (int nlid = 1; nlid < NumNode(); ++nlid)
  {
    if (numdofpernode_ != ParentMasterElement()->NumDofPerNode(*Nodes()[nlid]))
      dserror(
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
      distype_(DRT::Element::dis_none),
      numdofpernode_(-1),
      gaussrule_(DRT::UTILS::intrule2D_undefined)
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
  DRT::ELEMENTS::StructuralSurface* newelement = new DRT::ELEMENTS::StructuralSurface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::StructuralSurface::Shape() const
{
  return distype_;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
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
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class DRT::FaceElement
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::FaceElement::Unpack(basedata);

  // distype_
  distype_ = static_cast<DRT::Element::DiscretizationType>(ExtractInt(position, data));
  // numdofpernode_
  numdofpernode_ = ExtractInt(position, data);
  // gaussrule_
  gaussrule_ = static_cast<DRT::UTILS::GaussRule2D>(ExtractInt(position, data));

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);

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
  return DRT::UTILS::ElementBoundaryFactory<DRT::ELEMENTS::StructuralLine,
      DRT::ELEMENTS::StructuralSurface>(DRT::UTILS::buildLines, this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::StructuralSurface::NumLine() const
{
  return DRT::UTILS::getNumberOfElementLines(Shape());
}

/*------------------------------------------------------------------------*
 |  Set Discretization Type of the Surface Element              ager 12/16|
 *-----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::SetDistype()
{
  // if NURBS elements:
  if (ParentMasterElement()->Shape() == nurbs8)
    distype_ = DRT::Element::nurbs4;
  else if (ParentMasterElement()->Shape() == nurbs27)
    distype_ = DRT::Element::nurbs9;
  // Lagrange elements:
  else
  {
    switch (NumNode())
    {
      case 3:
        distype_ = DRT::Element::tri3;
        break;
      case 6:
      {
        if (ParentMasterElement()->Shape() == tet10)
          distype_ = DRT::Element::tri6;
        else if (ParentMasterElement()->Shape() == hex18)
          distype_ = DRT::Element::quad6;
        else
        {
          dserror("what other surface element has 6 nodes???");
          distype_ = DRT::Element::dis_none;
        }
        break;
      }
      case 4:
        distype_ = DRT::Element::quad4;
        break;
      case 8:
        distype_ = DRT::Element::quad8;
        break;
      case 9:
        distype_ = DRT::Element::quad9;
        break;
      default:
        dserror("Unknown shape of surface element (unknown number of nodes)");
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
    case tri3:
      gaussrule_ = DRT::UTILS::intrule_tri_3point;
      break;
    case tri6:
      gaussrule_ = DRT::UTILS::intrule_tri_6point;
      break;
    case quad4:
      gaussrule_ = DRT::UTILS::intrule_quad_4point;
      break;
    case quad8:
      gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
    case quad9:
      gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
    case nurbs9:
      gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
    case quad6:
      gaussrule_ = DRT::UTILS::intrule_quad_6point;
      break;
    default:
      dserror("shape type unknown!\n");
  }
}
