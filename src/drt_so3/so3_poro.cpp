/*----------------------------------------------------------------------*/
/*! \file
\brief implementation of the 3D solid-poro element

\maintainer Johannes Kremheller
            kremheller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249

\level 2

*----------------------------------------------------------------------*/

#include "so3_poro.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_utils_factory.H"

#include "so3_poro_eletypes.H"

#include "so_surface.H"
#include "so_line.H"

// for ReadElement()
#include "../drt_mat/structporo.H"
// for secondDerivativesZero
#include "../drt_fem_general/drt_utils_shapefunctions_service.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Poro<so3_ele, distype>::So3_Poro(int id, int owner)
    : so3_ele(id, owner),
      data_(),
      intpoints_(distype),
      init_(false),
      scatracoupling_(false),
      isNurbs_(false),
      weights_(true),
      myknots_(numdim_),
      fluidmat_(Teuchos::null),
      fluidmultimat_(Teuchos::null),
      structmat_(Teuchos::null)
// numscal_(0)
{
  numgpt_ = intpoints_.NumPoints();
  // ishigherorder_ = DRT::UTILS::secondDerivativesZero<distype>();

  invJ_.resize(numgpt_, LINALG::Matrix<numdim_, numdim_>(true));
  detJ_.resize(numgpt_, 0.0);
  xsi_.resize(numgpt_, LINALG::Matrix<numdim_, 1>(true));

  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Poro<so3_ele, distype>::So3_Poro(
    const DRT::ELEMENTS::So3_Poro<so3_ele, distype>& old)
    : so3_ele(old),
      data_(old.data_),
      invJ_(old.invJ_),
      detJ_(old.detJ_),
      xsi_(old.xsi_),
      intpoints_(distype),
      // ishigherorder_(old.ishigherorder_),
      init_(old.init_),
      scatracoupling_(old.scatracoupling_),
      isNurbs_(old.isNurbs_),
      weights_(old.weights_),
      myknots_(old.myknots_),
      fluidmat_(old.fluidmat_),
      fluidmultimat_(old.fluidmultimat_),
      structmat_(old.structmat_)
// numscal_(old.numscal_)
{
  numgpt_ = intpoints_.NumPoints();
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Clone() const
{
  DRT::ELEMENTS::So3_Poro<so3_ele, distype>* newelement =
      new DRT::ELEMENTS::So3_Poro<so3_ele, distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data, type);

  // data_
  so3_ele::AddtoPack(data, data_);

  // detJ_
  so3_ele::AddtoPack(data, detJ_);

  // invJ_
  int size = (int)invJ_.size();
  so3_ele::AddtoPack(data, size);
  for (int i = 0; i < size; ++i) so3_ele::AddtoPack(data, invJ_[i]);

  // xsi_
  size = (int)xsi_.size();
  so3_ele::AddtoPack(data, size);
  for (int i = 0; i < size; ++i) so3_ele::AddtoPack(data, xsi_[i]);

  // scatracoupling_
  so3_ele::AddtoPack(data, scatracoupling_);

  // isNurbs_
  so3_ele::AddtoPack(data, isNurbs_);

  // so3_ele::AddtoPack(data,numscal_);

  // add base class Element
  so3_ele::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  so3_ele::ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // data_
  std::vector<char> tmp(0);
  so3_ele::ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);

  // detJ_
  so3_ele::ExtractfromPack(position, data, detJ_);

  // invJ_
  int size = 0;
  so3_ele::ExtractfromPack(position, data, size);
  invJ_.resize(size, LINALG::Matrix<numdim_, numdim_>(true));
  for (int i = 0; i < size; ++i) so3_ele::ExtractfromPack(position, data, invJ_[i]);

  // xsi_
  size = 0;
  so3_ele::ExtractfromPack(position, data, size);
  xsi_.resize(size, LINALG::Matrix<numdim_, 1>(true));
  for (int i = 0; i < size; ++i) so3_ele::ExtractfromPack(position, data, xsi_[i]);

  // scatracoupling_
  scatracoupling_ = (bool)(so3_ele::ExtractInt(position, data));

  // isNurbs_
  isNurbs_ = (bool)(so3_ele::ExtractInt(position, data));

  // numscal_ = so3_ele::ExtractInt(position,data);

  // extract base class Element
  std::vector<char> basedata(0);
  so3_ele::ExtractfromPack(position, data, basedata);
  so3_ele::Unpack(basedata);

  init_ = true;

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                 vuong 11/13|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Volumes()
{
  std::vector<Teuchos::RCP<Element>> volumes(1);
  volumes[0] = Teuchos::rcp(this, false);
  return volumes;
}

/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                           vuong 11/13|
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface, DRT::Element>(
      DRT::UTILS::buildSurfaces, this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             vuong 11/13|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine, DRT::Element>(
      DRT::UTILS::buildLines, this);
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Print(std::ostream& os) const
{
  os << "So3_poro ";
  os << DRT::DistypeToString(distype).c_str() << " ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, DRT::INPUT::LineDefinition* linedef)
{
  // read base element
  so3_ele::ReadElement(eletype, eledistype, linedef);

  // setup poro material
  Teuchos::RCP<MAT::StructPoro> poromat = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(Material());
  if (poromat == Teuchos::null) dserror("no poro material assigned to poro element!");
  poromat->PoroSetup(numgpt_, linedef);

  return true;
}

/*----------------------------------------------------------------------*
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::VisNames(std::map<std::string, int>& names)
{
  so3_ele::VisNames(names);
}

/*----------------------------------------------------------------------*
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::So3_Poro<so3_ele, distype>::VisData(
    const std::string& name, std::vector<double>& data)
{
  return so3_ele::VisData(name, data);
}

/*----------------------------------------------------------------------*
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Poro<so3_ele, distype>::UniqueParObjectId() const
{
  switch (distype)
  {
    case DRT::Element::tet4:
      return So_tet4PoroType::Instance().UniqueParObjectId();
      break;
    case DRT::Element::tet10:
      return So_tet10PoroType::Instance().UniqueParObjectId();
      break;
    case DRT::Element::hex8:
      return So_hex8PoroType::Instance().UniqueParObjectId();
      break;
    case DRT::Element::hex27:
      return So_hex27PoroType::Instance().UniqueParObjectId();
      break;
    case DRT::Element::nurbs27:
      return So_nurbs27PoroType::Instance().UniqueParObjectId();
      break;
    default:
      dserror("unknown element type!");
      break;
  }
  return -1;
}

/*----------------------------------------------------------------------*
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ElementType& DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ElementType() const
{
  switch (distype)
  {
    case DRT::Element::tet4:
      return So_tet4PoroType::Instance();
    case DRT::Element::tet10:
      return So_tet10PoroType::Instance();
    case DRT::Element::hex8:
      return So_hex8PoroType::Instance();
    case DRT::Element::hex27:
      return So_hex27PoroType::Instance();
    case DRT::Element::nurbs27:
      return So_nurbs27PoroType::Instance();
    default:
      dserror("unknown element type!");
      break;
  }
  return So_hex8PoroType::Instance();
};

/*----------------------------------------------------------------------*
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
inline DRT::Node** DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Nodes()
{
  return so3_ele::Nodes();
}

/*----------------------------------------------------------------------*
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
inline Teuchos::RCP<MAT::Material> DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Material() const
{
  return so3_ele::Material();
}

/*----------------------------------------------------------------------*
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
inline int DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Id() const
{
  return so3_ele::Id();
}

#include "so3_poro_fwd.hpp"
