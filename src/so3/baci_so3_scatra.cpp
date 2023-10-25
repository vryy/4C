/*----------------------------------------------------------------------*/
/*! \file

\brief Solid-scatra elements base class

\level 2


*----------------------------------------------------------------------*/


#include "baci_so3_scatra.H"

#include "baci_inpar_ssi.H"
#include "baci_io_linedefinition.H"
#include "baci_lib_element_integration_select.H"
#include "baci_lib_globalproblem.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::So3_Scatra(int id, int owner)
    : so3_ele(id, owner),
      impltype_(INPAR::SCATRA::impltype_undefined),
      intpoints_(distype == CORE::FE::CellType::tet4
                     ? CORE::DRT::UTILS::GaussRule3D::tet_1point
                     : (distype == CORE::FE::CellType::tet10
                               ? CORE::DRT::UTILS::GaussRule3D::tet_4point
                               : DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule)),
      numgpt_(intpoints_.nquad),
      xsi_(0.0),
      invJ_(0.0),
      detJ_(0.0)
{
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::So3_Scatra(
    const DRT::ELEMENTS::So3_Scatra<so3_ele, distype>& old)
    : so3_ele(old),
      impltype_(old.impltype_),
      intpoints_(old.intpoints_),
      numgpt_(old.numgpt_),
      xsi_(old.xsi_),
      invJ_(old.invJ_),
      detJ_(old.detJ_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
DRT::Element* DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::So3_Scatra<so3_ele, distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data, type);

  // pack scalar transport impltype
  so3_ele::AddtoPack(data, impltype_);
  // data_
  // so3_ele::AddtoPack(data,data_);

  // detJ_
  so3_ele::AddtoPack(data, detJ_);

  // invJ_
  auto size = (int)invJ_.size();
  so3_ele::AddtoPack(data, size);
  for (int i = 0; i < size; ++i) so3_ele::AddtoPack(data, invJ_[i]);

  // xsi_
  size = (int)xsi_.size();
  so3_ele::AddtoPack(data, size);
  for (int i = 0; i < size; ++i) so3_ele::AddtoPack(data, xsi_[i]);


  // add base class Element
  so3_ele::Pack(data);



  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  so3_ele::ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract scalar transport impltype_
  impltype_ = static_cast<INPAR::SCATRA::ImplType>(so3_ele::ExtractInt(position, data));

  // data_
  // vector<char> tmp(0);
  // so3_ele::ExtractfromPack(position,data,tmp);
  // data_.Unpack(tmp);

  // detJ_
  so3_ele::ExtractfromPack(position, data, detJ_);

  // invJ_
  int size = 0;
  so3_ele::ExtractfromPack(position, data, size);
  invJ_.resize(size, CORE::LINALG::Matrix<numdim_, numdim_>(true));
  for (int i = 0; i < size; ++i) so3_ele::ExtractfromPack(position, data, invJ_[i]);

  // xsi_
  size = 0;
  so3_ele::ExtractfromPack(position, data, size);
  xsi_.resize(size, CORE::LINALG::Matrix<numdim_, 1>(true));
  for (int i = 0; i < size; ++i) so3_ele::ExtractfromPack(position, data, xsi_[i]);

  // extract base class Element
  std::vector<char> basedata(0);
  so3_ele::ExtractfromPack(position, data, basedata);

  so3_ele::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::Print(std::ostream& os) const
{
  os << "So3_scatra ";
  os << " Discretization type: " << DRT::DistypeToString(distype).c_str();
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             schmidt 09/17|
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
bool DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, DRT::INPUT::LineDefinition* linedef)
{
  so3_ele::ReadElement(eletype, eledistype, linedef);

  // read scalar transport implementation type
  std::string impltype;
  linedef->ExtractString("TYPE", impltype);

  if (impltype == "Undefined")
    impltype_ = INPAR::SCATRA::impltype_undefined;
  else if (impltype == "AdvReac")
    impltype_ = INPAR::SCATRA::impltype_advreac;
  else if (impltype == "CardMono")
    impltype_ = INPAR::SCATRA::impltype_cardiac_monodomain;
  else if (impltype == "Chemo")
    impltype_ = INPAR::SCATRA::impltype_chemo;
  else if (impltype == "ChemoReac")
    impltype_ = INPAR::SCATRA::impltype_chemoreac;
  else if (impltype == "ElchDiffCond")
    impltype_ = INPAR::SCATRA::impltype_elch_diffcond;
  else if (impltype == "ElchElectrode")
    impltype_ = INPAR::SCATRA::impltype_elch_electrode;
  else if (impltype == "Loma")
    impltype_ = INPAR::SCATRA::impltype_loma;
  else if (impltype == "RefConcReac")
    impltype_ = INPAR::SCATRA::impltype_refconcreac;
  else if (impltype == "Std")
    impltype_ = INPAR::SCATRA::impltype_std;
  else
    dserror("Invalid implementation type for So3_Scatra elements!");

  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
inline DRT::Node** DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::Nodes()
{
  return so3_ele::Nodes();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
inline Teuchos::RCP<MAT::Material> DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::Material() const
{
  return so3_ele::Material();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
inline int DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::Id() const
{
  return so3_ele::Id();
}


/*--------------------------------------------------------------------------*
 | set the material  (public)                                 schmidt 10/17 |
 *                                                                          */
template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::SetMaterial(int matnum)
{
  // call base class
  so3_ele::SetMaterial(matnum);

  // get the scatra structure control parameter list
  const Teuchos::ParameterList& ssicontrol = DRT::Problem::Instance()->SSIControlParams();
  // get the material
  Teuchos::RCP<MAT::Material> mat = Material();

  if ((Teuchos::getIntegralValue<INPAR::SSI::SolutionSchemeOverFields>(ssicontrol, "COUPALGO") ==
          INPAR::SSI::SolutionSchemeOverFields::ssi_Monolithic) and
      (mat->MaterialType() != INPAR::MAT::m_multiplicative_split_defgrad_elasthyper))
    dserror(
        "When you use the 'COUPALGO' 'ssi_Monolithic' from the 'SSI CONTROL' section,"
        " you need to use the material 'MAT_MultiplicativeSplitDefgradElastHyper'! If you want to "
        "use another material feel free to implement it! ;-)");

  return;
}


template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar, CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27, CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10, CORE::FE::CellType::tet10>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, CORE::FE::CellType::wedge6>;
