/*----------------------------------------------------------------------*/
/*!
 \file so_base.cpp

\maintainer Anh-Tu Vuong

 *----------------------------------------------------------------------*/


#include "so_base.H"
#include "../drt_mat/so3_material.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_base::So_base(int id, int owner) :
DRT::Element(id,owner),
kintype_(INPAR::STR::kinem_vague),
interface_ptr_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_base::So_base(const DRT::ELEMENTS::So_base& old) :
DRT::Element(old),
kintype_(old.kintype_),
interface_ptr_(old.interface_ptr_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           vuong 03/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_base::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);
  // kintype_
  AddtoPack(data,kintype_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_base::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // kintype_
  kintype_ = static_cast<INPAR::STR::KinemType>( ExtractInt(position,data) );

  return;
}

/*----------------------------------------------------------------------*
 |  return solid material                                      (public) |
 |                                                           seitz 03/15|
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::So3Material> DRT::ELEMENTS::So_base::SolidMaterial(int nummat) const
{
  return Teuchos::rcp_dynamic_cast<MAT::So3Material>(DRT::Element::Material(nummat),true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_base::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ =
        Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>
        (p.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> >("interface"));
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> DRT::ELEMENTS::So_base::ParamsInterfacePtr()
{
  return interface_ptr_;
}
