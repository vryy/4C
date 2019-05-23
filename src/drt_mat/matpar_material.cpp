/*!-----------------------------------------------------------------------
\brief A condition of any kind

\level 1

\maintainer Fabian Braeu

*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_material.H"


MAT::PAR::ParMaterialType MAT::PAR::ParMaterialType::instance_;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Material::Material(
    const int id, const INPAR::MAT::MaterialType type, const std::string name)
    : Container(), id_(id), type_(type), name_(name), comm_(Teuchos::null), params_(Teuchos::null)
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Material::Material()
    : Container(),
      id_(-1),
      type_(INPAR::MAT::m_none),
      name_(""),
      comm_(Teuchos::null),
      params_(Teuchos::null)
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Material::Material(const MAT::PAR::Material& old)
    : Container(old), id_(old.id_), type_(old.type_), comm_(old.comm_), params_(old.params_)
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Material::~Material() { return; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const MAT::PAR::Material& cond)
{
  cond.Print(os);
  return os;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PAR::Material::Print(std::ostream& os) const
{
  os << "MAT " << Id() << " " << Name() << " :: ";

  DRT::Container::Print(os);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PAR::Material::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class container
  DRT::Container::Pack(data);
  // id_
  AddtoPack(data, id_);
  // type_
  AddtoPack(data, type_);
  // name_
  AddtoPack(data, name_);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PAR::Material::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Container
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Container::Unpack(basedata);
  // id_
  ExtractfromPack(position, data, id_);
  // type_
  type_ = static_cast<INPAR::MAT::MaterialType>(ExtractInt(position, data));
  // name_
  ExtractfromPack(position, data, name_);

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
  return;
}


/*----------------------------------------------------------------------*/
