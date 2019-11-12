/*----------------------------------------------------------------------*/
/*! \file

\brief Implements a container class containing contact constitutivelaw parameters from the input
file as well as a quick access wrapper for those parameters

\level 1

\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/

#include "contactconstitutivelaw_parameter.H"

#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_dserror.H"


CONTACT::CONSTITUTIVELAW::Parameter::Parameter(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container>
        coconstlawdata  ///< read and validate contactconstitutivelaw data (of 'slow' access)
    )
    : id_(coconstlawdata->Id()),
      offset_(coconstlawdata->GetDouble("Offset")),
      type_(coconstlawdata->Type()),
      name_(coconstlawdata->Name()){};
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::ConstitutiveLawType
    CONTACT::CONSTITUTIVELAW::ConstitutiveLawType::instance_;
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::Container::Container(
    const int id, const INPAR::CONTACT::ConstitutiveLawType type, const std::string name)
    : DRT::Container(), id_(id), type_(type), name_(name), params_(Teuchos::null)
{
  return;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::Container::Container()
    : DRT::Container(),
      id_(-1),
      type_(INPAR::CONTACT::ConstitutiveLawType::colaw_none),
      name_(""),
      params_(Teuchos::null)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::Container::Container(const CONTACT::CONSTITUTIVELAW::Container& old)
    : DRT::Container(old), id_(old.id_), type_(old.type_), params_(old.params_)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CONTACT::CONSTITUTIVELAW::Container& cond)
{
  cond.Print(os);
  return os;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::Container::Print(std::ostream& os) const
{
  os << "ContactConstitutiveLaw " << Id() << " " << Name() << " :: ";

  DRT::Container::Print(os);

  return;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::Container::Pack(DRT::PackBuffer& data) const
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
  AddtoPack(data, static_cast<int>(type_));
  // name_
  AddtoPack(data, name_);

  return;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::Container::Unpack(const std::vector<char>& data)
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
  type_ = static_cast<INPAR::CONTACT::ConstitutiveLawType>(ExtractInt(position, data));
  // name_
  ExtractfromPack(position, data, name_);

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
  return;
}
