/*----------------------------------------------------------------------*/
/*! \file
\brief A material parameter container

\level 1

*/
/*---------------------------------------------------------------------*/


#include "4C_material_input_base.hpp"

FOUR_C_NAMESPACE_OPEN


CORE::MAT::PAR::Material::Material(
    const int id, const CORE::Materials::MaterialType type, const std::string name)
    : InputParameterContainer(), id_(id), type_(type), name_(name), params_(Teuchos::null)
{
}

CORE::MAT::PAR::Material::Material(const CORE::MAT::PAR::Material& old)
    : InputParameterContainer(old), id_(old.id_), type_(old.type_), params_(old.params_)
{
}

std::ostream& operator<<(std::ostream& os, const CORE::MAT::PAR::Material& cond)
{
  cond.Print(os);
  return os;
}

void CORE::MAT::PAR::Material::Print(std::ostream& os) const
{
  os << "MAT " << Id() << " " << Name() << " :: ";

  IO::InputParameterContainer::Print(os);
}

FOUR_C_NAMESPACE_CLOSE
