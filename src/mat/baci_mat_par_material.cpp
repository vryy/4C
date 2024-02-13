/*----------------------------------------------------------------------*/
/*! \file
\brief A material parameter container

\level 1

*/
/*---------------------------------------------------------------------*/


#include "baci_mat_par_material.hpp"

BACI_NAMESPACE_OPEN


MAT::PAR::Material::Material(
    const int id, const INPAR::MAT::MaterialType type, const std::string name)
    : InputParameterContainer(), id_(id), type_(type), name_(name), params_(Teuchos::null)
{
}

MAT::PAR::Material::Material(const MAT::PAR::Material& old)
    : InputParameterContainer(old), id_(old.id_), type_(old.type_), params_(old.params_)
{
}

std::ostream& operator<<(std::ostream& os, const MAT::PAR::Material& cond)
{
  cond.Print(os);
  return os;
}

void MAT::PAR::Material::Print(std::ostream& os) const
{
  os << "MAT " << Id() << " " << Name() << " :: ";

  INPAR::InputParameterContainer::Print(os);
}

BACI_NAMESPACE_CLOSE
