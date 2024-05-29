/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0


*/
/*---------------------------------------------------------------------*/

#include "4C_comm_utils_factory.hpp"

#include "4C_comm_parobjectfactory.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  allocate an instance of a specific impl. of ParObject (public) mwgee 12/06|
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* CORE::COMM::Factory(const std::vector<char>& data)
{
  return ParObjectFactory::Instance().Create(data);
}

/*----------------------------------------------------------------------*
 |  allocate an element of a specific type (public)          mwgee 03|07|
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::Elements::Element> CORE::COMM::Factory(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  return ParObjectFactory::Instance().Create(eletype, eledistype, id, owner);
}

FOUR_C_NAMESPACE_CLOSE
