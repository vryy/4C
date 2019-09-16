/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#include "drt_utils_factory.H"
#include "drt_parobjectfactory.H"

/*----------------------------------------------------------------------*
 |  allocate an instance of a specific impl. of ParObject (public) mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::UTILS::Factory(const std::vector<char>& data)
{
  return ParObjectFactory::Instance().Create(data);
}

/*----------------------------------------------------------------------*
 |  allocate an element of a specific type (public)          mwgee 03|07|
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::UTILS::Factory(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  return ParObjectFactory::Instance().Create(eletype, eledistype, id, owner);
}
