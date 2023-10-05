/*----------------------------------------------------------------------*/
/*! \file
\brief Base class for components of input lines
\level 0
*/
/*----------------------------------------------------------------------*/

#include "baci_lib_linecomponent.H"

#include <iterator>

namespace INPUT
{

  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  LineComponent::LineComponent(std::string name, bool optional)
      : optional_(optional), name_(std::move(name))
  {
  }


  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  Teuchos::RCP<std::stringstream> LineComponent::PushBack(
      const std::string& token, const Teuchos::RCP<std::stringstream>& stream)
  {
    Teuchos::RCP<std::stringstream> out = Teuchos::rcp(new std::stringstream());
    (*out) << token << " ";
    std::copy(std::istream_iterator<std::string>(*stream), std::istream_iterator<std::string>(),
        std::ostream_iterator<std::string>(*out, " "));
    return out;
  }
}  // namespace INPUT
