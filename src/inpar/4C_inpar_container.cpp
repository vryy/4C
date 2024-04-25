/*---------------------------------------------------------------------*/
/*! \file

\brief A data storage container

\level 0

*/
/*---------------------------------------------------------------------*/

#include "4C_inpar_container.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN


std::ostream& operator<<(std::ostream& os, const INPAR::InputParameterContainer& cont)
{
  cont.Print(os);
  return os;
}

namespace
{
  //! Print various types that occurr in the Container
  struct PrintHelper
  {
    //! Base case: print the object directly.
    template <typename T>
    void operator()(const T& object)
    {
      os << object << " ";
    }

    //! Print elements of a vector.
    template <typename T>
    void operator()(const std::vector<T>& vector)
    {
      for (const auto& v : vector)
      {
        (*this)(v);
      }
    }

    //! Print elements of a map.
    template <typename Key, typename Value>
    void operator()(const std::map<Key, Value>& map)
    {
      for (const auto& [key, value] : map)
      {
        os << key << " : ";
        (*this)(value);
      }
    }

    //! Print any data.
    void operator()(const std::any& /*unused*/) { os << "non-printable data of type std::any"; }

    std::ostream& os;
  };
}  // namespace


void INPAR::InputParameterContainer::Print(std::ostream& os) const
{
  PrintHelper printer{os};
  printer(intdata_);
  printer(doubledata_);
  printer(booldata_);
  printer(vecintdata_);
  printer(vecdoubledata_);
  printer(mapdata_);
  printer(stringdata_);
  printer(matdata_);
  printer(anydata_);
}


FOUR_C_NAMESPACE_CLOSE
