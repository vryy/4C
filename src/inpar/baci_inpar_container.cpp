/*---------------------------------------------------------------------*/
/*! \file

\brief A data storage container

\level 0

*/
/*---------------------------------------------------------------------*/

#include "baci_inpar_container.hpp"

#include "baci_linalg_serialdensematrix.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_Vector.h>

BACI_NAMESPACE_OPEN


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
    //! Base case: print the object directly
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

    std::ostream& os;
  };
}  // namespace


void INPAR::InputParameterContainer::Print(std::ostream& os) const
{
  PrintHelper printer{os};
  printer(intdata_);
  printer(doubledata_);
  printer(mapdata_);
  printer(stringdata_);
  printer(matdata_);
}


void INPAR::InputParameterContainer::Add(
    const std::string& name, const int* data, const std::size_t num)
{
  intdata_[name] = std::vector<int>(data, data + num);
}

void INPAR::InputParameterContainer::Add(
    const std::string& name, const double* data, const std::size_t num)
{
  doubledata_[name] = std::vector<double>(data, data + num);
}

void INPAR::InputParameterContainer::Add(
    const std::string& name, const std::map<int, std::vector<double>>& data)
{
  mapdata_[name] = data;
}

void INPAR::InputParameterContainer::Add(const std::string& name, const std::string& data)
{
  stringdata_[name] = data;
}

void INPAR::InputParameterContainer::Add(
    const std::string& name, const CORE::LINALG::SerialDenseMatrix& matrix)
{
  matdata_[name] = matrix;
}

void INPAR::InputParameterContainer::Add(const std::string& name, const std::vector<int>& data)
{
  Add(name, data.data(), data.size());
}

void INPAR::InputParameterContainer::Add(const std::string& name, const std::vector<double>& data)
{
  Add(name, data.data(), data.size());
}

void INPAR::InputParameterContainer::Add(const std::string& name, const int& data)
{
  Add(name, &data, 1);
}

void INPAR::InputParameterContainer::Add(const std::string& name, const double& data)
{
  Add(name, &data, 1);
}

int INPAR::InputParameterContainer::GetInt(const std::string& name) const
{
  const std::vector<int>* vecptr = Get<std::vector<int>>(name);
  if (vecptr == nullptr) dserror("Integer %s cannot be read from the container.", name.c_str());
  if (vecptr->size() != 1) dserror("Trying to read integer from vector of wrong length.");
  return (*vecptr)[0];
}


double INPAR::InputParameterContainer::GetDouble(const std::string& name) const
{
  const std::vector<double>* vecptr = Get<std::vector<double>>(name);
  if (vecptr == nullptr) dserror("Double %s cannot be read from the container.", name.c_str());
  if (vecptr->size() != 1) dserror("Trying to read double from vector of wrong length.");
  return (*vecptr)[0];
}

BACI_NAMESPACE_CLOSE
