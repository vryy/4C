/*---------------------------------------------------------------------*/
/*! \file

\brief A data storage container

\level 0

*/
/*---------------------------------------------------------------------*/

#include "baci_lib_container.hpp"

#include "baci_linalg_serialdensematrix.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_Vector.h>

BACI_NAMESPACE_OPEN


DRT::ContainerType DRT::ContainerType::instance_;


CORE::COMM::ParObject* DRT::ContainerType::Create(const std::vector<char>& data)
{
  DRT::Container* object = new DRT::Container();
  object->Unpack(data);
  return object;
}


std::ostream& operator<<(std::ostream& os, const DRT::Container& cont)
{
  cont.Print(os);
  return os;
}


void DRT::Container::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  const auto pack_map = [](auto& data, const auto& map)
  {
    AddtoPack(data, static_cast<int>(map.size()));
    for (const auto& [key, value] : map)
    {
      AddtoPack(data, key);
      AddtoPack(data, value);
    }
  };

  // pack type of this instance of ParObject
  AddtoPack(data, UniqueParObjectId());

  pack_map(data, intdata_);
  pack_map(data, doubledata_);
  pack_map(data, mapdata_);
  pack_map(data, stringdata_);
  pack_map(data, matdata_);
}


void DRT::Container::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  const auto unpack_map = [this](auto& position, const auto& data, auto& target_map)
  {
    int target_size;
    ExtractfromPack(position, data, target_size);

    for (int i = 0; i < target_size; ++i)
    {
      std::string key;
      ExtractfromPack(position, data, key);
      typename std::decay_t<decltype(target_map)>::mapped_type value;
      ExtractfromPack(position, data, value);
      Add(key, value);
    }
  };

  unpack_map(position, data, intdata_);
  unpack_map(position, data, doubledata_);
  unpack_map(position, data, mapdata_);
  unpack_map(position, data, stringdata_);
  unpack_map(position, data, matdata_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
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


void DRT::Container::Print(std::ostream& os) const
{
  PrintHelper printer{os};
  printer(intdata_);
  printer(doubledata_);
  printer(mapdata_);
  printer(stringdata_);
  printer(matdata_);
}


void DRT::Container::Add(const std::string& name, const int* data, const std::size_t num)
{
  intdata_[name] = std::vector<int>(data, data + num);
}

void DRT::Container::Add(const std::string& name, const double* data, const std::size_t num)
{
  doubledata_[name] = std::vector<double>(data, data + num);
}

void DRT::Container::Add(const std::string& name, const std::map<int, std::vector<double>>& data)
{
  mapdata_[name] = data;
}

void DRT::Container::Add(const std::string& name, const std::string& data)
{
  stringdata_[name] = data;
}

void DRT::Container::Add(const std::string& name, const CORE::LINALG::SerialDenseMatrix& matrix)
{
  matdata_[name] = matrix;
}

void DRT::Container::Add(const std::string& name, const std::vector<int>& data)
{
  Add(name, data.data(), data.size());
}

void DRT::Container::Add(const std::string& name, const std::vector<double>& data)
{
  Add(name, data.data(), data.size());
}

void DRT::Container::Add(const std::string& name, const int& data) { Add(name, &data, 1); }

void DRT::Container::Add(const std::string& name, const double& data) { Add(name, &data, 1); }


void DRT::Container::Delete(const std::string& name)
{
  [[maybe_unused]] const unsigned n_erased_entries =
      intdata_.erase(name) + doubledata_.erase(name) + mapdata_.erase(name) +
      stringdata_.erase(name) + mapdata_.erase(name);
  dsassert(n_erased_entries == 1, "Internal error: key was present in more than one map.");
}


int DRT::Container::GetInt(const std::string& name) const
{
  const std::vector<int>* vecptr = Get<std::vector<int>>(name);
  if (vecptr == nullptr) dserror("Integer %s cannot be read from the container.", name.c_str());
  if (vecptr->size() != 1) dserror("Trying to read integer from vector of wrong length.");
  return (*vecptr)[0];
}


double DRT::Container::GetDouble(const std::string& name) const
{
  const std::vector<double>* vecptr = Get<std::vector<double>>(name);
  if (vecptr == nullptr) dserror("Double %s cannot be read from the container.", name.c_str());
  if (vecptr->size() != 1) dserror("Trying to read double from vector of wrong length.");
  return (*vecptr)[0];
}


int DRT::Container::UniqueParObjectId() const
{
  return ContainerType::Instance().UniqueParObjectId();
}


void DRT::Container::Clear()
{
  intdata_.clear();
  doubledata_.clear();
  mapdata_.clear();
  stringdata_.clear();
  matdata_.clear();
}

BACI_NAMESPACE_CLOSE
