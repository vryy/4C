/*---------------------------------------------------------------------*/
/*! \file

\brief A data storage container

\level 0

*/
/*---------------------------------------------------------------------*/

#ifndef BACI_INPAR_CONTAINER_HPP
#define BACI_INPAR_CONTAINER_HPP


#include "baci_config.hpp"

#include "baci_linalg_serialdensematrix.hpp"
#include "baci_utils_demangle.hpp"
#include "baci_utils_exceptions.hpp"

#include <Teuchos_ENull.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <any>
#include <string>

class Epetra_MultiVector;
class Epetra_Vector;

BACI_NAMESPACE_OPEN

namespace INPAR
{
  namespace INTERNAL
  {
    template <typename T>
    const T* TryGetAnyData(const std::string& name, const std::any& data)
    {
      if (typeid(T) == data.type())
      {
        const T* any_ptr = std::any_cast<T>(&data);
        dsassert(any_ptr != nullptr, "Implementation error.");
        return any_ptr;
      }
      else
      {
        dserror(
            "You tried to get the data named %s from the container as type '%s'.\n"
            "Actually, it has type '%s'.",
            name.c_str(), CORE::UTILS::TryDemangle(typeid(T).name()).c_str(),
            CORE::UTILS::TryDemangle(data.type().name()).c_str());
      }
    }
  }  // namespace INTERNAL

  /*!
  \brief A data storage container

  You can store various types of data in this container and access the data by keys.

  The intention of this class is to store rather 'small' units of data. Though possible,
  it is not meant to be used at the system level to store huge data sets such as sparse matrices or
  vectors of system length. It does therefore not support any Epetra_Vector
  or Epetra_CrsMatrix objects and is not supposed to in the future either.
  */
  class InputParameterContainer
  {
   public:
    /*!
     * \brief Default constructor.
     */
    InputParameterContainer() = default;

    /*!
     * \brief Destructor.
     */
    virtual ~InputParameterContainer() = default;

    /*!
     * \brief Print this container.
     */
    virtual void Print(std::ostream& os) const;


    /*!
     * \brief Add @data to the container at the given key @name.
     *
     * If a record with name name already exists, it will be overwritten.
     */
    template <typename T>
    void Add(const std::string& name, const T& data)
    {
      if constexpr (std::is_same_v<T, int>)
        intdata_[name] = data;
      else if constexpr (std::is_same_v<T, double>)
        doubledata_[name] = data;
      else if constexpr (std::is_same_v<T, bool>)
        booldata_[name] = data;
      else if constexpr (std::is_same_v<T, std::vector<int>>)
        vecintdata_[name] = data;
      else if constexpr (std::is_same_v<T, std::vector<double>>)
        vecdoubledata_[name] = data;
      else if constexpr (std::is_same_v<T, std::map<int, std::vector<double>>>)
        mapdata_[name] = data;
      else if constexpr (std::is_same_v<T, std::string>)
        stringdata_[name] = data;
      else if constexpr (std::is_same_v<T, CORE::LINALG::SerialDenseMatrix>)
        matdata_[name] = data;
      else
        anydata_[name] = data;
    }


    /*!
     * \brief Get the const pointer to the data stored at the key @name from the container. An error
     * is thrown in case the key @name is not existent in the specific data map.
     */
    template <typename T>
    const T* Get(const std::string& name) const
    {
      // access function for known types
      [[maybe_unused]] const auto access = [](const auto& map, const std::string& name) -> const T*
      {
        const auto it = map.find(name);
        if (it != map.end())
        {
          return std::addressof(it->second);
        }
        else
        {
          dserror("Key %s cannot be found in the container's map!", name.c_str());
        }
      };

      // access function for any type
      [[maybe_unused]] const auto access_any = [](const auto& map,
                                                   const std::string& name) -> const T*
      {
        const auto it = map.find(name);
        if (it != map.end())
        {
          TryGetAnyData(name, it->second);
        }
        else
        {
          dserror("Key %s cannot be found in the container's map!", name.c_str());
        }
      };

      // get data from the container
      if constexpr (std::is_same_v<T, int>)
        return access(intdata_, name);
      else if constexpr (std::is_same_v<T, double>)
        return access(doubledata_, name);
      else if constexpr (std::is_same_v<T, bool>)
        return access(booldata_, name);
      else if constexpr (std::is_same_v<T, std::vector<int>>)
        return access(vecintdata_, name);
      else if constexpr (std::is_same_v<T, std::vector<double>>)
        return access(vecdoubledata_, name);
      else if constexpr (std::is_same_v<T, std::map<int, std::vector<double>>>)
        return access(mapdata_, name);
      else if constexpr (std::is_same_v<T, std::string>)
        return access(stringdata_, name);
      else if constexpr (std::is_same_v<T, CORE::LINALG::SerialDenseMatrix>)
        return access(matdata_, name);
      else
        return access_any(anydata_, name);
    };

    /*!
     * \brief Get the pointer to the data stored at the key @name from the container. An error is
     * thrown in case the key @name is not existent in the specific data map.
     */
    template <typename T>
    T* Get(const std::string& name)
    {
      return const_cast<T*>(const_cast<const InputParameterContainer*>(this)->Get<T>(name));
    };


    /*!
     * \brief Get the const pointer to the data stored at the key @name from the container. Return
     * the @default_value if the @name is not existent in the specific data map.
     */
    template <typename T>
    const T* GetOr(const std::string& name, T* default_value) const
    {
      // access function for known types
      [[maybe_unused]] const auto access = [](const auto& map, const std::string& name,
                                               T* default_value) -> const T*
      {
        const auto it = map.find(name);
        if (it != map.end())
        {
          return std::addressof(it->second);
        }
        else
        {
          return default_value;
        }
      };

      // access function for any type
      [[maybe_unused]] const auto access_any = [](const auto& map, const std::string& name,
                                                   T* default_value) -> const T*
      {
        const auto it = map.find(name);
        if (it != map.end())
        {
          TryGetAnyData(name, it->second);
        }
        else
        {
          return default_value;
        }
      };

      // get data from the container
      if constexpr (std::is_same_v<T, int>)
        return access(intdata_, name, default_value);
      else if constexpr (std::is_same_v<T, double>)
        return access(doubledata_, name, default_value);
      else if constexpr (std::is_same_v<T, bool>)
        return access(booldata_, name, default_value);
      else if constexpr (std::is_same_v<T, std::vector<int>>)
        return access(vecintdata_, name, default_value);
      else if constexpr (std::is_same_v<T, std::vector<double>>)
        return access(vecdoubledata_, name, default_value);
      else if constexpr (std::is_same_v<T, std::map<int, std::vector<double>>>)
        return access(mapdata_, name, default_value);
      else if constexpr (std::is_same_v<T, std::string>)
        return access(stringdata_, name, default_value);
      else if constexpr (std::is_same_v<T, CORE::LINALG::SerialDenseMatrix>)
        return access(matdata_, name, default_value);
      else
        return access_any(anydata_, name, default_value);
    };

    /*!
     * \brief Get the pointer to the data stored at the key @name from the container. Return
     * the @default_value if the @name is not existent in the specific data map.
     */
    template <typename T>
    T* GetOr(const std::string& name, T* default_value)
    {
      return const_cast<T*>(
          const_cast<const InputParameterContainer*>(this)->GetOr<T>(name, default_value));
    };


    /*!
     * \brief Get the const pointer to the data stored at the key @name from the container. Return a
     * nullpointer if the @name is not existent in the specific data map.
     */
    template <typename T>
    const T* GetIf(const std::string& name) const
    {
      return GetOr<T>(name, nullptr);
    };

    /*!
     * \brief Get the pointer to the data stored at the key @name from the container. Return a
     * nullpointer if the @name is not existent in the specific data map.
     */
    template <typename T>
    T* GetIf(const std::string& name)
    {
      return const_cast<T*>(const_cast<const InputParameterContainer*>(this)->GetIf<T>(name));
    };


   private:
    //! a map to store integer data in
    std::map<std::string, int> intdata_;

    //! a map to store double data in
    std::map<std::string, double> doubledata_;

    //! a map to store bool data in
    std::map<std::string, bool> booldata_;

    //! a map to store vector integer data in
    std::map<std::string, std::vector<int>> vecintdata_;

    //! a map to store vector double data in
    std::map<std::string, std::vector<double>> vecdoubledata_;

    //! a map to store maps of stl vector data
    std::map<std::string, std::map<int, std::vector<double>>> mapdata_;

    //! a map to store string data in
    std::map<std::string, std::string> stringdata_;

    //! a map to store matrices in
    std::map<std::string, CORE::LINALG::SerialDenseMatrix> matdata_;

    //! a map to store anything
    std::map<std::string, std::any> anydata_;
  };  // class InputParameterContainer
}  // namespace INPAR

// << operator
std::ostream& operator<<(std::ostream& os, const INPAR::InputParameterContainer& cont);


BACI_NAMESPACE_CLOSE

#endif  // INPAR_CONTAINER_H
