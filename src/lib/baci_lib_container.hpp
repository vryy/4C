/*---------------------------------------------------------------------*/
/*! \file

\brief A data storage container

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LIB_CONTAINER_HPP
#define BACI_LIB_CONTAINER_HPP


#include "baci_config.hpp"

#include "baci_comm_parobject.hpp"
#include "baci_comm_parobjectfactory.hpp"
#include "baci_linalg_serialdensematrix.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <memory>
#include <string>

class Epetra_MultiVector;
class Epetra_Vector;

BACI_NAMESPACE_OPEN

namespace DRT
{
  class ContainerType : public CORE::COMM::ParObjectType
  {
   public:
    [[nodiscard]] std::string Name() const override { return "ContainerType"; }

    static ContainerType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ContainerType instance_;
  };

  /*!
  \brief A data storage container

  You can store arrays of integer, double, and strings in this container
  and access the data by keys. As the Container implements ParObject
  everything stored in here can easily be communicated using the CORE::COMM::Exporter.

  The intention of this class is to store rather 'small' units of data. Though possible,
  it is not meant to be used at the system level to store huge data sets such as sparse matrices or
  vectors of system length. It does therefore not support any Epetra_Vector
  or Epetra_CrsMatrix objects and is not supposed to in the future either.

  */
  class Container : public CORE::COMM::ParObject
  {
   public:
    [[nodiscard]] int UniqueParObjectId() const override;

    void Pack(CORE::COMM::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    //! @name Query methods

    /*!
    \brief Print this element
    */
    virtual void Print(std::ostream& os) const;

    //@}

    //! @name Construction methods

    /*!
    \brief Add an integer to the container

    \param name : Name of data to store data with
    \param data : the integer to add

    \note If record with name name already exists, it will be overwritten
    */
    void Add(const std::string& name, const int& data);

    /*!
    \brief Add vector of int to the container

    \param name : Name of data to store data with
    \param data : ptr to beginning of data
    \param num  : Number of objects in data

    \note If record with name name already exists, it will be overwritten
    */
    void Add(const std::string& name, const int* data, const std::size_t num);

    /*!
    \brief Add vector of int to the container

    \param name : Name of data to store data with
    \param data : vector of integers to be added

    \note If record with name name already exists, it will be overwritten
    */
    void Add(const std::string& name, const std::vector<int>& data);

    /*!
    \brief Add a double to the container

    \param name : Name of data to store data with
    \param data : one double to be added to the container

    \note If record with name name already exists, it will be overwritten
    */
    void Add(const std::string& name, const double& data);

    /*!
    \brief Add vector of double to the container

    \param name : Name of data to store data with
    \param data : ptr to beginning of data
    \param num  : Number of objects in data

    \note If record with name name already exists, it will be overwritten
    */
    void Add(const std::string& name, const double* data, const std::size_t num);

    /*!
    \brief Add vector of double to the container

    \param name : Name of data to store data with
    \param data : vector of doubles to be added

    \note If record with name name already exists, it will be overwritten
    */
    void Add(const std::string& name, const std::vector<double>& data);

    /*!
    \brief Add a map of vectors of doubles to the container

    \param name : Name of data to store data with
    \param data : map of vector of doubles to be added

    \note If record with name name already exists, it will be overwritten
    */
    void Add(const std::string& name, const std::map<int, std::vector<double>>& data);

    /*!
    \brief Add a string to the container

    \param name : Name of data to store data with
    \param data : string to store

    \note If record with name name already exists, it will be overwritten
    */
    void Add(const std::string& name, const std::string& data);

    /*!
    \brief Add a CORE::LINALG::SerialDenseMatrix to the container

    \param name : Name of data to store data with
    \param matrix : matrix to store

    \note If record with name name already exists, it will be overwritten
    */
    void Add(const std::string& name, const CORE::LINALG::SerialDenseMatrix& matrix);

    //@}

    //! @name Access methods

    /*!
    \brief Get data of type int, double, string or CORE::LINALG::SerialDenseMatrix, read only

    Receive <int>, <double>, <string>, <CORE::LINALG::SerialDenseMatrix> data of a given name.

    \param name (in): Name of data to receive

    Usage is<br>
      \code
      const std::vector<int>* ifool = Container::Get<std::vector<int> >("ifool_name");
      \endcode
      or <br>
      \code
      const std::vector<double>* dfool = Container::Get<std::vector<double> >("dfool_name");
      \endcode
      or <br>
      \code
      const string* sfool = Container::Get<string>("sfool_name");
      \endcode
      or <br>
      \code
      const CORE::LINALG::SerialDenseMatrix* mfool =
    Container::Get<CORE::LINALG::SerialDenseMatrix>("mfool_name"); \endcode

    \note This template will let you experience one of the most kryptic error messages
          (at link time) you have ever seen if you try to use it with other data types than
          mentioned above. The same holds for the case where you did not get the
          syntax absolutely correct.

    \note (Advanced users:) This is a template with no general definition but with four
                            specializations only. It will compile for any data type
                            but will discover at link time that it is unable to generate
                            a definition other than one of the four specializations.

    \return const Ptr to data if data with that name exists, nullptr otherwise
    */
    template <typename T>
    const T* Get(const std::string& name) const
    {
      const auto access = [](const auto& map, const std::string& name) -> const T*
      {
        const auto it = map.find(name);
        if (it != map.end())
          return std::addressof(it->second);
        else
          return nullptr;
      };

      if constexpr (std::is_same_v<T, std::vector<int>>)
        return access(intdata_, name);
      else if constexpr (std::is_same_v<T, std::vector<double>>)
        return access(doubledata_, name);
      else if constexpr (std::is_same_v<T, std::map<int, std::vector<double>>>)
        return access(mapdata_, name);
      else if constexpr (std::is_same_v<T, std::string>)
        return access(stringdata_, name);
      else if constexpr (std::is_same_v<T, CORE::LINALG::SerialDenseMatrix>)
        return access(matdata_, name);
      else
        static_assert(!sizeof(T), "The Get function does not support this type.");
    };

    /*!
    \brief Get a single integer back without all that hassle

    Receive <int>

    \param name (in): Name of data to receive

    Usage is<br>
    \code
    int mynumber = Container::Getint("name_of_my_number")
    \endcode
    */
    [[nodiscard]] int GetInt(const std::string& name) const;

    /*!
    \brief Get a single double back without all that hassle

    Receive <int>

    \param name (in): Name of data to receive

    Usage is<br>
    \code
    double mynumber = Container::GetDouble("name_of_my_number")
    \endcode
    */
    [[nodiscard]] double GetDouble(const std::string& name) const;

    /*!
    \brief Get data of type int, double, string or CORE::LINALG::SerialDenseMatrix

    Receive <int>, <double>, <string>, <CORE::LINALG::SerialDenseMatrix> data of a given name.

    \param name (in): Name of data to receive

    Usage is<br>
      \code
      std::vector<int>* ifool = Container::Get<std::vector<int> >("ifool_name");
      \endcode
      or <br>
      \code
      std::vector<double>* dfool = Container::Get<std::vector<double> >("dfool_name");
      \endcode
      or <br>
      \code
      string* sfool = Container::Get<string>("sfool_name");
      \endcode
      or <br>
      \code
      CORE::LINALG::SerialDenseMatrix* mfool =
    Container::Get<CORE::LINALG::SerialDenseMatrix>("mfool_name"); \endcode

    \note This template will let you experience one of the most kryptic error messages
          (at link time) you have ever seen if you try to use it with other data types than
          mentioned above. The same holds for the case where you did not get the
          syntax absolutely correct.

    \note (Advanced users:) This is a template with no general definition but with four
                            specializations only. It will compile for any data type
                            but will discover at link time that it is unable to generate
                            a definition other than one of the four specializations.

    \return non-constant Ptr to data if data with that name exists, nullptr otherwise
    */
    template <typename T>
    T* Get(const std::string& name)
    {
      return const_cast<T*>(const_cast<const Container*>(this)->Get<T>(name));
    };
    //@}

    //! @name Erase methods

    /*!
    \brief Delete a record from the container

    \param name : Name of data to be deleted from container

    */
    void Delete(const std::string& name);

    /*!
    \brief Delete everything in this container
    */
    void Clear();

    //@}


   private:
    //! a map to store integer data in
    std::map<std::string, std::vector<int>> intdata_;

    //! a map to store double data in
    std::map<std::string, std::vector<double>> doubledata_;

    //! a map to store maps of stl vector data
    std::map<std::string, std::map<int, std::vector<double>>> mapdata_;

    //! a map to store string data in
    std::map<std::string, std::string> stringdata_;

    //! a map to store matrices in
    std::map<std::string, CORE::LINALG::SerialDenseMatrix> matdata_;
  };  // class Container
}  // namespace DRT

// << operator
std::ostream& operator<<(std::ostream& os, const DRT::Container& cont);



BACI_NAMESPACE_CLOSE

#endif  // LIB_CONTAINER_H
