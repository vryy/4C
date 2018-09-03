/*----------------------------------------------------------------------*/
/**
\file any_data_container.H

\brief Any data container based on vectors of Teuchos::any data.

\maintainer Michael Hiermeier

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef SRC_HEADERS_ANY_DATA_CONTAINER_H_
#define SRC_HEADERS_ANY_DATA_CONTAINER_H_

#include <Teuchos_any.hpp>
#include "../drt_lib/drt_dserror.H"
#include <unordered_map>


namespace GEN
{
  /** \brief Data container of any content
   *
   *  The AnyDataContainer is meant as a container class of any content
   *  which can be used to exchange data of any type and any quantity over
   *  different function calls and/or class hierarchies without the need to
   *  adapt or extend your code. The container is NOT meant as a collective data
   *  storage but more as a vehicle to transport your data to the place you
   *  need it.
   *
   *  \author hiermeier \date 12/17 */
  class AnyDataContainer
  {
    //! alias templates
    //! @{

    //! alias for an unordered_map
    template <typename... Ts>
    using UMap = std::unordered_map<Ts...>;

    //! alias for a vector
    template <typename... Ts>
    using Vec = std::vector<Ts...>;

    //! @}

    struct AnyData
    {
      AnyData() = default;

      AnyData(const Teuchos::any& data) : data_(data){/* empty */};

      ~AnyData() = default;

      const Teuchos::any& get() const
      {
        if (data_.empty()) dserror("The data is empty!");

        return data_;
      }

      template <typename T>
      void set(const T& data)
      {
        data_ = data;
      }

     private:
      Teuchos::any data_;
    };

   public:
    enum class DataType
    {
      vague,
      any,
      vector,
      unordered_map
    };

   public:
    AnyDataContainer() = default;
    ~AnyDataContainer() = default;

    /// @name any data
    /// @{

    template <typename T>
    void Set(const T* data, const unsigned id = 0)
    {
      SetData<T, DataType::any>(data, id);
    }

    template <typename T>
    T* Get(const unsigned id = 0)
    {
      return const_cast<T*>(GetData<T, DataType::any>(id));
    }

    template <typename T>
    const T* Get(const unsigned id = 0) const
    {
      return GetData<T, DataType::any>(id);
    }

    /// @}

    /// @name std::vector
    /// @{

    template <typename... Ts>
    void SetVector(const Vec<Ts...>* unordered_map, const unsigned id = 0)
    {
      SetData<Vec<Ts...>, DataType::vector>(unordered_map, id);
    }

    template <typename... Ts>
    Vec<Ts...>* GetVector(const unsigned id = 0)
    {
      return const_cast<Vec<Ts...>*>(GetData<Vec<Ts...>, DataType::vector>(id));
    }

    template <typename... Ts>
    const Vec<Ts...>* GetVector(const unsigned id = 0) const
    {
      return GetData<Vec<Ts...>, DataType::vector>(id);
    }

    /// @}

    /// @name std::unordered_map
    /// @{

    template <typename... Ts>
    void SetUnorderedMap(const UMap<Ts...>* unordered_map, const unsigned id = 0)
    {
      SetData<UMap<Ts...>, DataType::unordered_map>(unordered_map, id);
    }

    template <typename... Ts>
    UMap<Ts...>* GetUnorderedMap(const unsigned id = 0)
    {
      return const_cast<UMap<Ts...>*>(GetData<UMap<Ts...>, DataType::unordered_map>(id));
    }

    template <typename... Ts>
    const UMap<Ts...>* GetUnorderedMap(const unsigned id = 0) const
    {
      return GetData<UMap<Ts...>, DataType::unordered_map>(id);
    }

    /// @}

    /// @name general methods
    /// @{

    void ClearEntry(const DataType type, const int id)
    {
      switch (type)
      {
        case DataType::vector:
        {
          Clear(vector_data_, id);

          break;
        }
        case DataType::unordered_map:
        {
          Clear(unordered_map_data_, id);

          break;
        }
        default:
        {
          dserror("Unsupported DataType!");
          exit(EXIT_FAILURE);
        }
      }
    }

    void ClearAll(const DataType type) { ClearEntry(type, -1); }

    /// @}

   private:
    void Clear(std::vector<AnyData>& any_data_vec, const int id)
    {
      // clear all entries
      if (id < 0)
      {
        any_data_vec.clear();
        return;
      }

      // clear only one entry
      if (id >= static_cast<int>(any_data_vec.size()))
        dserror("Requested ID #%d exceeds the AnyData vector size (=%d).", id, any_data_vec.size());
      else if (id == (static_cast<int>(any_data_vec.size()) - 1))
        any_data_vec.pop_back();
      else
        any_data_vec[id] = Teuchos::any();
    }

    template <typename T, DataType type>
    void SetData(const T* data, const unsigned id)
    {
      Teuchos::any any_data(data);
      SetAnyData<type>(any_data, id);
    }

    template <DataType type>
    void SetAnyData(const Teuchos::any& any_data, const unsigned id)
    {
      switch (type)
      {
        case DataType::vector:
        {
          AddToAnyDataVec(any_data, id, vector_data_);

          break;
        }
        case DataType::unordered_map:
        {
          AddToAnyDataVec(any_data, id, unordered_map_data_);

          break;
        }
        case DataType::any:
        {
          AddToAnyDataVec(any_data, id, any_data_);

          break;
        }
        default:
          dserror("Unsupported DataType!");
      }
    }

    template <typename T, DataType type>
    const T* GetData(const unsigned id) const
    {
      const Teuchos::any& any_data = GetAnyData<type>(id);

      //      any_data.print( std::cout );
      //      std::cout << any_data.typeName() << std::endl;
      return Teuchos::any_cast<const T*>(any_data);
    }

    template <DataType type>
    const Teuchos::any& GetAnyData(const unsigned id) const
    {
      switch (type)
      {
        case DataType::vector:
        {
          return GetFromAnyDataVec(id, vector_data_);
        }
        case DataType::unordered_map:
        {
          return GetFromAnyDataVec(id, unordered_map_data_);
        }
        case DataType::any:
        {
          return GetFromAnyDataVec(id, any_data_);
        }
        default:
        {
          dserror("Unsupported DataType!");
          exit(EXIT_FAILURE);
        }
      }
    }

    void AddToAnyDataVec(
        const Teuchos::any& any_data, const unsigned id, std::vector<AnyData>& any_data_vec) const
    {
      if (any_data_vec.size() <= id) any_data_vec.resize(id + 1);

      AnyData& data_id = any_data_vec[id];
      data_id.set(any_data);
    }

    inline const Teuchos::any& GetFromAnyDataVec(
        const unsigned id, const std::vector<AnyData>& any_data_vec) const
    {
      if (id >= any_data_vec.size())
        dserror("Requested ID #%d exceeds the AnyData vector size (=%d).", id, any_data_vec.size());

      return any_data_vec[id].get();
    }

   private:
    std::vector<AnyData> vector_data_;

    std::vector<AnyData> unordered_map_data_;

    std::vector<AnyData> any_data_;
  };
}  // namespace GEN


#endif /* SRC_HEADERS_ANY_DATA_CONTAINER_H_ */
