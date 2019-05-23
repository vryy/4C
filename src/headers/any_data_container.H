/*----------------------------------------------------------------------*/
/*!

\brief Any data container based on vectors of Teuchos::any data.

\maintainer Martin Kronbichler

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef SRC_HEADERS_ANY_DATA_CONTAINER_H_
#define SRC_HEADERS_ANY_DATA_CONTAINER_H_

#include <Teuchos_any.hpp>
#include "../drt_lib/drt_dserror.H"
#include <unordered_map>

// forward declaration
namespace CONTACT
{
  namespace AUG
  {
    template <typename enum_class>
    class TimeMonitor;
  }  // namespace AUG
}  // namespace CONTACT

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

    //! alias for the time monitor
    template <typename... Ts>
    using TimeMonitor = CONTACT::AUG::TimeMonitor<Ts...>;

    //! @}

    /// internal data container
    struct AnyData
    {
      /// default constructor
      AnyData() = default;

      /// copy constructor
      AnyData(const Teuchos::any& data) : data_(data){/* empty */};

      /// destructor
      ~AnyData() = default;

      /// accessor
      const Teuchos::any& get() const
      {
        if (data_.empty()) dserror("The data is empty!");

        return data_;
      }

      /// set any data and perform a sanity check
      template <typename T>
      void set(const T& data)
      {
        if (data_.empty())
          data_ = data;
        else
          dserror(
              "There are already data:\n%s\nAre you sure, that you want "
              "to overwrite them? If yes: Clear the content first. Best practice "
              "is to clear the content right after you finished the "
              "respective task.",
              data_.typeName().c_str());
      }

     private:
      /// actual stored data (of any type)
      Teuchos::any data_;
    };

   public:
    /// supported (more specific) data types
    enum class DataType
    {
      vague,
      any,
      vector,
      unordered_map,
      time_monitor
    };

   public:
    /// constructor
    AnyDataContainer() = default;

    /// destructor
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

    /// @name Time monitoring
    /// @{

    template <typename... Ts>
    void SetTimer(const CONTACT::AUG::TimeMonitor<Ts...>* timer, const unsigned id = 0)
    {
      SetData<TimeMonitor<Ts...>, DataType::time_monitor>(timer, id);
    }

    template <typename... Ts>
    TimeMonitor<Ts...>* GetTimer(const unsigned id = 0)
    {
      return const_cast<TimeMonitor<Ts...>*>(
          GetData<TimeMonitor<Ts...>, DataType::time_monitor>(id));
    }

    template <typename... Ts>
    const TimeMonitor<Ts...>* GetTimer(const unsigned id = 0) const
    {
      return GetData<TimeMonitor<Ts...>, DataType::time_monitor>(id);
    }

    /// @}

    /// @name general methods
    /// @{

    /// clear an entry in the respective container
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
        case DataType::time_monitor:
        {
          Clear(time_monitor_data_, id);

          break;
        }
        case DataType::any:
        {
          Clear(any_data_, id);

          break;
        }
        default:
        {
          dserror("Unsupported DataType!");
          exit(EXIT_FAILURE);
        }
      }
    }

    // clear all entries in the respective container
    void ClearAll(const DataType type) { ClearEntry(type, -1); }

    /// @}

   private:
    /// helper function to clear content
    void Clear(std::vector<AnyData>& any_data_vec, const int id)
    {
      // clear all entries
      if (id < 0)
      {
        for (auto& any_data : any_data_vec) any_data = Teuchos::any();

        return;
      }

      // direct return if the id exceeds the vector size
      if (id >= static_cast<int>(any_data_vec.size())) return;

      // clear only one entry
      any_data_vec[id] = Teuchos::any();
    }

    /// pack and set the data pointer
    template <typename T, DataType type>
    void SetData(const T* data, const unsigned id)
    {
      Teuchos::any any_data(data);
      SetAnyData<type>(any_data, id);
    }

    /// set the data in the respective container
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
        case DataType::time_monitor:
        {
          AddToAnyDataVec(any_data, id, time_monitor_data_);

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

    /// access the data and cast the any pointer
    template <typename T, DataType type>
    const T* GetData(const unsigned id) const
    {
      const Teuchos::any& any_data = GetAnyData<type>(id);

      //      any_data.print( std::cout );
      //      std::cout << any_data.typeName() << std::endl;
      return Teuchos::any_cast<const T*>(any_data);
    }

    /// access the data
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
        case DataType::time_monitor:
        {
          return GetFromAnyDataVec(id, time_monitor_data_);
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

    /// add to any data vector
    void AddToAnyDataVec(
        const Teuchos::any& any_data, const unsigned id, std::vector<AnyData>& any_data_vec) const
    {
      if (any_data_vec.size() <= id) any_data_vec.resize(id + 1);

      AnyData& data_id = any_data_vec[id];
      data_id.set(any_data);
    }

    /// access content of any data vector
    inline const Teuchos::any& GetFromAnyDataVec(
        const unsigned id, const std::vector<AnyData>& any_data_vec) const
    {
      if (id >= any_data_vec.size())
        dserror("Requested ID #%d exceeds the AnyData vector size (=%d).", id, any_data_vec.size());

      return any_data_vec[id].get();
    }

   private:
    /// specific container for vector data
    std::vector<AnyData> vector_data_;

    /// specific container for unordered map data
    std::vector<AnyData> unordered_map_data_;

    /// specific container for time monitor data
    std::vector<AnyData> time_monitor_data_;

    /// container for any data
    std::vector<AnyData> any_data_;
  };
}  // namespace GEN


#endif /* SRC_HEADERS_ANY_DATA_CONTAINER_H_ */
