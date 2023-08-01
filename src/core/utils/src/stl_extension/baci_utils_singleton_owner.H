/*----------------------------------------------------------------------------*/
/*! \file

\brief  SingletonOwner that manages singleton instances

\level 0
*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_UTILS_SINGLETON_OWNER_H
#define BACI_UTILS_SINGLETON_OWNER_H

#include <functional>
#include <map>
#include <memory>
#include <type_traits>

namespace CORE::UTILS
{
  //! Action types for the singleton owner
  enum class SingletonAction
  {
    create,   //!< Create an instance (if not already created)
    destruct  //!< Destruct an existing instance
  };

  /*!
   * @brief Owner class for singletons
   *
   * Manages (creates, destructs, or returns) the singleton instance of any type.
   *
   * @tparam T  Type of the singleton
   * @tparam CreationArgs A template parameter pack which contains the types of all arguments that
   *   are necessary to create the owned singleton object.
   */
  template <typename T, typename... CreationArgs>
  class SingletonOwner
  {
   public:
    /*!
     * @brief Constructor.
     *
     * Construction is preferably done via MakeSingletonOwner() which comes with powerful template
     * argument deduction.
     *
     * @param [in] creator  Function (object) that can create the singleton instance and return a
     * unique_ptr to it. The function's signature must be `std::unique_ptr<T>(CreationArgs...)`.
     */
    template <typename Fn, typename E = std::enable_if_t<!std::is_same_v<Fn, SingletonOwner>, Fn>>
    SingletonOwner(Fn&& creator);

    //! Deleted copy constructor. Only one object can own the singleton.
    SingletonOwner(const SingletonOwner& other) = delete;

    //! Allow move construction.
    SingletonOwner(SingletonOwner&& other) noexcept = default;

    //! Deleted copy assignment. Only one object can own the singleton.
    SingletonOwner& operator=(const SingletonOwner& other) = delete;

    //! Allow move assignment.
    SingletonOwner& operator=(SingletonOwner&& other) noexcept = default;

    /*!
     * @brief Return pointer to singleton instance
     *
     * If the @p action is SingletonAction::create, create the singleton instance (if not already
     * created) and return the pointer to it. If the @p action is SingletonAction::destruct,
     * destroy the instance and return a `nullptr`.
     *
     * In case any parameters are needed to call the constructor of the singleton object, they have
     * to be passed as well. Unfortunately, this is necessary regardless of the @p action flag.
     *
     */
    T* Instance(SingletonAction action, CreationArgs... args);

   private:
    //! singleton instance
    std::unique_ptr<T> instance_;

    //! Function that creates a singleton object
    std::function<std::unique_ptr<T>(CreationArgs...)> creator_;
  };

  /**
   * @brief Store multiple SingletonOwner objects by a given Key.
   *
   * This class is a useful extension to SingletonOwner in places where multiple differently
   * parametrized objects of the same Singleton class should be used. This classes Instance() method
   * takes an additional Key argument to distinguish between the singletons and apart from this
   * difference behaves the same as SingletonOwner::Instance().
   *
   * @tparam Key The type used to access the different instances.
   * @tparam T The type of the singleton.
   * @tparam CreationArgs A template parameter pack which contains the types of all arguments that
   *   are necessary to create the owned singleton object.
   */
  template <typename Key, typename T, typename... CreationArgs>
  class SingletonMap
  {
   public:
    /*!
     * @brief Constructor
     *
     * Construction is preferably done via MakeSingletonMap() which comes with powerful template
     * argument deduction.
     *
     * @param [in] creator  Function (object) that can create the singleton instance and return a
     * unique_ptr to it. The function's signature must be `std::unique_ptr<T>(CreationArgs...)`.
     */
    template <typename Fn, typename E = std::enable_if_t<!std::is_same_v<Fn, SingletonMap>, Fn>>
    SingletonMap(Fn&& creator);

    /**
     * Return a SingletonOwner for the given @p key. If it does not exist, one is created the first
     * time this function is called with the @p key. This call is normally directly followed by a
     * call to SingletonOwner::Instance() e.g.
     *
     * @code
     *   singleton_map[disname].Instance(CORE::UTILS::SingletonAction::create, numdofpernode,
     * disname);
     * @endoce
     */
    SingletonOwner<T, CreationArgs...>& operator[](const Key& key);

   private:
    //! Function that creates a singleton object
    std::function<std::unique_ptr<T>(CreationArgs...)> creator_;

    //! All SingletonOwner objects that are stored internally.
    std::map<Key, SingletonOwner<T, CreationArgs...>> map_;
  };


  /**
   * Convenience function to create a SingletonOwner object. This function mainly exists to deduce
   * the template arguments of more complicated @p creator functions and pass them on to
   * SingletonOwner. This allows to write code like:
   *
   * @code
   *   // static initialization upon first call
   *   static auto singleton_owner = MakeSingletonOwner(
   *     [](double a, int b) { return std::make_unique<MyClass>(a,b); });
   *
   *   // access an instance constructed from given parameters (if it does not exist already)
   *   singleton_owner.Instance(SingletonAction::create, 1.0, 2);
   * @endcode
   *
   * @param creator Function (object) that can create the singleton instance and return a
   *   unique_ptr to it. The function's signature must be `std::unique_ptr<T>(CreationArgs...)`
   * @return A SingletonOwner object that should probably be stored as a static variable.
   */
  template <typename Fn>
  auto MakeSingletonOwner(Fn&& creator);

  /**
   * Similar to MakeSingletonOwner(), this function helps to deduce the template arguments for a
   * SingletonMap. Note that you need to specify the first template argument, which defines the type
   * of the key that is used to access the map entries.
   *
   * @code
   *
   *   // static initialization upon first call
   *  static auto singleton_map = ::UTILS::make_singleton_map<std::string>(
   *    [](const int numdofpernode, const int numscal, const std::string& disname)
   *    {
   *      return std::unique_ptr<ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>>(
   *         new ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>(
   *           numdofpernode, numscal, disname));
   *    });
   *
   *  // Use the map as
   *  singleton_map[disname].Instance(SingletonAction::create, numdofpernode, numscal, disname);
   *
   * @endcode
   *
   * @param creator Function (object) that can create the singleton instance and return a
   *   unique_ptr to it. The function's signature must be `std::unique_ptr<T>(CreationArgs...)`
   * @return A SingletonMap object that should probably be stored as a static variable.
   */
  template <typename KeyType, typename Fn>
  auto MakeSingletonMap(Fn&& creator);



  // --- template and inline functions --- //


  template <typename T, typename... CreationArgs>
  template <typename Fn, typename E>
  SingletonOwner<T, CreationArgs...>::SingletonOwner(Fn&& creator)
      : creator_(std::forward<Fn>(creator))
  {
  }


  template <typename T, typename... CreationArgs>
  T* SingletonOwner<T, CreationArgs...>::Instance(SingletonAction action, CreationArgs... args)
  {
    if (action == SingletonAction::create and !instance_)
    {
      instance_ = creator_(args...);
    }
    else if (action == SingletonAction::destruct)
    {
      instance_.reset();
    }
    return instance_.get();
  }


  template <typename Key, typename T, typename... CreationArgs>
  template <typename Fn, typename E>
  SingletonMap<Key, T, CreationArgs...>::SingletonMap(Fn&& creator)
      : creator_(std::forward<Fn>(creator))
  {
  }


  template <typename Key, typename T, typename... CreationArgs>
  SingletonOwner<T, CreationArgs...>& SingletonMap<Key, T, CreationArgs...>::operator[](
      const Key& key)
  {
    auto it = map_.find(key);
    if (it == map_.end())
    {
      it = map_.emplace(key, creator_).first;
    }

    return it->second;
  }


  /**
   * Implementation details.
   */
  namespace detail  // NOLINT
  {
    template <class F>
    struct DeduceSingletonOwner
    {
    };

    template <class Ret, class C, class... Args>
    struct DeduceSingletonOwner<Ret (C::*)(Args...)>
    {
      using ReturnType = typename Ret::element_type;
      using type = SingletonOwner<ReturnType, Args...>;
    };

    template <class Ret, class C, class... Args>
    struct DeduceSingletonOwner<Ret (C::*)(Args...) const>
    {
      using ReturnType = typename Ret::element_type;
      using type = SingletonOwner<ReturnType, Args...>;
    };

    template <typename KeyType, class F>
    struct DeduceSingletonMap
    {
    };

    template <typename KeyType, class Ret, class C, class... Args>
    struct DeduceSingletonMap<KeyType, Ret (C::*)(Args...)>
    {
      using ReturnType = typename Ret::element_type;
      using type = SingletonMap<KeyType, ReturnType, Args...>;
    };

    template <typename KeyType, class Ret, class C, class... Args>
    struct DeduceSingletonMap<KeyType, Ret (C::*)(Args...) const>
    {
      using ReturnType = typename Ret::element_type;
      using type = SingletonMap<KeyType, ReturnType, Args...>;
    };

  }  // namespace detail

  template <typename Fn>
  auto MakeSingletonOwner(Fn&& creator)
  {
    using T = typename detail::DeduceSingletonOwner<decltype(&Fn::operator())>::type;
    return T(std::forward<Fn>(creator));
  }


  template <typename KeyType, typename Fn>
  auto MakeSingletonMap(Fn&& creator)
  {
    using T = typename detail::DeduceSingletonMap<KeyType, decltype(&Fn::operator())>::type;
    return T(std::forward<Fn>(creator));
  }
}  // namespace CORE::UTILS


#endif
