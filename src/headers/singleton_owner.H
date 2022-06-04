/*----------------------------------------------------------------------------*/
/*! \file

\brief  SingletonOwner that manages singleton instances

\level 0
*/
/*----------------------------------------------------------------------------*/

#ifndef SINGLETON_OWNER_H
#define SINGLETON_OWNER_H

#include <functional>
#include <memory>

namespace UTILS
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
   */
  template <typename T>
  class SingletonOwner
  {
   public:
    /*!
     * @brief Constructor
     *
     * @param [in] creator  Function that can create the singleton instance and return a unique_ptr
     * to it
     */
    SingletonOwner(std::function<std::unique_ptr<T>()> creator) : creator_(creator) {}

    /*!
     * @brief Return pointer to singleton instance
     *
     * If the @p action is SingletonAction::create, create the singleton instance (if not already
     * created) and return the pointer to it. If the @p action is SingletonAction::destruct,
     * destroy the instance and return a `nullptr`.
     *
     */
    T* Instance(SingletonAction action)
    {
      if (action == SingletonAction::create and !instance_)
      {
        instance_ = creator_();
      }
      else if (action == SingletonAction::destruct)
      {
        instance_.reset();
      }
      return instance_.get();
    }

   private:
    //! singleton instance
    std::unique_ptr<T> instance_;

    //! Function that creates a singleton object
    std::function<std::unique_ptr<T>()> creator_;
  };
}  // namespace UTILS


#endif  // SINGLETON_OWNER_H
