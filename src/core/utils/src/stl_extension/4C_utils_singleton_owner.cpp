/*----------------------------------------------------------------------------*/
/*! \file

\brief  SingletonOwner that manages singleton instances

\level 0
    */
/*----------------------------------------------------------------------------*/

#include "4C_utils_singleton_owner.hpp"

void FourC::Core::UTILS::SingletonOwnerRegistry::finalize()
{
  for (const auto& [_, deleter] : instance().deleters_)
  {
    deleter();
  }
}

void FourC::Core::UTILS::SingletonOwnerRegistry::register_deleter(
    void* owner, std::function<void()> deleter)
{
  instance().deleters_.emplace(owner, std::move(deleter));
}

void FourC::Core::UTILS::SingletonOwnerRegistry::unregister(void* owner)
{
  instance().deleters_.erase(owner);
}

void FourC::Core::UTILS::SingletonOwnerRegistry::initialize() { instance(); }

FourC::Core::UTILS::SingletonOwnerRegistry& FourC::Core::UTILS::SingletonOwnerRegistry::instance()
{
  static SingletonOwnerRegistry instance;
  return instance;
}
