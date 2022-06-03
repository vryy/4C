/*----------------------------------------------------------------------*/
/*! \file

\brief unit tests for SingletonOwner

\level 1

*-----------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include "src/headers/singleton_owner.H"

namespace
{
  class DummySingleton
  {
   public:
    DummySingleton() = default;
  };

  TEST(TestSingletonOwner, CreatesSingleton)
  {
    ::UTILS::SingletonOwner<DummySingleton> singleton_owner(
        []() { return std::unique_ptr<DummySingleton>(new DummySingleton()); });

    // Expect that the returned object is of DummySingleton type
    EXPECT_TRUE(
        dynamic_cast<DummySingleton*>(singleton_owner.Instance(::UTILS::SingletonAction::create)));
  }

  TEST(TestSingletonOwner, DestructsSingleton)
  {
    ::UTILS::SingletonOwner<DummySingleton> singleton_owner(
        []() { return std::unique_ptr<DummySingleton>(new DummySingleton()); });

    // Create a singleton to destruct it in the following
    singleton_owner.Instance(::UTILS::SingletonAction::create);

    // Expect that a nullptr is returned at destruction
    EXPECT_EQ(singleton_owner.Instance(::UTILS::SingletonAction::destruct), nullptr);
  }

  TEST(TestSingletonOwner, ReturnsExistingInstance)
  {
    ::UTILS::SingletonOwner<DummySingleton> singleton_owner(
        []() { return std::unique_ptr<DummySingleton>(new DummySingleton()); });

    DummySingleton* ptr_1 = singleton_owner.Instance(::UTILS::SingletonAction::create);
    DummySingleton* ptr_2 = singleton_owner.Instance(::UTILS::SingletonAction::create);

    // Expect that both pointers point to the same object
    EXPECT_EQ(ptr_1, ptr_2);
  }
}  // namespace