/*----------------------------------------------------------------------*/
/*! \file

\brief unit tests for SingletonOwner

\level 1

*-----------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "baci_utils_singleton_owner.H"

namespace
{
  class DummySingleton
  {
   private:
    //! Private constructor to mimic typical use case of singletons.
    DummySingleton() = default;

    FRIEND_TEST(TestSingletonOwner, CreatesSingleton);
    FRIEND_TEST(TestSingletonOwner, DestructsSingleton);
    FRIEND_TEST(TestSingletonOwner, ReturnsExistingInstance);
    FRIEND_TEST(TestSingletonMap, DifferentKeys);
  };

  TEST(TestSingletonOwner, CreatesSingleton)
  {
    auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
        []() { return std::unique_ptr<DummySingleton>(new DummySingleton()); });

    // Expect that the returned object is of DummySingleton type
    EXPECT_TRUE(dynamic_cast<DummySingleton*>(
        singleton_owner.Instance(CORE::UTILS::SingletonAction::create)));
  }

  TEST(TestSingletonOwner, DestructsSingleton)
  {
    CORE::UTILS::SingletonOwner<DummySingleton> singleton_owner(
        []() { return std::unique_ptr<DummySingleton>(new DummySingleton()); });

    // Create a singleton to destruct it in the following
    singleton_owner.Instance(CORE::UTILS::SingletonAction::create);

    // Expect that a nullptr is returned at destruction
    EXPECT_EQ(singleton_owner.Instance(CORE::UTILS::SingletonAction::destruct), nullptr);
  }

  TEST(TestSingletonOwner, ReturnsExistingInstance)
  {
    struct Creator
    {
      MOCK_METHOD((std::unique_ptr<DummySingleton>), create, (), (const));
    };
    Creator creator;

    // Return a new DummySingleton exactly once, otherwise the test fails
    EXPECT_CALL(creator, create)
        .WillOnce([]() { return std::unique_ptr<DummySingleton>(new DummySingleton()); });

    CORE::UTILS::SingletonOwner<DummySingleton> singleton_owner{
        [&creator]() { return creator.create(); }};

    DummySingleton* ptr_1 = singleton_owner.Instance(CORE::UTILS::SingletonAction::create);
    DummySingleton* ptr_2 = singleton_owner.Instance(CORE::UTILS::SingletonAction::create);

    // Expect that both pointers point to the same object
    EXPECT_EQ(ptr_1, ptr_2);
  }

  TEST(TestSingletonMap, DifferentKeys)
  {
    struct Creator
    {
      MOCK_METHOD((std::unique_ptr<DummySingleton>), create, (), (const));
    };
    Creator creator;

    // Return a new DummySingleton exactly twice (for the two keys tested below)
    EXPECT_CALL(creator, create)
        .Times(2)
        .WillRepeatedly([]() { return std::unique_ptr<DummySingleton>(new DummySingleton()); });

    auto singleton_map =
        CORE::UTILS::MakeSingletonMap<std::string>([&creator]() { return creator.create(); });


    auto* a = singleton_map["a"].Instance(CORE::UTILS::SingletonAction::create);
    auto* b = singleton_map["b"].Instance(CORE::UTILS::SingletonAction::create);

    EXPECT_NE(a, b);
    EXPECT_EQ(singleton_map["a"].Instance(CORE::UTILS::SingletonAction::create), a);
    EXPECT_EQ(singleton_map["b"].Instance(CORE::UTILS::SingletonAction::create), b);
  }

  TEST(TestSingletonMap, ForwardConstructorArgs)
  {
    struct DummyWithArgs
    {
      DummyWithArgs(int a) : a(a) {}

      int a;
    };
    auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
        [](int input) { return std::make_unique<DummyWithArgs>(input); });

    EXPECT_EQ(singleton_map["a"].Instance(CORE::UTILS::SingletonAction::create, 2)->a, 2);
  }
}  // namespace