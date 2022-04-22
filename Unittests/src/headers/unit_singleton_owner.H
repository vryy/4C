/*! \file

\brief unit tests for SingletonOwner

\level 3
*/

#ifndef BACI_UNIT_SINGLETON_OWNER_H
#define BACI_UNIT_SINGLETON_OWNER_H

#include "Unittests/src/common/unit_cxx_test_wrapper.H"
#include "src/headers/singleton_owner.H"
#include <memory>


namespace UTILS
{
  class SingletonOwner_TestSuite;

  class DummySingleton
  {
   private:
    DummySingleton() { ; }

    friend class SingletonOwner_TestSuite;
  };
}  // namespace UTILS

class UTILS::SingletonOwner_TestSuite : public BACICxxTestWrapper
{
 public:
  void Setup() override { ; }

  void TearDown() override { ; }

  void TestCreation()
  {
    SingletonOwner<DummySingleton> singleton_owner(
        []() { return std::unique_ptr<DummySingleton>(new DummySingleton()); });

    // Assert that the returned object is of DummySingleton type
    TS_ASSERT(
        dynamic_cast<DummySingleton*>(singleton_owner.Instance(::UTILS::SingletonAction::create)));
  }

  void TestDestruction()
  {
    SingletonOwner<DummySingleton> singleton_owner(
        []() { return std::unique_ptr<DummySingleton>(new DummySingleton()); });

    // Assert that a nullptr is returned
    TS_ASSERT(not singleton_owner.Instance(::UTILS::SingletonAction::destruct));
  }

  void TestReturnExistingInstance()
  {
    SingletonOwner<DummySingleton> singleton_owner(
        []() { return std::unique_ptr<DummySingleton>(new DummySingleton()); });

    DummySingleton* ptr_1 = singleton_owner.Instance(::UTILS::SingletonAction::create);
    DummySingleton* ptr_2 = singleton_owner.Instance(::UTILS::SingletonAction::create);

    // Assert that both pointers point to the same object
    TS_ASSERT(ptr_1 == ptr_2);
  }
};

#endif  // BACI_UNIT_SINGLETON_OWNER_H
