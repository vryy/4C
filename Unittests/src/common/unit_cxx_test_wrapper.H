/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper that is the interface between BACI unit tests and CXX test

\level 0

*/

#ifndef UNIT_CXX_TEST_WRAPPER_H
#define UNIT_CXX_TEST_WRAPPER_H

#include <cxxtest/TestSuite.h>

//! This class serves as an interface between the Cxx test suite and BACI.
//! All classes in Unittests/src that provide functions beginning with 'Test' should derive from
//! this class and override the Setup() and TearDown() call
//! This class aims to use BACI's case style.
class BACICxxTestWrapper : public CxxTest::TestSuite
{
 public:
  void setUp() final { Setup(); }
  void tearDown() final { TearDown(); }

  virtual void Setup(){};
  virtual void TearDown(){};
};

#endif  // UNIT_CXX_TEST_WRAPPER_H