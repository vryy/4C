/*----------------------------------------------------------------------*/
/*! \file
\brief special assertions for BACI code
\level 1
*----------------------------------------------------------------------*/
#ifndef UNITTESTS_COMMON_ASSERTIONS_H
#define UNITTESTS_COMMON_ASSERTIONS_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>


/*!
 * Extension of EXPECT_THROW which also checks for a substring in the what() expression.
 */
#define BACI_EXPECT_THROW_WITH_MESSAGE(statement, expectedException, messageSubString) \
  std::function<void(void)> find_the_statement_below([&]() {                           \
    try                                                                                \
    {                                                                                  \
      statement;                                                                       \
    }                                                                                  \
    catch (const expectedException& caughtException)                                   \
    {                                                                                  \
      using ::testing::HasSubstr;                                                      \
      EXPECT_THAT(caughtException.what(), HasSubstr(messageSubString))                 \
          << "Caught the expected exception type but message has wrong substring.";    \
      throw;                                                                           \
    }                                                                                  \
  });                                                                                  \
  EXPECT_THROW(find_the_statement_below(), expectedException) << "statement: " << #statement;

#endif  // UNITTESTS_COMMON_ASSERTIONS_H
