/*----------------------------------------------------------------------------*/
/*! \file

\brief Result tests for pure ALE problems

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef FOUR_C_ALE_RESULTTEST_HPP
#define FOUR_C_ALE_RESULTTEST_HPP

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_config.hpp"

#include "4C_adapter_ale.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_result_test.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace ALE
{
  class Ale;
}

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace ALE
{
  /// Result test subclass for linear ale algorithm
  class AleResultTest : public Core::UTILS::ResultTest
  {
   public:
    AleResultTest(ALE::Ale& ale);

    /// our version of nodal value tests
    /*!
      Possible position flags are "velx", "vely", "velz" and
      "pressure". With the obvious meaning.
     */
    void test_node(
        const Core::IO::InputParameterContainer& container, int& nerr, int& test_count) override;

   private:
    /// pointer to ALE discretization
    Teuchos::RCP<const Core::FE::Discretization> aledis_;

    /// pointer to ALE displacement result vector
    Teuchos::RCP<const Core::LinAlg::Vector> dispnp_;
  };

}  // namespace ALE

FOUR_C_NAMESPACE_CLOSE

#endif
