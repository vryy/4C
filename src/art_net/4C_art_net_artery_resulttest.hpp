/*----------------------------------------------------------------------*/
/*! \file

\brief testing of artery calculation results


\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ART_NET_ARTERY_RESULTTEST_HPP
#define FOUR_C_ART_NET_ARTERY_RESULTTEST_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_utils_result_test.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Arteries
{
  // forward declaration
  class ArtNetExplicitTimeInt;
  class ArtNetImplStationary;

  /*!
    \brief artnet specific result test class

    \author Mahmoud Ismail
    \date 11/11
  */
  class ArteryResultTest : public Core::Utils::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    ArteryResultTest(ArtNetExplicitTimeInt& art_net);

    /*!
    \brief constructor
    */
    ArteryResultTest(ArtNetImplStationary& art_net);


    /// our version of nodal value tests
    /*!
      Possible position flags is only "phi"
     */
    void test_node(
        const Core::IO::InputParameterContainer& container, int& nerr, int& test_count) override;

    /// our version of element value tests
    void test_element(
        const Core::IO::InputParameterContainer& container, int& nerr, int& test_count) override;

   private:
    /// Teuchos::RCP to scalar transport discretization
    Teuchos::RCP<Core::FE::Discretization> dis_;
    /// Teuchos::RCP to solution vector
    Teuchos::RCP<const Core::LinAlg::Vector<double>> mysol_;
    /// Teuchos::RCP to element volumetric flow
    Teuchos::RCP<const Core::LinAlg::Vector<double>> myelevolflow_;
    /// Teuchos::RCP to element radius
    Teuchos::RCP<const Core::LinAlg::Vector<double>> myeleradius_;
  };

}  // namespace Arteries

FOUR_C_NAMESPACE_CLOSE

#endif
