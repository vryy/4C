#ifndef FOUR_C_THERMO_RESULTTEST_HPP
#define FOUR_C_THERMO_RESULTTEST_HPP


/*----------------------------------------------------------------------*
 | headers                                                   dano 08/09 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_thermo_timint.hpp"
#include "4C_utils_result_test.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | belongs to thermal dynamics namespace                     dano 08/09 |
 *----------------------------------------------------------------------*/
namespace Thermo
{
  //!
  //! \brief Thermo specific result test class
  //!
  //! \author cd
  class ResultTest : public Core::Utils::ResultTest
  {
   public:
    //! Constructor for time integrators of general kind
    //! \author bborn \date 06/08
    ResultTest(TimInt& tintegrator);

    //! \brief thermo version of nodal value tests
    //!
    //! Possible position flags are "temp",
    //!                             "rate",
    void test_node(
        const Core::IO::InputParameterContainer& container, int& nerr, int& test_count) override;

   private:
    //! our discretisation
    Teuchos::RCP<Core::FE::Discretization> thrdisc_;
    // our solution
    //! global temperature DOFs
    Teuchos::RCP<Core::LinAlg::Vector<double>> temp_;
    //! global temperature rate DOFs
    Teuchos::RCP<Core::LinAlg::Vector<double>> rate_;
    //! global temperature DOFs
    Teuchos::RCP<Core::LinAlg::Vector<double>> flux_;
    //! NOTE: these have to be present explicitly
    //! as they are not part of the problem instance like in fluid3

  };  // Core::Utils::ResultTest

}  // namespace Thermo


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
