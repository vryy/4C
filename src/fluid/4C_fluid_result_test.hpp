/*-----------------------------------------------------------*/
/*! \file

\brief testing of fluid calculation results


\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_RESULT_TEST_HPP
#define FOUR_C_FLUID_RESULT_TEST_HPP


#include "4C_config.hpp"

#include "4C_utils_result_test.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace FLD
{
  // forward declarations
  class FluidImplicitTimeInt;
  class FluidGenAlphaIntegration;


  /*!
  \brief Fluid specific result test class

  \author u.kue
  */
  class FluidResultTest : public Core::UTILS::ResultTest
  {
   public:
    /*!
    \brief not documented yet
    */
    FluidResultTest(FluidImplicitTimeInt& fluid);

    /// our version of nodal value tests
    /*!
    Possible position flags are "velx", "vely", "velz", "pressure",
    "tractionx", "tractiony", "tractionz", "errvel", "errpre" and "divu".
    With the obvious meaning.
    */
    void test_node(
        const Core::IO::InputParameterContainer& container, int& nerr, int& test_count) override;

   private:
    /// pointer to fluid discretization
    Teuchos::RCP<Core::FE::Discretization> fluiddis_;
    /// pointer to unknown vector with nodal values
    Teuchos::RCP<Epetra_Vector> mysol_;
    /// pointer to traction vector with values
    Teuchos::RCP<Epetra_Vector> mytraction_;
    /// pointer to traction vector with values
    Teuchos::RCP<Epetra_Vector> mywss_;
    /// pointer to error evaluation
    Teuchos::RCP<std::vector<double>> myerror_;
    /// pointer to div u evaluation
    Teuchos::RCP<double> mydivu_;
  };

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
