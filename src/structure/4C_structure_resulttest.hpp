/*----------------------------------------------------------------------*/
/*! \file
\brief testing of structure calculation results


\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_RESULTTEST_HPP
#define FOUR_C_STRUCTURE_RESULTTEST_HPP

#include "4C_config.hpp"

#include "4C_lib_resulttest.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace STR
{
  class TimInt;
}

/*!
  \brief Structure specific result test class
*/
class StruResultTest : public DRT::ResultTest
{
 public:
  //! Constructor for time integrators of general kind
  //! \author bborn \date 06/08
  StruResultTest(STR::TimInt& tintegrator);

  /*!
  \brief Test nodal values in solid/structure field and time integration

  Possible position flags are
  - displacements "dispx", "dispy", "dispz"
  - velocitites "velx", "vely", "velz",
  - accelerations "accx", "accy", "accz"
  */
  void TestNode(INPUT::LineDefinition& res, int& nerr, int& test_count) final;

  /*!
  \brief Test special quantities of structure discretization / time integration

  See GetSpecialResultForTesting() for a list and extraction mechanisms of special results to be
  tested.
  */
  void TestSpecial(INPUT::LineDefinition& res, int& nerr, int& test_count) final;

 private:
  //! Extract specified quantity from the structure field
  double GetSpecialResultForTesting(const std::string& quantity);

  //! Structure time integrator
  Teuchos::RCP<STR::TimInt> timeintegrator_;

  //! Structure discretisation
  Teuchos::RCP<DRT::Discretization> strudisc_;

  //! @name Solution
  //!@{

  //! global displacement DOFs
  Teuchos::RCP<const Epetra_Vector> dis_;

  //! global material displacement DOFs
  Teuchos::RCP<const Epetra_Vector> dism_;

  //! global velocity DOFs
  Teuchos::RCP<const Epetra_Vector> vel_;

  //! global acceleration DOFs
  Teuchos::RCP<const Epetra_Vector> acc_;

  //!@}
};

FOUR_C_NAMESPACE_CLOSE

#endif
