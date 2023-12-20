/*----------------------------------------------------------------------*/
/*! \file

\brief Monolithic coupling of 3D structural dynamics and 0D cardiovascular flow models

\level 2

*----------------------------------------------------------------------*/

#ifndef BACI_CARDIOVASCULAR0D_RESULTTEST_H
#define BACI_CARDIOVASCULAR0D_RESULTTEST_H

#include "baci_config.H"

#include "baci_lib_resulttest.H"

#include <Epetra_Vector.h>

BACI_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace UTILS
{
  class Cardiovascular0D;
  class Cardiovascular0DManager;
}  // namespace UTILS

namespace IO
{
  class DiscretizationWriter;
}

/*!
  \brief Structure specific result test class
*/
class Cardiovascular0DResultTest : public DRT::ResultTest
{
 public:
  Cardiovascular0DResultTest(
      UTILS::Cardiovascular0DManager& cardvasc0dman, Teuchos::RCP<DRT::Discretization> discr);

  void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;



 private:
  Teuchos::RCP<DRT::Discretization> actdisc_;  ///< standard discretization

  const Teuchos::RCP<Epetra_Vector> cardvasc0d_dof_;

  const bool havecardio_4elementwindkessel_;
  const bool havecardio_arterialproxdist_;
  const bool havecardio_syspulcirculation_;
  const bool havecardiorespir_syspulperiphcirculation_;
};

BACI_NAMESPACE_CLOSE

#endif  // CARDIOVASCULAR0D_RESULTTEST_H
