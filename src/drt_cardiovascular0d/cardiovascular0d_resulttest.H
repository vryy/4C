/*!----------------------------------------------------------------------
\file cardiovascular0d_resulttest.H

\brief Monolithic coupling of 3D structural dynamics and 0D cardiovascular flow models

\level 2

<pre>
\maintainer Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>
*----------------------------------------------------------------------*/

#ifndef CARDIOVASCULAR0DRESULTTEST_H
#define CARDIOVASCULAR0DRESULTTEST_H

#include <Epetra_Vector.h>
#include "../drt_lib/drt_resulttest.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace UTILS
{
  class Cardiovascular0D;
  class Cardiovascular0DManager;
}

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

  Cardiovascular0DResultTest(UTILS::Cardiovascular0DManager& cardvasc0dman, Teuchos::RCP<DRT::Discretization> discr);

  virtual void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);



private:

  Teuchos::RCP<DRT::Discretization> actdisc_; ///< standard discretization

  Teuchos::RCP<Epetra_Vector> cardvasc0d_dof_;

  bool havecardio_4elementwindkessel_;
  bool havecardio_arterialproxdist_;
  bool havecardio_syspulcirculation_;

};

#endif /* #ifndef CARDIOVASCULAR0DRESULTTEST_H */
