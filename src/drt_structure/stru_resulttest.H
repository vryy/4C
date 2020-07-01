/*----------------------------------------------------------------------*/
/*! \file
\brief testing of structure calculation results


\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef STRURESULTTEST_H
#define STRURESULTTEST_H

#include "../drt_lib/drt_resulttest.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class Solver;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace STR
{
  class TimInt;
}

class Epetra_Vector;

/*!
  \brief Structure specific result test class
*/
class StruResultTest : public DRT::ResultTest
{
 public:
  //! Constructor for time integrators of general kind
  //! \author bborn \date 06/08
  StruResultTest(STR::TimInt& tintegrator);

  //! \brief structure version of nodal value tests
  //!
  //! Possible position flags are "dispx", "dispy", "dispz",
  //!                             "velx", "vely", "velz",
  //!                             "accx", "accy", "accz"
  virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

 private:
  //! our discretisation
  Teuchos::RCP<DRT::Discretization> strudisc_;
  // our solution
  //! global displacement DOFs
  Teuchos::RCP<const Epetra_Vector> dis_;
  //! global material displacement DOFs
  Teuchos::RCP<const Epetra_Vector> dism_;
  //! global velocity DOFs
  Teuchos::RCP<const Epetra_Vector> vel_;
  //! global acceleration DOFs
  Teuchos::RCP<const Epetra_Vector> acc_;
  /* NOTE: these have to be present explicitly
   * as they are not part of the problem instance like in fluid3
   */
};

#endif /* #ifndef STRURESULTTEST_H */
