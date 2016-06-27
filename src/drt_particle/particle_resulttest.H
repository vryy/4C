/*----------------------------------------------------------------------*/
/*!
\file particle_resulttest.H

\brief testing of particle calculation results

\level 2

\maintainer Georg Hammerl
*/
/*----------------------------------------------------------------------*/

#ifndef PARTICLERESULTTEST_H
#define PARTICLERESULTTEST_H

#include <Epetra_Vector.h>
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

namespace PARTICLE
{
  class TimInt;
}

/*!
  \brief Structure specific result test class
*/
class PartResultTest : public DRT::ResultTest
{
public:
  //! Constructor for particle time integrators of general kind
  PartResultTest(PARTICLE::TimInt& tintegrator);

  //! \brief particle version of nodal value tests
  //!
  //! Possible position flags are "dispx", "dispy", "dispz",
  //!                             "velx", "vely", "velz",
  //!                             "accx", "accy", "accz"
  //!                             "radius", "temperature"
  virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

private:

  //! our discretisation
  Teuchos::RCP<DRT::Discretization> partdisc_;
  // our solution
  //! global displacement DOFs
  Teuchos::RCP<const Epetra_Vector> dis_;
  //! global velocity DOFs
  Teuchos::RCP<const Epetra_Vector> vel_;
  //! global acceleration DOFs
  Teuchos::RCP<const Epetra_Vector> acc_;
  //! global radius DOFs
  Teuchos::RCP<const Epetra_Vector> radius_;
  //! global temperature DOFs
  Teuchos::RCP<const Epetra_Vector> temperature_;
};

#endif /* #ifndef PARTICLERESULTTEST_H */
