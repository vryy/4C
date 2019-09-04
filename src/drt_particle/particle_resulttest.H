/*----------------------------------------------------------------------*/
/*! \file

\brief testing of particle calculation results

\level 2

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/

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
  //!                             "radius"
  virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

  //! test special quantity not associated with a particular element or node
  void TestSpecial(
      DRT::INPUT::LineDefinition& res,  //!< input file line containing result test specification
      int& nerr,                        //!< number of failed result tests
      int& test_count                   ///< number of result tests
  );

 private:
  //! get special result to be tested
  double ResultSpecial(const std::string& quantity  //! name of quantity to be tested
      ) const;

  //! our discretisation
  Teuchos::RCP<DRT::Discretization> partdisc_;
  // our solution
  //! global displacement DOFs
  Teuchos::RCP<const Epetra_Vector> dis_;
  //! global velocity DOFs
  Teuchos::RCP<const Epetra_Vector> vel_;
  //! global acceleration DOFs
  Teuchos::RCP<const Epetra_Vector> acc_;
  //! global angular velocity DOFs
  Teuchos::RCP<const Epetra_Vector> angvel_;
  //! global angular acceleration DOFs
  Teuchos::RCP<const Epetra_Vector> angacc_;
  //! global radius
  Teuchos::RCP<const Epetra_Vector> radius_;
  //! kinetic energy
  const double kineticenergy_;
  //! internal energy
  const double internalenergy_;
  //! external energy
  const double externalenergy_;
  //! maximum particle-particle or particle-wall penetration
  const double maxpenetration_;
};

#endif /* #ifndef PARTICLERESULTTEST_H */
