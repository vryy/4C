/*----------------------------------------------------------------------*/
/*!
\file particleMeshFree_rendering_resulttest.H

\brief testing of meshfree rendering calculation results

\level 2

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                             sfuchs 06/17 |
 *----------------------------------------------------------------------*/
#ifndef PARTICLEMESHFREERENDERINGRESULTTEST_H
#define PARTICLEMESHFREERENDERINGRESULTTEST_H

/*----------------------------------------------------------------------*
 | headers                                                 sfuchs 06/17 |
 *----------------------------------------------------------------------*/
#include <Epetra_Vector.h>
#include "../drt_lib/drt_resulttest.H"

/*----------------------------------------------------------------------*
 | forward declarations                                    sfuchs 06/17 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
} // namespace DRT

namespace PARTICLE
{
  class Rendering;
} // namespace PARTICLE

/*----------------------------------------------------------------------*
 | meshfree rendering result test class                    sfuchs 06/17 |
 *----------------------------------------------------------------------*/
class ParticleMeshfreeRenderingResultTest : public DRT::ResultTest
{
public:
  //! Constructor for particle time integrators of general kind
  ParticleMeshfreeRenderingResultTest(PARTICLE::Rendering& rendering);

  //! nodal value tests
  //!
  //! Possible position flags are "velx", "vely", "velz",
  //!                             "accx", "accy", "accz",
  //!                             "density", "specEnthalpy", "temperature", "pressure"
  virtual void TestNode(
      DRT::INPUT::LineDefinition&   res,         //!< input file line containing result test specification
      int&                          nerr,        //!< number of failed result tests
      int&                          test_count   ///< number of result tests
      );

private:

  //! rendering discretisation
  Teuchos::RCP<const DRT::Discretization> discret_;

  //! rendering solution vectors
  Teuchos::RCP<const Epetra_Vector> vel_; //! rendering velocity

  Teuchos::RCP<const Epetra_Vector> acc_; //! rendering acceleration

  Teuchos::RCP<const Epetra_Vector> density_; //! rendering density

  Teuchos::RCP<const Epetra_Vector> specEnthalpy_; //! rendering specific enthalpy

  Teuchos::RCP<const Epetra_Vector> temperature_; //! rendering temperature

  Teuchos::RCP<const Epetra_Vector> pressure_; //! rendering pressure

};

/*----------------------------------------------------------------------*/
#endif // PARTICLEMESHFREERENDERINGRESULTTEST_H
