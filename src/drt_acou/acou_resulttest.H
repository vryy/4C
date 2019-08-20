/*----------------------------------------------------------------------*/
/*! \file

\brief testing of acoustical calculation results

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
*/
/*----------------------------------------------------------------------*/

#ifndef ACOU_RESULTTEST_H
#define ACOU_RESULTTEST_H

#include "../drt_lib/drt_resulttest.H"
#include <Epetra_Vector.h>

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace ACOU
{
  class AcouTimeInt;

  class AcouResultTest : public DRT::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    AcouResultTest(AcouTimeInt& acoualgo);


    /// our version of nodal value tests
    /*!
      Possible position flags is only "pressure"
     */
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

   private:
    /// Teuchos::RCP to acoustical discretization
    Teuchos::RCP<DRT::Discretization> dis_;
    /// Teuchos::RCP to solution vector
    Teuchos::RCP<Epetra_Vector> mysol_;

  };  // class AcouResultTest

}  // namespace ACOU

#endif /* ACOU_RESULTTEST_H */
