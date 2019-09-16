/*----------------------------------------------------------------------*/
/*! \file

\brief testing of electromagnetic calculation results

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            089 - 289-15244
*/
/*----------------------------------------------------------------------*/

#ifndef ELEMAG_RESULTTEST_H
#define ELEMAG_RESULTTEST_H

#include "../drt_lib/drt_resulttest.H"
#include <Epetra_Vector.h>
#include <Epetra_SerialDenseVector.h>

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace ELEMAG
{
  class ElemagTimeInt;

  class ElemagResultTest : public DRT::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    ElemagResultTest(ElemagTimeInt& elemagalgo);


    /// nodal value tests
    /*!
      Possible position flags is only "pressure"
     */
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

   private:
    /// Teuchos::RCP to elemagstical discretization
    Teuchos::RCP<DRT::Discretization> dis_;
    /// Teuchos::RCP to solution vector
    Teuchos::RCP<Epetra_Vector> mysol_;
    /// Error vector
    Teuchos::RCP<Epetra_SerialDenseVector> error_;

  };  // class ElemagResultTest

}  // namespace ELEMAG

#endif /* ELEMAG_RESULTTEST_H */