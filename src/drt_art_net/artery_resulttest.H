/*----------------------------------------------------------------------*/
/*!
\file artery_resulttest.H

\brief testing of artery calculation results

\maintainer Lena Yoshihara

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef ART_NET_RESULTTEST_H
#define ART_NET_RESULTTEST_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>

#include "../drt_lib/drt_resulttest.H"

namespace DRT
{
  class Discretization;
}

namespace ART
{
  // forward declaration
  class ArtNetExplicitTimeInt;
  class ArtNetImplStationary;

  /*!
    \brief artnet specific result test class

    \author Mahmoud Ismail
    \date 11/11
  */
  class ArteryResultTest : public DRT::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    ArteryResultTest(ArtNetExplicitTimeInt& art_net);

    /*!
    \brief constructor
    */
    ArteryResultTest(ArtNetImplStationary& art_net);


    /// our version of nodal value tests
    /*!
      Possible position flags is only "phi"
     */
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

   private:
    /// Teuchos::RCP to scalar transport discretization
    Teuchos::RCP<DRT::Discretization> dis_;
    /// Teuchos::RCP to solution vector
    Teuchos::RCP<const Epetra_Vector> mysol_;
  };

}  // namespace ART

#endif  // ART_NET_RESULTTEST_H
