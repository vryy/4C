/*----------------------------------------------------------------------*/
/*! \file

\brief testing of artery calculation results


\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ART_NET_ARTERY_RESULTTEST_HPP
#define FOUR_C_ART_NET_ARTERY_RESULTTEST_HPP

#include "baci_config.hpp"

#include "baci_lib_resulttest.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

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
    void TestNode(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

    /// our version of element value tests
    void TestElement(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    /// Teuchos::RCP to scalar transport discretization
    Teuchos::RCP<DRT::Discretization> dis_;
    /// Teuchos::RCP to solution vector
    Teuchos::RCP<const Epetra_Vector> mysol_;
    /// Teuchos::RCP to element volumetric flow
    Teuchos::RCP<const Epetra_Vector> myelevolflow_;
    /// Teuchos::RCP to element radius
    Teuchos::RCP<const Epetra_Vector> myeleradius_;
  };

}  // namespace ART

BACI_NAMESPACE_CLOSE

#endif
