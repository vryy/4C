/*--------------------------------------------------------------------------*/
/*!
\file lubrication_resulttest.H

\brief testing of lubrication calculation results

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_LUBRICATION_LUBRICATION_RESULTTEST_H_
#define SRC_DRT_LUBRICATION_LUBRICATION_RESULTTEST_H_


#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>

#include "../drt_lib/drt_resulttest.H"

namespace DRT
{
  class Discretization;
  class Node;
}  // namespace DRT

namespace LUBRICATION
{
  // forward declaration
  class TimIntImpl;

  /*!
    \brief lubrication specific result test class

    \author wirtz
    \date 11/15
  */
  class ResultTest : public DRT::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    ResultTest(Teuchos::RCP<TimIntImpl> lubrication);


    /// our version of nodal value tests
    /*!
      Possible position flags is only "pre"
     */
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    //! test special quantity not associated with a particular element or node
    virtual void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

   protected:
    //! get nodal result to be tested
    double ResultNode(const std::string quantity,  //! name of quantity to be tested
        DRT::Node* node                            //! node carrying the result to be tested
        ) const;

    //! get special result to be tested
    virtual double ResultSpecial(const std::string quantity  //! name of quantity to be tested
        ) const;

   private:
    /// Teuchos::RCP to lubrication discretization
    Teuchos::RCP<DRT::Discretization> dis_;
    /// Teuchos::RCP to solution vector
    Teuchos::RCP<Epetra_Vector> mysol_;
    /// number of iterations in last newton iteration
    int mynumiter_;
  };
}  // namespace LUBRICATION


#endif /* SRC_DRT_LUBRICATION_LUBRICATION_RESULTTEST_H_ */
