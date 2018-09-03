/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_resulttest.H

 \brief result test for multiphase porous fluid

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_RESULTTEST_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_RESULTTEST_H_



#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>

#include "../drt_lib/drt_resulttest.H"

// forward declaration
namespace DRT
{
  class Discretization;
  class Node;
}  // namespace DRT

namespace POROFLUIDMULTIPHASE
{
  // forward declaration
  class TimIntImpl;

  /*!
    \brief POROFLUIDMULTIPHASE specific result test class

    \author vuong
    \date 08/16
  */
  class ResultTest : public DRT::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    ResultTest(TimIntImpl& porotimint);


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
    //! time integrator
    const TimIntImpl& porotimint_;
  };
}  // namespace POROFLUIDMULTIPHASE


#endif /* SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASERESULTTEST_H_ */
