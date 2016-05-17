/*----------------------------------------------------------------------*/
/*!
\file scatra_resulttest.H

\brief testing of scalar transport calculation results

<pre>
\level 1

\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_RESULTTEST_H
#define SCATRA_RESULTTEST_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>

#include "../drt_lib/drt_resulttest.H"

namespace DRT
{
  class Discretization;
  class Node;
}

namespace SCATRA{

// forward declaration
class ScaTraTimIntImpl;

/*!
  \brief scalar-transport specific result test class

  \author gjb
  \date 07/08
*/
class ScaTraResultTest : public DRT::ResultTest
{
public:

  /*!
  \brief constructor
  */
  ScaTraResultTest(Teuchos::RCP<ScaTraTimIntImpl> scatratimint);


  /// our version of nodal value tests
  /*!
    Possible position flags is only "phi"
   */
  virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

  //! test special quantity not associated with a particular element or node
  virtual void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

protected:

  //! get nodal result to be tested
  virtual double ResultNode(
      const std::string   quantity,   //! name of quantity to be tested
      DRT::Node*          node        //! node carrying the result to be tested
      ) const;

  //! get special result to be tested
  virtual double ResultSpecial(
      const std::string   quantity   //! name of quantity to be tested
      ) const;

  //! time integrator
  const Teuchos::RCP<const ScaTraTimIntImpl> scatratimint_;

private:
};
} // namespace SCATRA
#endif
