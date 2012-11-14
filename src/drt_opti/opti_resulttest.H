/*!------------------------------------------------------------------------------------------------*
\file opti_resulttest.H

\brief Result test for optimization algorithms

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef OPTI_RESULTTEST_H_
#define OPTI_RESULTTEST_H_


#include "../drt_lib/drt_resulttest.H"


class Epetra_Vector;

namespace DRT
{
  class Discretization;
}

namespace OPTI
{
  class GCMMA;

/*!
\brief Optimization specific result test class

\author winklmaier
*/
class OptiResultTest : public DRT::ResultTest
{
public:

    /*!
    \brief constructor
    */
  OptiResultTest(const OPTI::GCMMA& optimizer);


    /// our version of nodal value tests
    /*!
    Possible position flags are "velx", "vely", "velz" and
    "pressure". With the obvious meaning.
    */
    virtual void TestNode(
        DRT::INPUT::LineDefinition& res,
        int& nerr,
        int& test_count
    );

private:

    /// pointer to fluid discretization
    const Teuchos::RCP<const DRT::Discretization> optidis_;
    /// pointer to vector with values
    const Teuchos::RCP<const Epetra_Vector> sol_;
};
}


#endif /* OPTI_RESULTTEST_H_ */
