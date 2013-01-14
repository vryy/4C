/*----------------------------------------------------------------------*/
/*!
\file fsi_resulttest.H

\brief testing of fsi specific calculation results

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-15262
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FSIRESULTTEST_H
#define FSIRESULTTEST_H

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include "../drt_lib/drt_resulttest.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace FSI
{
  class Monolithic;

  class MonolithicNoNOX;

  /*!
    \brief FSI specific result test class

    Here, additional result tests of quantities that do not belong to a single
    field are tested. Basically, this should work for monolithic and partitioned
    schemes as well.

    Feel free to add further testing functionalities!

    \sa ResultTest
    \author mayr.mt
    \date 11/2012
  */
  class FSIResultTest : public DRT::ResultTest
  {
  public:
    //! constructor
    FSIResultTest(
      Teuchos::RCP<FSI::Monolithic> fsi,
      const Teuchos::ParameterList& fsidyn
    );

    FSIResultTest(
      Teuchos::RCP<FSI::MonolithicNoNOX> fsi,
      const Teuchos::ParameterList& fsidyn
    );

    //! \brief fsi version of nodal value tests
    //!
    //! Possible position flags are "lambdax", "lambday", "lambdaz"
    virtual void TestNode(DRT::INPUT::LineDefinition& res,  ///< the line definition from input file
                          int& nerr,                        ///< number of tests with errors
                          int& test_count);                 ///< number of tests performed

  private:

    //! slave discretisation
    Teuchos::RCP<DRT::Discretization> slavedisc_;
    //! Lagrange multiplier living on the slave discretization
    Teuchos::RCP<Epetra_Vector> fsilambda_;
  };
}
#endif /* #ifndef FSIRESULTTEST_H */
