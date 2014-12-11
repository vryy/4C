/*----------------------------------------------------------------------*/
/*!
\file xfluidresulttest.H

\brief tesing of fluid calculation results

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef XFLUIDFLUIDRESULTTEST_H_
#define XFLUIDFLUIDRESULTTEST_H_


#include "../drt_lib/drt_resulttest.H"
#include <Epetra_Vector.h>

namespace DRT
{
class Discretization;
}

namespace FLD
{
// forward declarations
class XFluid;

/*!
  ResultTest class for XFluid
 */
class XFluidResultTest : public DRT::ResultTest
{
public:

  XFluidResultTest( XFluid * xfluid );

  /// our version of nodal value tests
  /*!
    Possible position flags are "velx", "vely", "velz" and
    "pressure". With the obvious meaning.
   */
  virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

private:
  DRT::Discretization & discret_;
  Teuchos::RCP<Epetra_Vector> velnp_;

};

}

#endif /* XFLUIDRESULTTEST_H_ */
