/*----------------------------------------------------------------------*/
/*!
\file xfluidresulttest.H

\brief tesing of fluid calculation results

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef XFLUIDRESULTTEST_H
#define XFLUIDRESULTTEST_H

#ifdef CCADISCRET

#include "../drt_lib/drt_resulttest.H"
#include "xfluidimplicitintegration.H"

namespace FLD
{

/*!
  \brief XFluid specific result test class

  \author u.kue
*/
class XFluidResultTest : public DRT::ResultTest
{
public:

  /*!
  \brief not documented yet
  */
  explicit XFluidResultTest(XFluidImplicitTimeInt& fluid);

  /// our version of nodal value tests
  /*!
    Possible position flags are "velx", "vely", "velz" and
    "pressure". With the obvious meaning.
   */
  virtual void TestNode(const RESULTDESCR* res, int& nerr, int& test_count);

  /// return true if this is a fluid test
  virtual bool Match(const RESULTDESCR* res);

private:

  /// pointer to fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis_;
  /// pointer to unknown vector with values 
  Teuchos::RCP<Epetra_Vector>       mysol_;
};

}

#endif /* CCADISCRET       */

#endif
