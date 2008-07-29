/*----------------------------------------------------------------------*/
/*!
\file scatra_resulttest.H

\brief tesing of scalar transport calculation results

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_RESULTTEST_H
#define SCATRA_RESULTTEST_H

#ifdef CCADISCRET

#include "../drt_lib/drt_resulttest.H"
#include "scatra_timint_implicit.H"


namespace SCATRA{

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
  ScaTraResultTest(ScaTraImplicitTimeInt& scatra);


  /// our version of nodal value tests
  /*!
    Possible position flags is only "phi"
   */
  virtual void TestNode(const RESULTDESCR* res, int& nerr, int& test_count);

  /// return true if this is a condif test
  virtual bool Match(const RESULTDESCR* res);

private:

  /// RCP to scalar transport discretization
  RefCountPtr<DRT::Discretization> dis_;
  /// RCP to solution vector
  RefCountPtr<Epetra_Vector>       mysol_;
  /// RCP to mass/heat flux vector
  RefCountPtr<Epetra_MultiVector>  myflux_;
};

} // namespace SCATRA

#endif /* CCADISCRET       */

#endif
