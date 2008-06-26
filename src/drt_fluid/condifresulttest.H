/*----------------------------------------------------------------------*/
/*!
\file condifresulttest.H

\brief tesing of convection-diffusion calculation results

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef CONDIFRESULTTEST_H
#define CONDIFRESULTTEST_H

#ifdef CCADISCRET

#include "../drt_lib/drt_resulttest.H"
#include "condifimplicitintegration.H"
#include "condif_genalpha_integration.H"

/*!
  \brief convection-diffusion specific result test class

  \author gjb
*/
class ConDifResultTest : public DRT::ResultTest
{
public:

  /*!
  \brief not documented yet
  */
  ConDifResultTest(CondifImplicitTimeInt& condif);

  /*!
  \brief not documented yet
  */
  ConDifResultTest(CondifGenAlphaIntegration& condif);

  /// our version of nodal value tests
  /*!
    Possible position flags is only "phi"
   */
  virtual void TestNode(const RESULTDESCR* res, int& nerr, int& test_count);

  /// return true if this is a condif test
  virtual bool Match(const RESULTDESCR* res);

private:

  RefCountPtr<DRT::Discretization> condifdis_;
  RefCountPtr<Epetra_Vector>       mysol_;
  RefCountPtr<Epetra_MultiVector>  myflux_;
};


#endif /* CCADISCRET       */

#endif
