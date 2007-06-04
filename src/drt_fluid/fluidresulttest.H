/*----------------------------------------------------------------------*/
/*!
\file fluidresulttest.H

\brief tesing of fluid calculation results

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUIDRESULTTEST_H
#define FLUIDRESULTTEST_H

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifdef RESULTTEST
#ifdef D_FLUID

#include "drt_resulttest.H"
#include "fluidimplicitintegration.H"

/*!
  \class FluidResultTest

  \brief Fluid specific result test class

  \author u.kue
*/
class FluidResultTest : public DRT::ResultTest
{
public:

  /*!
  \brief not documented yet
  */
  FluidResultTest(FluidImplicitTimeInt& fluid);

  /// our version of nodal value tests
  /*!
    Possible position flags are "velx", "vely", "velz" and
    "pressure". With the obvious meaning.
   */
  virtual void TestNode(RESULTDESCR* res, int& nerr, int& test_count);

  /// return true if this is a fluid test
  virtual bool Match(RESULTDESCR* res);

private:

  /// our fluid algorithmus
  FluidImplicitTimeInt& fluid_;
};


#endif
#endif
#endif /* TRILINOS_PACKAGE */
#endif /* CCADISCRET       */

#endif
