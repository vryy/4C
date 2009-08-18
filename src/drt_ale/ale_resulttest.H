/*----------------------------------------------------------------------*/
/*!
\file ale_resulttest.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ALE_RESULTTEST_H
#define ALE_RESULTTEST_H

#include "../drt_lib/drt_resulttest.H"
#include "ale.H"
#include "ale_lin.H"
#include "ale_springs.H"


namespace ALE {

class AleLinear;

/// Result test subclass for linear ale algorithm
class AleResultTest : public DRT::ResultTest
{
public:
  AleResultTest(AleLinear& ale);

  /// our version of nodal value tests
  /*!
    Possible position flags are "velx", "vely", "velz" and
    "pressure". With the obvious meaning.
   */
  virtual void TestNode(const RESULTDESCR* res, int& nerr, int& test_count);

  /// return true if this is a ale test
  virtual bool Match(const RESULTDESCR* res);

private:

  /// our ale algorithmus
  AleLinear& ale_;

};

}

#endif
#endif
