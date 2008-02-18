#ifdef CCADISCRET

#ifndef FSI_ALE_RESULTTEST_H
#define FSI_ALE_RESULTTEST_H

#include "../drt_lib/drt_resulttest.H"
#include "fsi_ale.H"

namespace FSI {

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
  virtual void TestNode(RESULTDESCR* res, int& nerr, int& test_count);

  /// return true if this is a ale test
  virtual bool Match(RESULTDESCR* res);

private:

  /// our ale algorithmus
  AleLinear& ale_;

};

}

#endif
#endif
