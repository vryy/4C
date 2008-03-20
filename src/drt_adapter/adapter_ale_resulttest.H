#ifdef CCADISCRET

#ifndef ADAPTER_ALE_RESULTTEST_H
#define ADAPTER_ALE_RESULTTEST_H

#include "../drt_lib/drt_resulttest.H"
#include "../drt_adapter/adapter_ale.H"

namespace ADAPTER {

class AleLinear;

/// Result test subclass for linear ale algorithm
class AleResultTest : public DRT::ResultTest
{
public:
  AleResultTest(ADAPTER::AleLinear& ale);

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
