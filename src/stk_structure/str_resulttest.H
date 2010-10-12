#ifdef STKADAPTIVE

#ifndef STRUCTURE_RESULTTEST_H
#define STRUCTURE_RESULTTEST_H

#include "../drt_lib/drt_resulttest.H"
#include "str_genalpha.H"

namespace STK
{
  class Mesh;

  namespace STR
  {

/*!
\brief Structure specific result test class

\author u.kue
*/
class StructureResultTest : public DRT::ResultTest
{
public:

    StructureResultTest(Structure& structure);

    /// our version of nodal value tests
    /*!
    Possible position flags are "velx", "vely", "velz" and
    "pressure". With the obvious meaning.
    */
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    /// return true if this is a structure test
    virtual bool Match(DRT::INPUT::LineDefinition& res);

private:

  Mesh & mesh_;
};

  }
}

#endif
#endif
