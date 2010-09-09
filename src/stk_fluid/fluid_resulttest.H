#ifdef STKADAPTIVE

#ifndef FLUID_RESULTTEST_H
#define FLUID_RESULTTEST_H

#include "../drt_lib/drt_resulttest.H"
#include "fluid_implicit.H"

namespace STK
{
  class Mesh;

  namespace FLD
  {

/*!
\brief Fluid specific result test class

\author u.kue
*/
class FluidResultTest : public DRT::ResultTest
{
public:

    FluidResultTest(Fluid& fluid);

    /// our version of nodal value tests
    /*!
    Possible position flags are "velx", "vely", "velz" and
    "pressure". With the obvious meaning.
    */
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    /// return true if this is a fluid test
    virtual bool Match(DRT::INPUT::LineDefinition& res);

private:

  Mesh & mesh_;
};

  }
}

#endif
#endif
