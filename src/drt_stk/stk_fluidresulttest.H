#if 0
#ifdef STKADAPTIVE

#ifndef STK_FLUIDRESULTTEST_H
#define STK_FLUIDRESULTTEST_H

#include "../drt_lib/drt_resulttest.H"
#include "stk_fluid.H"

namespace STK
{

  class Mesh;

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

#endif
#endif
#endif
