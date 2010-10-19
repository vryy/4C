/*----------------------------------------------------------------------*/
/*!
\file fluidxfluidresulttest.H

\brief tesing of fluid calculation results

<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef FLUIDXFLUIDRESULTTEST_H_
#define FLUIDXFLUIDRESULTTEST_H_

#ifdef CCADISCRET

#include "../drt_lib/drt_resulttest.H"
#include "fluidxfluidimplicitintegration.H"
#include "fluid_genalpha_integration.H"

namespace FLD
{
  
class FluidXFluidResultTest : public DRT::ResultTest
{
public:

    /*!
    \brief not documented yet
    */
    FluidXFluidResultTest(FluidXFluidImplicitTimeInt& fluid);
    
    /// our version of nodal value tests
    /*!
    Possible position flags are "velx", "vely", "velz" and
    "pressure". With the obvious meaning.
    */
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    /// return true if this is a fluid test
    virtual bool Match(DRT::INPUT::LineDefinition& res);

private:

    /// pointer to fluid discretization
    RefCountPtr<DRT::Discretization> fluiddis_;
    /// pointer to unknown vector with values
    RefCountPtr<Epetra_Vector>       mysol_;
    /// pointer to traction vector with values
    RefCountPtr<Epetra_Vector>       mytraction_;
};

}

#endif /* CCADISCRET       */

#endif /* FLUIDXFLUIDRESULTTEST_H_ */
