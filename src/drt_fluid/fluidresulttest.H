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


#include "../drt_lib/drt_resulttest.H"
#include <Epetra_Vector.h>

namespace DRT
{
class Discretization;
}

namespace FLD
{
// forward declarations
class FluidImplicitTimeInt;
class FluidGenAlphaIntegration;


/*!
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

    /*!
    \brief not documented yet
    */
    FluidResultTest(FluidGenAlphaIntegration& fluid);

    /// our version of nodal value tests
    /*!
    Possible position flags are "velx", "vely", "velz" and
    "pressure". With the obvious meaning.
    */
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

private:

    /// pointer to fluid discretization
    Teuchos::RCP<DRT::Discretization> fluiddis_;
    /// pointer to unknown vector with values
    Teuchos::RCP<Epetra_Vector>       mysol_;
    /// pointer to traction vector with values
    Teuchos::RCP<Epetra_Vector>       mytraction_;
};

}


#endif
