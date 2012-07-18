
/*----------------------------------------------------------------------*/
/*!
\file xfluidfluidresulttest.H

\brief tesing of fluid calculation results

<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef XFLUIDFLUIDRESULTTEST_H_
#define XFLUIDFLUIDRESULTTEST_H_


#include "../drt_lib/drt_resulttest.H"
#include <Epetra_Vector.h>

namespace DRT
{
class Discretization;
}

namespace FLD
{
// forward declarations
class XFluidFluid;


class XFluidFluidResultTest : public DRT::ResultTest
{
public:

    /*!
    \brief not documented yet
    */
    XFluidFluidResultTest(XFluidFluid& fluid);

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
    Teuchos::RCP<DRT::Discretization> embfluiddis_;
    /// pointer to xfluid discretization
    Teuchos::RCP<DRT::Discretization> bgfluiddis_;
    /// pointer to unknown vector with values
    Teuchos::RCP<Epetra_Vector>       embfluidsol_;
    /// pointer to unknown vector with values
    Teuchos::RCP<Epetra_Vector>       bgfluidsol_;
    /// pointer to traction vector with values
    //Teuchos::RCP<Epetra_Vector>       mytraction_;
};

}


#endif /* XFLUIDFLUIDRESULTTEST_H_ */
