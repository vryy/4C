/*----------------------------------------------------------------------*/
/*!
\file fluidresulttest.H

\brief testing of fluid calculation results

\level 1

<pre>
\maintainer Benjamin Krank & Martin Kronbichler
             {krank,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15252/-235
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

    /// our version of nodal value tests
    /*!
    Possible position flags are "velx", "vely", "velz", "pressure",
    "tractionx", "tractiony", "tractionz", "errvel", "errpre" and "divu".
    With the obvious meaning.
    */
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    /// our version of element value tests
    /*!
    Possible flag is "fluidfraction". With the obvious meaning.
    */
    virtual void TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

private:

    /// pointer to fluid discretization
    Teuchos::RCP<DRT::Discretization> fluiddis_;
    /// pointer to unknown vector with nodal values
    Teuchos::RCP<Epetra_Vector>       mysol_;
    /// pointer to unknown vector with values at nodes
    Teuchos::RCP<Epetra_Vector>       myvan_;
    /// pointer to traction vector with values
    Teuchos::RCP<Epetra_Vector>       mytraction_;
    /// pointer to traction vector with values
    Teuchos::RCP<Epetra_Vector>       mywss_;
    /// pointer to error evaluation
    Teuchos::RCP<std::vector <double> > myerror_;
    /// pointer to div u evaluation
    Teuchos::RCP<double> mydivu_;
    /// pointer to density scaling values
    Teuchos::RCP<Epetra_Vector>       mydensity_scaling_;

};

}


#endif
