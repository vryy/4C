/*----------------------------------------------------------------------*/
/*!
\file invana_resulttest.H

\brief tesing of inverse analysis results

<pre>
Maintainer: Sebastian Kehl
kehl@lnm.mw.tum.de
089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef INVANA_RESULTTEST_H
#define INVANA_RESULTTEST_H


#include "../drt_lib/drt_resulttest.H"
#include <Epetra_Vector.h>

namespace DRT
{
class Discretization;
}

namespace STR
{
namespace INVANA
{
// forward declarations
class StatInvAnalysis;
class MatParManager;


/*!
\brief Invana specific result test class

\author keh
*/
class InvAnaResultTest : public DRT::ResultTest
{
public:

    /*!
    \brief not documented yet
    */
    InvAnaResultTest(StatInvAnalysis& ia);

    /*!
    Possible flag is the name of the material parameter as used by the material, e.g AAAstopro: YOUNG or BETA
    */
    virtual void TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    /*
    Possible flags are "gradient" and "error"
    */

    virtual void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

private:

    //! pointer to discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! pointer to unknown vector with values
    Teuchos::RCP<Epetra_MultiVector> mysol_;

    //! pointer to manager of the material parameters to be optimized
    Teuchos::RCP<MatParManager> matman_;

    //! the inverse analysis algorithm object
    StatInvAnalysis& ia_;



};
} //namespace INVANA
} //namespace STR


#endif
