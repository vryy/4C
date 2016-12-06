/*----------------------------------------------------------------------*/
/*!
\file combust_fluidresulttest.H

\brief tesing of fluid calculation results

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST_FLUIDRESULTTEST_H
#define COMBUST_FLUIDRESULTTEST_H


#include "../drt_lib/drt_resulttest.H"


class Epetra_Vector;
namespace DRT
{
  class Discretization;
  class DofSet;
}

namespace FLD
{
  class CombustFluidImplicitTimeInt;

/*!
  \brief CombustFluid specific result test class

  \author u.kue
*/
class CombustFluidResultTest : public DRT::ResultTest
{
public:

  /*!
  \brief not documented yet
  */
  explicit CombustFluidResultTest(CombustFluidImplicitTimeInt& fluid);

  /// our version of nodal value tests
  /*!
    Possible position flags are "velx", "vely", "velz" and
    "pressure". With the obvious meaning.
   */
  virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

private:

  /// pointer to fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis_;
  /// pointer to standard (no XFEM) fluid dofset
  Teuchos::RCP<const DRT::DofSet>   fluidstddofset_;
  /// pointer to unknown vector with values
  Teuchos::RCP<Epetra_Vector>       mysol_;
};

}

#endif
