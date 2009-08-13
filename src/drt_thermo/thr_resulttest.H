/*----------------------------------------------------------------------*/
/*!
\file thr_resulttest.H

\brief testing of structure calculation results

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef THR_RESULTTEST_H
#define THR_RESULTTEST_H

#ifdef CCADISCRET

#include "thrtimint.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"


namespace THR
{

/*!
  \brief Thermo specific result test class

  \author cd
*/
class ResultTest : public DRT::ResultTest
{
public:

  //! Constructor for time integrators of general kind
  //! \author bborn \date 06/08
  ResultTest(TimInt& tintegrator);

  //! \brief thermo version of nodal value tests
  //!
  //! Possible position flags are "temp",
  //!                             "rate",
  virtual void TestNode(const _RESULTDESCR* res, int& nerr, int& test_count);

  //! return true if this is a thermo test
  virtual bool Match(const _RESULTDESCR*);

  /* _RESULTDESCR has been found working, but its exact function is unknown due to
   * missing documentation
   * REMARK: _RESULTDESCR is documented in headers/standardtypes.h ,  but should be changed to OO
   * and moved to a propper place :)
   */

private:

  //! our discretisation
  Teuchos::RCP<DRT::Discretization> thrdisc_;
  // our solution
  //! global temperature DOFs
  Teuchos::RCP<Epetra_Vector> temp_;
  //! global temperature rate DOFs
  Teuchos::RCP<Epetra_Vector> rate_;
  /* NOTE: these have to be present explicitly
   * as they are not part of the problem instance like in fluid3
   */
};

}  // namespace THR

#endif // #ifdef CCADISCRET

#endif // #ifndef THR_RESULTTEST_H
