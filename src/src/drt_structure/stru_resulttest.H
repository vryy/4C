/*----------------------------------------------------------------------*/
/*!
\file stru_resulttest.H

\brief testing of structure calculation results

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef STRURESULTTEST_H
#define STRURESULTTEST_H

#ifdef CCADISCRET

#include "stru_static_drt.H"
#include "strugenalpha.H"
#include "strtimint.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"


/*!
  \brief Structure specific result test class

  \author a.vlf
*/
class StruResultTest : public DRT::ResultTest
{
public:
  //! constructor
  StruResultTest(
    RCP<DRT::Discretization> strudis_in,  //!< structure discretisation
    RCP<Epetra_Vector> dis,  //!< displacements to test
    RCP<Epetra_Vector> vel,  //!< velocities to test
    RCP<Epetra_Vector> acc  //!< accelerations to test
  );

  //! constructor to use when timeintegrator is generalized alpha!
  StruResultTest(StruGenAlpha& tintegrator);

  //! Constructor for time integrators of general kind
  //! \author bborn \date 06/08
  StruResultTest(STR::TimInt& tintegrator);

  //! \brief structure version of nodal value tests
  //!
  //! Possible position flags are "dispx", "dispy", "dispz",
  //!                             "velx", "vely", "velz",
  //!                             "accx", "accy", "accz"
  virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

  //! return true if this is a structure test
  virtual bool Match(DRT::INPUT::LineDefinition& res);

private:

  //! our discretisation
  RCP<DRT::Discretization> strudisc_;
  // our solution
  //! global displacement DOFs
  RCP<Epetra_Vector> dis_;
  //! global velocity DOFs
  RCP<Epetra_Vector> vel_;
  //! global acceleration DOFs
  RCP<Epetra_Vector> acc_;
  /* NOTE: these have to be present explicitly
   * as they are not part of the problem instance like in fluid3
   */
};

#endif /* #ifdef CCADISCRET        */

#endif /* #ifndef STRURESULTTEST_H */
