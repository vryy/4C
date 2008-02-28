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
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../io/io_drt.H"


/*!
  \class StruResultTest

  \brief Structure specific result test class

  \author a.vlf
*/
class StruResultTest : public DRT::ResultTest
{
public:
  //! constructor
  StruResultTest(
    RefCountPtr<DRT::Discretization> strudis_in,  //!< structure discretisation
    RefCountPtr<Epetra_Vector> dis,  //!< displacements to test
    RefCountPtr<Epetra_Vector> vel,  //!< velocities to test
    RefCountPtr<Epetra_Vector> acc  //!< accelerations to test
  );

  //! \brief structure version of nodal value tests
  //!
  //! Possible position flags are "dispx", "dispy", "dispz", 
  //!                             "velx", "vely", "velz",
  //!                             "accx", "accy", "accz"
  virtual void TestNode(_RESULTDESCR* res, int& nerr, int& test_count);

  //! return true if this is a structure test
  virtual bool Match(_RESULTDESCR*);

  /* _RESULTDESCR has been found working, but its exact function is unknown due to
   * missing documentation
   * REMARK: _RESULTDESCR is documented in headers/standardtypes.h ,  but should be changed to OO
   * and moved to a propper place :)
   */

private:

  //! our discretisation
  RefCountPtr<DRT::Discretization> strudisc_;
  // our solution
  //! global displacement DOFs
  RefCountPtr<Epetra_Vector> dis_;
  //! global velocity DOFs
  RefCountPtr<Epetra_Vector> vel_;
  //! global acceleration DOFs
  RefCountPtr<Epetra_Vector> acc_;
  /* NOTE: these have to be present explicitly
   * as they are not part of the problem instance like in fluid3
   */
};

#endif /* #ifdef CCADISCRET        */

#endif /* #ifndef STRURESULTTEST_H */
