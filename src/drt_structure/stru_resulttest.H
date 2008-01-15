/*----------------------------------------------------------------------*/
/*!
\file stru_resulttest.H

\brief tesing of structure calculation results

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
#ifdef RESULTTEST

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
  /*!
  \brief constructor 
  */
  StruResultTest(RefCountPtr<DRT::Discretization> strudis_in, RefCountPtr<Epetra_Vector> dis_);

  // structure version of nodal value tests
  /*!
   Possible position flags are "dispx", "dispy", "dispz"
  */
  virtual void TestNode(_RESULTDESCR* res, int& nerr, int& test_count);

  // return true if this is a structure test
  virtual bool Match(_RESULTDESCR*);
  
  /* _RESULTDESCR has been found working, but its exact function is unknown due to 
   * missing documentation
   * REMARK: _RESULTDESCR is documented in headers/standardtypes.h ,  but should be changed to OO
   * and moved to a propper place :)
   */
   
private:

  // our discretisation 
  RefCountPtr<DRT::Discretization> stru_dis_;
  // our solution
  RefCountPtr<Epetra_Vector> mysol_;
  /* NOTE: these have to be present explicitly 
   * as they are not part of the problem instance like in fluid3
   */
};

#endif /* #ifdef RESULTTEST        */
#endif /* #ifdef CCADISCRET        */

#endif /* #ifndef STRURESULTTEST_H */
