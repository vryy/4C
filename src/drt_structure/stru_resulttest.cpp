/*----------------------------------------------------------------------*/
/*!
\file stru_resulttest.cpp

\brief tesing of structure calculation results

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <string>

#include "stru_resulttest.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
//#include "../headers/standardtypes.h"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
StruResultTest::StruResultTest(
	RefCountPtr<DRT::Discretization> strudis_in,
	RefCountPtr<Epetra_Vector> dis_)
{
  stru_dis_ =strudis_in;
  mysol_   =dis_;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StruResultTest::TestNode(_RESULTDESCR* res, int& nerr, int& test_count)
{
  if (res->dis != 0)
    dserror("fix me: only one structure discretization supported for testing");

  /* this implementation does not allow testing of stresses
   */
  if (stru_dis_->HaveGlobalNode(res->node))
  {
    DRT::Node* actnode = stru_dis_->gNode(res->node);

    // Strange! It seems we might actually have a global node around
    // even if it does not belong to us. But here we are just
    // interested in our nodes!
    if (actnode->Owner() != stru_dis_->Comm().MyPID())
      return;

    double result = 0.;

    const Epetra_BlockMap& velnpmap = mysol_->Map();

    string position = res->position;

    //verbose output
    cout << "TESTING STRUCTURE RESULTS with StruResultTest::TestNode(..)" << endl;

    if (position=="dispx")
    {
      result = (*mysol_)[velnpmap.LID(stru_dis_->Dof(actnode,0))];
    }
  	else if (position=="dispy")
    {
      result = (*mysol_)[velnpmap.LID(stru_dis_->Dof(actnode,1))];
    }
    else if (position=="dispz")
    {
      result = (*mysol_)[velnpmap.LID(stru_dis_->Dof(actnode,2))];
    }
    else
    {
      dserror("position '%s' not supported in structure testing", position.c_str());
    }

	//verbose output
    cout.precision(18);
	  cout << "RESULT IS " << result << endl;

    nerr += CompareValues(result, res);
    test_count++;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool StruResultTest::Match(_RESULTDESCR* res)
{
  /* res.field is a enum of type _FIELDTYP and can be found in headers/enums.h
   */
  return res->field==structure;
}


#endif /* CCADISCRET       */
