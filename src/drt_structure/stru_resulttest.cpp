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
StruResultTest::StruResultTest(RefCountPtr<DRT::Discretization> strudis_in,
                               RefCountPtr<Epetra_Vector> dis,
                               RefCountPtr<Epetra_Vector> vel,
                               RefCountPtr<Epetra_Vector> acc)
{
  strudisc_ = strudis_in;
  dis_ = dis;
  vel_ = vel;
  acc_ = acc;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StruResultTest::TestNode(_RESULTDESCR* res, int& nerr, int& test_count)
{
  if (res->dis != 0)
    dserror("fix me: only one structure discretization supported for testing");

  /* this implementation does not allow testing of stresses
   */
  if (strudisc_->HaveGlobalNode(res->node))
  {
    DRT::Node* actnode = strudisc_->gNode(res->node);

    // Strange! It seems we might actually have a global node around
    // even if it does not belong to us. But here we are just
    // interested in our nodes!
    if (actnode->Owner() != strudisc_->Comm().MyPID())
      return;

    // verbose output
    //cout << "TESTING STRUCTURE RESULTS with StruResultTest::TestNode(..)" << endl;

    double result = 0;  // will hold the actual result of run
    string position = res->position;  // type of result value
    bool unknownpos = true;  // make sure the result value string can be handled

    // test displacements
    if (dis_ != null)
    {
      const Epetra_BlockMap& disnpmap = dis_->Map();

      if (position=="dispx")
      {
        unknownpos = false;
        result = (*dis_)[disnpmap.LID(strudisc_->Dof(actnode,0))];
      }
      else if (position=="dispy")
      {
        unknownpos = false;
        result = (*dis_)[disnpmap.LID(strudisc_->Dof(actnode,1))];
      }
      else if (position=="dispz")
      {
        unknownpos = false;
        result = (*dis_)[disnpmap.LID(strudisc_->Dof(actnode,2))];
      }
    }

    // test velocities
    if (vel_ != null)
    {
      const Epetra_BlockMap& velnpmap = vel_->Map();

      if (position=="velx")
      {
        unknownpos = false;
        result = (*vel_)[velnpmap.LID(strudisc_->Dof(actnode,0))];
      }
      else if (position=="vely")
      {
        unknownpos = false;
        result = (*vel_)[velnpmap.LID(strudisc_->Dof(actnode,1))];
      }
      else if (position=="velz")
      {
        unknownpos = false;
        result = (*vel_)[velnpmap.LID(strudisc_->Dof(actnode,2))];
      }
    }

    // test accelerations
    if (acc_ != null)
    {
      const Epetra_BlockMap& accnpmap = acc_->Map();

      if (position=="accx")
      {
        unknownpos = false;
        result = (*acc_)[accnpmap.LID(strudisc_->Dof(actnode,0))];
      }
      else if (position=="accy")
      {
        unknownpos = false;
        result = (*acc_)[accnpmap.LID(strudisc_->Dof(actnode,1))];
      }
      else if (position=="accz")
      {
        unknownpos = false;
        result = (*acc_)[accnpmap.LID(strudisc_->Dof(actnode,2))];
      }
    }

    // catch position strings, which are not handled by structure result test
    if (unknownpos)
      dserror("position '%s' not supported in structure testing", position.c_str());

    // verbose output
    cout.precision(18);
    cout << "RESULT IS " << std::scientific << result << endl;

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
