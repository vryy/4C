/*-----------------------------------------------------------*/
/*!
\file str_resulttest.cpp

\maintainer Michael Hiermeier

\date Dec 21, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_resulttest.H"
#include "str_timint_base.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::ResultTest::ResultTest()
    : DRT::ResultTest("STRUCTURE"),
      isinit_(false),
      issetup_(false),
      strudisc_(Teuchos::null),
      disn_(Teuchos::null),
      dismatn_(Teuchos::null),
      veln_(Teuchos::null),
      accn_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ResultTest::Init(const STR::TIMINT::BaseDataGlobalState& gstate)
{
  issetup_ = false;

  disn_  = gstate.GetDisN();
  veln_  = gstate.GetVelN();
  accn_  = gstate.GetAccN();
  strudisc_ = gstate.GetDiscret();
  if (DRT::Problem::Instance()->ProblemType() == prb_struct_ale and
      (DRT::Problem::Instance()->WearParams()).get<double>("WEARCOEFF")>0.0)
    dserror("Material displ. are not yet considered!");
  else
    dismatn_=Teuchos::null;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ResultTest::Setup()
{
  CheckInit();
  // currently unused
  issetup_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STR::ResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  CheckInitSetup();
  // this implementation does not allow testing of stresses !

  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS",dis);
  if (dis != strudisc_->Name())
    return;

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  int havenode(strudisc_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  strudisc_->Comm().SumAll(&havenode,&isnodeofanybody,1);

  if (isnodeofanybody==0)
  {
    dserror("Node %d does not belong to discretization %s",node+1,strudisc_->Name().c_str());
  }
  else
  {
    if (strudisc_->HaveGlobalNode(node))
    {
      const DRT::Node* actnode = strudisc_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != strudisc_->Comm().MyPID())
        return;

      std::string position;
      res.ExtractString("QUANTITY",position);
      bool unknownpos = true;  // make sure the result value std::string can be handled
      double result = 0.0;  // will hold the actual result of run

      // test displacements or pressure
      if (disn_ != Teuchos::null)
      {
        const Epetra_BlockMap& disnpmap = disn_->Map();
        int idx = -1;
        if (position == "dispx")
          idx = 0;
        else if (position == "dispy")
          idx = 1;
        else if (position == "dispz")
          idx = 2;
        else if (position == "press")
          idx = 3;

        if (idx >= 0)
        {
          unknownpos = false;
          int lid = disnpmap.LID(strudisc_->Dof(0,actnode,idx));
          if (lid < 0)
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx, actnode->Id());
          result = (*disn_)[lid];
        }
      }

      // test material displacements
      if (!dismatn_.is_null())
      {
        const Epetra_BlockMap& dismpmap = dismatn_->Map();
        int idx = -1;
        if (position == "dispmx")
          idx = 0;
        else if (position == "dispmy")
          idx = 1;
        else if (position == "dispmz")
          idx = 2;

        if (idx >= 0)
        {
          unknownpos = false;
          int lid = dismpmap.LID(strudisc_->Dof(0,actnode,idx));
          if (lid < 0)
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx, actnode->Id());
          result = (*dismatn_)[lid];
        }
      }

      // test velocities
      if (veln_ != Teuchos::null)
      {
        const Epetra_BlockMap& velnpmap = veln_->Map();
        int idx = -1;
        if (position == "velx")
          idx = 0;
        else if (position == "vely")
          idx = 1;
        else if (position == "velz")
          idx = 2;

        if (idx >= 0)
        {
          unknownpos = false;
          int lid = velnpmap.LID(strudisc_->Dof(0,actnode,idx));
          if (lid < 0)
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx, actnode->Id());
          result = (*veln_)[lid];
        }
      }

      // test accelerations
      if (accn_ != Teuchos::null)
      {
        const Epetra_BlockMap& accnpmap = accn_->Map();
        int idx = -1;
        if (position == "accx")
          idx = 0;
        else if (position == "accy")
          idx = 1;
        else if (position == "accz")
          idx = 2;

        if (idx >= 0)
        {
          unknownpos = false;
          int lid = accnpmap.LID(strudisc_->Dof(0,actnode,idx));
          if (lid < 0)
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx, actnode->Id());
          result = (*accn_)[lid];
        }
      }

      // catch position std::strings, which are not handled by structure result test
      if (unknownpos)
        dserror("Quantity '%s' not supported in structure testing", position.c_str());

      // compare values
      const int err = CompareValues(result, "NODE", res);
      nerr += err;
      test_count++;
    }
  }
}
