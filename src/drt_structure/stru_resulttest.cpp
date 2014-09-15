/*----------------------------------------------------------------------*/
/*!
\file stru_resulttest.cpp

\brief tesing of structure calculation results

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#include <string>
#include "stru_resulttest.H"
#include "strtimint.H"
#include "../drt_particle/particle_timint.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
StruResultTest::StruResultTest(STR::TimInt& tintegrator)
  : DRT::ResultTest("STRUCTURE")
{
  dis_  = tintegrator.Dis();
  vel_  = tintegrator.Vel();
  acc_  = tintegrator.Acc();
  strudisc_ = tintegrator.Discretization();
  if (DRT::Problem::Instance()->ProblemType() == prb_struct_ale and
      (DRT::Problem::Instance()->WearParams()).get<double>("WEARCOEFF")>0.0)
    dism_ = tintegrator.Dismat();
  else
    dism_=Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
StruResultTest::StruResultTest(PARTICLE::TimInt& tintegrator)
  : DRT::ResultTest("PARTICLE")
{
  dis_  = tintegrator.Dispnp();
  if (DRT::Problem::Instance()->ProblemType() != prb_level_set and
      DRT::Problem::Instance()->ProblemType() != prb_combust and
      DRT::Problem::Instance()->ProblemType() != prb_two_phase_flow)
  {
    vel_  = tintegrator.Veln();
    acc_  = tintegrator.Accn();
  }
  strudisc_ = tintegrator.Discretization();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StruResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
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
      if (dis_ != Teuchos::null)
      {
        const Epetra_BlockMap& disnpmap = dis_->Map();
        if (position=="dispx")
        {
          unknownpos = false;
          result = (*dis_)[disnpmap.LID(strudisc_->Dof(0,actnode,0))];
        }
        else if (position=="dispy")
        {
          unknownpos = false;
          result = (*dis_)[disnpmap.LID(strudisc_->Dof(0,actnode,1))];
        }
        else if (position=="dispz")
        {
          unknownpos = false;
          result = (*dis_)[disnpmap.LID(strudisc_->Dof(0,actnode,2))];
        }
        else if (position=="press")
        {
          unknownpos = false;
          result = (*dis_)[disnpmap.LID(strudisc_->Dof(0,actnode,3))];
        }
      }

      // test material displacements
      if (dism_ != Teuchos::null)
      {
        const Epetra_BlockMap& dismpmap = dism_->Map();
        if (position=="dispmx")
        {
          unknownpos = false;
          result = (*dism_)[dismpmap.LID(strudisc_->Dof(0,actnode,0))];
        }
        else if (position=="dispmy")
        {
          unknownpos = false;
          result = (*dism_)[dismpmap.LID(strudisc_->Dof(0,actnode,1))];
        }
        else if (position=="dispmz")
        {
          unknownpos = false;
          result = (*dism_)[dismpmap.LID(strudisc_->Dof(0,actnode,2))];
        }
      }

      // test velocities
      if (vel_ != Teuchos::null)
      {
        const Epetra_BlockMap& velnpmap = vel_->Map();

        if (position=="velx")
        {
          unknownpos = false;
          result = (*vel_)[velnpmap.LID(strudisc_->Dof(0,actnode,0))];
        }
        else if (position=="vely")
        {
          unknownpos = false;
          result = (*vel_)[velnpmap.LID(strudisc_->Dof(0,actnode,1))];
        }
        else if (position=="velz")
        {
          unknownpos = false;
          result = (*vel_)[velnpmap.LID(strudisc_->Dof(0,actnode,2))];
        }
      }

      // test accelerations
      if (acc_ != Teuchos::null)
      {
        const Epetra_BlockMap& accnpmap = acc_->Map();

        if (position=="accx")
        {
          unknownpos = false;
          result = (*acc_)[accnpmap.LID(strudisc_->Dof(0,actnode,0))];
        }
        else if (position=="accy")
        {
          unknownpos = false;
          result = (*acc_)[accnpmap.LID(strudisc_->Dof(0,actnode,1))];
        }
        else if (position=="accz")
        {
          unknownpos = false;
          result = (*acc_)[accnpmap.LID(strudisc_->Dof(0,actnode,2))];
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
