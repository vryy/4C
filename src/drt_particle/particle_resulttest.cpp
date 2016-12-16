/*----------------------------------------------------------------------*/
/*!
\file particle_resulttest.cpp

\brief testing of particle calculation results

\level 2

<pre>
\maintainer Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/

#include <string>
#include "particle_resulttest.H"
#include "../drt_particle/particle_timint.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PartResultTest::PartResultTest(PARTICLE::TimInt& tintegrator)
  : DRT::ResultTest("PARTICLE")
{
  partdisc_ = tintegrator.Discretization();

    dis_  = tintegrator.Dispnp();
  if (tintegrator.Velnp() != Teuchos::null)
    vel_  = tintegrator.Velnp();
  if (tintegrator.Accnp() != Teuchos::null)
    acc_  = tintegrator.Accnp();

  if (tintegrator.Radiusnp() != Teuchos::null)
    radius_ = tintegrator.Radiusnp();
  else if (tintegrator.Radiusn() != Teuchos::null)
    radius_ = tintegrator.Radiusn();

  if (tintegrator.Densitynp() != Teuchos::null)
    density_  = tintegrator.Densitynp();
  if (tintegrator.Temperature() != Teuchos::null)
    temperature_ = tintegrator.Temperature();
  if (tintegrator.Temperature() != Teuchos::null)
    pressure_ = tintegrator.Pressure();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PartResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS",dis);
  if (dis != partdisc_->Name())
    return;

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  int havenode(partdisc_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  partdisc_->Comm().SumAll(&havenode,&isnodeofanybody,1);

  if (isnodeofanybody==0)
  {
    dserror("Node %d does not belong to discretization %s",node+1,partdisc_->Name().c_str());
  }
  else
  {
    if (partdisc_->HaveGlobalNode(node))
    {
      const DRT::Node* actnode = partdisc_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != partdisc_->Comm().MyPID())
        return;

      std::string position;
      res.ExtractString("QUANTITY",position);
      bool unknownpos = true;  // make sure the result value std::string can be handled
      double result = 0.0;  // will hold the actual result of run

      // test displacements or pressure
      if (dis_ != Teuchos::null)
      {
        const Epetra_BlockMap& disnpmap = dis_->Map();
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
          int lid = disnpmap.LID(partdisc_->Dof(0,actnode,idx));
          if (lid < 0)
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx, actnode->Id());
          result = (*dis_)[lid];
        }
      }

      // test velocities
      if (vel_ != Teuchos::null)
      {
        const Epetra_BlockMap& velnpmap = vel_->Map();
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
          int lid = velnpmap.LID(partdisc_->Dof(0,actnode,idx));
          if (lid < 0)
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx, actnode->Id());
          result = (*vel_)[lid];
        }
      }

      // test accelerations
      if (acc_ != Teuchos::null)
      {
        const Epetra_BlockMap& accnpmap = acc_->Map();
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
          int lid = accnpmap.LID(partdisc_->Dof(0,actnode,idx));
          if (lid < 0)
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx, actnode->Id());
          result = (*acc_)[lid];
        }
      }

      // test radius
      if (radius_ != Teuchos::null)
      {
        const Epetra_BlockMap& radiusmap = radius_->Map();
        int idx = -1;
        if (position == "radius")
          idx = 0;

        if (idx >= 0)
        {
          unknownpos = false;
          // node based vector
          int lid = radiusmap.LID(actnode->Id());
          if (lid < 0)
            dserror("You tried to test %s on nonexistent node %d", position.c_str(), actnode->Id());
          result = (*radius_)[lid];
        }
      }

      // test temperature
      if (temperature_ != Teuchos::null)
      {
        const Epetra_BlockMap& temperaturenpmap = temperature_->Map();
        int idx = -1;
        if (position == "temperature")
          idx = 0;

        if (idx >= 0)
        {
          unknownpos = false;
          // node based vector
          int lid = temperaturenpmap.LID(actnode->Id());
          if (lid < 0)
            dserror("You tried to test %s on nonexistent node %d", position.c_str(), actnode->Id());
          result = (*temperature_)[lid];
        }
      }

      // test density
      if (density_ != Teuchos::null)
      {
        const Epetra_BlockMap& densitynpmap = density_->Map();
        int idx = -1;
        if (position == "density")
          idx = 0;

        if (idx >= 0)
        {
          unknownpos = false;
          // node based vector
          int lid = densitynpmap.LID(actnode->Id());
          if (lid < 0)
            dserror("You tried to test %s on nonexistent node %d", position.c_str(), actnode->Id());
          result = (*density_)[lid];
        }
      }

      // catch position std::strings, which are not handled by particle result test
      if (unknownpos)
        dserror("Quantity '%s' not supported in particle testing", position.c_str());

      // compare values
      const int err = CompareValues(result, "NODE", res);
      nerr += err;
      test_count++;
    }
  }
}
