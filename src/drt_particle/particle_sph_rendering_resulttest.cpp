/*----------------------------------------------------------------------*/
/*!
\file particle_sph_rendering_resulttest.cpp

\brief testing of SPH rendering calculation results

\level 2

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                 sfuchs 06/17 |
 *----------------------------------------------------------------------*/
#include "particle_sph_rendering_resulttest.H"

#include <string>

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "particle_sph_rendering.H"

/*----------------------------------------------------------------------*
 | constructor                                             sfuchs 06/17 |
 *----------------------------------------------------------------------*/
ParticleSPHRenderingResultTest::ParticleSPHRenderingResultTest(PARTICLE::Rendering& rendering)
  : DRT::ResultTest("PARTICLE_RENDERING")
{
  discret_ = rendering.GetRenderingDiscret();

  vel_ = rendering.GetRenderingVelocity();
  acc_ = rendering.GetRenderingAcceleration();
  velmod_ = rendering.GetRenderingVelocityMod();
  density_ = rendering.GetRenderingDensity();
  pressure_ = rendering.GetRenderingPressure();

  return;
} // ParticleSPHRenderingResultTest::ParticleSPHRenderingResultTest

/*----------------------------------------------------------------------*
 |                                                         sfuchs 06/17 |
 *----------------------------------------------------------------------*/
void ParticleSPHRenderingResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS",dis);
  if (dis != discret_->Name())
    return;

  int node;
  res.ExtractInt("NODE",node);

  int havenode(discret_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  discret_->Comm().SumAll(&havenode,&isnodeofanybody,1);

  if (isnodeofanybody==0)
  {
    dserror("Node %d does not belong to discretization %s",node,discret_->Name().c_str());
  }
  else
  {
    if (discret_->HaveGlobalNode(node))
    {
      const DRT::Node* actnode = discret_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != discret_->Comm().MyPID())
        return;

      std::string position;
      res.ExtractString("QUANTITY",position);
      bool unknownpos = true;  // make sure the result value std::string can be handled
      double result = 0.0;  // will hold the actual result of run

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
          int lid = velnpmap.LID(discret_->Dof(0,actnode,idx));
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
          int lid = accnpmap.LID(discret_->Dof(0,actnode,idx));
          if (lid < 0)
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx, actnode->Id());
          result = (*acc_)[lid];
        }
      }

      // test modified velocities
      if (velmod_ != Teuchos::null)
      {
        const Epetra_BlockMap& velmodnpmap = velmod_->Map();
        int idx = -1;
        if (position == "velmodx")
          idx = 0;
        else if (position == "velmody")
          idx = 1;
        else if (position == "velmodz")
          idx = 2;

        if (idx >= 0)
        {
          unknownpos = false;
          int lid = velmodnpmap.LID(discret_->Dof(0,actnode,idx));
          if (lid < 0)
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx, actnode->Id());
          result = (*velmod_)[lid];
        }
      }

      // test density
      if (density_ != Teuchos::null)
      {
        const Epetra_BlockMap& densitymap = density_->Map();
        int idx = -1;
        if (position == "density")
          idx = 0;

        if (idx >= 0)
        {
          unknownpos = false;
          // node based vector
          int lid = densitymap.LID(actnode->Id());
          if (lid < 0)
            dserror("You tried to test %s on nonexistent node %d", position.c_str(), actnode->Id());
          result = (*density_)[lid];
        }
      }

      // test pressure
      if (pressure_ != Teuchos::null)
      {
        const Epetra_BlockMap& pressuremap = pressure_->Map();
        int idx = -1;
        if (position == "pressure")
          idx = 0;

        if (idx >= 0)
        {
          unknownpos = false;
          // node based vector
          int lid = pressuremap.LID(actnode->Id());
          if (lid < 0)
            dserror("You tried to test %s on nonexistent node %d", position.c_str(), actnode->Id());
          result = (*pressure_)[lid];
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

  return;
} // ParticleSPHRenderingResultTest::TestNode
