/*---------------------------------------------------------------------------*/
/*! \file
\brief particle wall result test for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_wall_result_test.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_particle_wall_datastate.hpp"
#include "4C_particle_wall_interface.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEWALL::WallResultTest::WallResultTest() : CORE::UTILS::ResultTest("PARTICLEWALL")
{
  // empty constructor
}

void PARTICLEWALL::WallResultTest::Init()
{
  // nothing to do
}

void PARTICLEWALL::WallResultTest::Setup(
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface)
{
  // set interface to particle wall handler
  particlewallinterface_ = particlewallinterface;

  // get wall discretization
  walldiscretization_ = particlewallinterface_->GetWallDiscretization();
}

void PARTICLEWALL::WallResultTest::TestNode(INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // extract and check discretization name
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != walldiscretization_->Name()) return;

  // extract node id
  int node;
  res.ExtractInt("NODE", node);
  node -= 1;

  int havenode(walldiscretization_->HaveGlobalNode(node));
  int havenodeonanyproc(0);
  walldiscretization_->Comm().SumAll(&havenode, &havenodeonanyproc, 1);

  // safety check
  if (not havenodeonanyproc)
    FOUR_C_THROW("node %d does not belong to discretization %s", node + 1,
        walldiscretization_->Name().c_str());


  if (walldiscretization_->HaveGlobalNode(node))
  {
    const DRT::Node* actnode = walldiscretization_->gNode(node);

    // node not owned on this processor
    if (actnode->Owner() != walldiscretization_->Comm().MyPID()) return;

    // get wall data state container
    std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
        particlewallinterface_->GetWallDataState();

    // extract test quantity
    std::string quantity;
    res.ExtractString("QUANTITY", quantity);

    // init actual result
    double actresult = 0.0;

    // position
    if (quantity == "posx" or quantity == "posy" or quantity == "posz")
    {
      // get wall displacements
      Teuchos::RCP<const Epetra_Vector> disp = walldatastate->GetDispCol();

      int idx = -1;
      if (quantity == "posx")
        idx = 0;
      else if (quantity == "posy")
        idx = 1;
      else if (quantity == "posz")
        idx = 2;

      if (idx >= 0)
      {
        actresult = actnode->X()[idx];

        if (disp != Teuchos::null)
        {
          const Epetra_BlockMap& disnpmap = disp->Map();
          int lid = disnpmap.LID(walldiscretization_->Dof(0, actnode, idx));
          if (lid < 0)
            FOUR_C_THROW("You tried to test %s on nonexistent dof %d on node %d", quantity.c_str(),
                idx, actnode->Id());
          actresult += (*disp)[lid];
        }
      }
    }
    // displacement
    else if (quantity == "dispx" or quantity == "dispy" or quantity == "dispz")
    {
      // get wall displacements
      Teuchos::RCP<const Epetra_Vector> disp = walldatastate->GetDispCol();

      if (disp == Teuchos::null) return;

      int idx = -1;
      if (quantity == "dispx")
        idx = 0;
      else if (quantity == "dispy")
        idx = 1;
      else if (quantity == "dispz")
        idx = 2;

      if (idx >= 0)
      {
        const Epetra_BlockMap& disnpmap = disp->Map();
        int lid = disnpmap.LID(walldiscretization_->Dof(0, actnode, idx));
        if (lid < 0)
          FOUR_C_THROW("You tried to test %s on nonexistent dof %d on node %d", quantity.c_str(),
              idx, actnode->Id());
        actresult = (*disp)[lid];
      }
    }
    else
      FOUR_C_THROW("result check failed with unknown quantity '%s'!", quantity.c_str());

    // compare values
    const int err = CompareValues(actresult, "NODE", res);
    nerr += err;
    test_count++;
  }
}

void PARTICLEWALL::WallResultTest::TestSpecial(
    INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // check results only for processor 0
  if (walldiscretization_->Comm().MyPID() != 0) return;

  // extract and check discretization name
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != walldiscretization_->Name()) return;

  // extract test quantity
  std::string quantity;
  res.ExtractString("QUANTITY", quantity);

  // init actual result
  double actresult = 0.0;

  // number of total wall elements
  if (quantity == "nwalleles") actresult = walldiscretization_->NumGlobalElements();
  // number of total wall nodes
  else if (quantity == "nwallnodes")
    actresult = walldiscretization_->NumGlobalNodes();
  else
    FOUR_C_THROW("result check failed with unknown quantity '%s'!", quantity.c_str());

  // compare values
  const int err = CompareValues(actresult, "SPECIAL", res);
  nerr += err;
  test_count++;
}

FOUR_C_NAMESPACE_CLOSE
