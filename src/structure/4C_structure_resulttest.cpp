/*----------------------------------------------------------------------*/
/*! \file
\brief Testing of structure calculation results


\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_structure_resulttest.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_structure_timint.hpp"

#include <Epetra_Vector.h>

#include <string>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
StruResultTest::StruResultTest(STR::TimInt& tintegrator)
    : CORE::UTILS::ResultTest("STRUCTURE"), timeintegrator_(Teuchos::rcpFromRef(tintegrator))
{
  dis_ = tintegrator.Dis();
  vel_ = tintegrator.Vel();
  acc_ = tintegrator.Acc();
  strudisc_ = tintegrator.discretization();

  if (tintegrator.DispMat() != Teuchos::null)
    dism_ = tintegrator.Dismat();
  else
    dism_ = Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StruResultTest::test_node(INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // this implementation does not allow testing of stresses !

  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != strudisc_->Name()) return;

  int node;
  res.ExtractInt("NODE", node);
  node -= 1;

  int havenode(strudisc_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  strudisc_->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    FOUR_C_THROW(
        "Node %d does not belong to discretization %s", node + 1, strudisc_->Name().c_str());
  }
  else
  {
    if (strudisc_->HaveGlobalNode(node))
    {
      const CORE::Nodes::Node* actnode = strudisc_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != strudisc_->Comm().MyPID()) return;

      std::string position;
      res.ExtractString("QUANTITY", position);
      bool unknownpos = true;  // make sure the result value std::string can be handled
      double result = 0.0;     // will hold the actual result of run

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
          int lid = disnpmap.LID(strudisc_->Dof(0, actnode, idx));
          if (lid < 0)
            FOUR_C_THROW("You tried to test %s on nonexistent dof %d on node %d", position.c_str(),
                idx, actnode->Id());
          result = (*dis_)[lid];
        }
      }

      // test material displacements
      if (dism_ != Teuchos::null)
      {
        const Epetra_BlockMap& dismpmap = dism_->Map();
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
          int lid = dismpmap.LID(strudisc_->Dof(0, actnode, idx));
          if (lid < 0)
            FOUR_C_THROW("You tried to test %s on nonexistent dof %d on node %d", position.c_str(),
                idx, actnode->Id());
          result = (*dism_)[lid];
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
          int lid = velnpmap.LID(strudisc_->Dof(0, actnode, idx));
          if (lid < 0)
            FOUR_C_THROW("You tried to test %s on nonexistent dof %d on node %d", position.c_str(),
                idx, actnode->Id());
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
          int lid = accnpmap.LID(strudisc_->Dof(0, actnode, idx));
          if (lid < 0)
            FOUR_C_THROW("You tried to test %s on nonexistent dof %d on node %d", position.c_str(),
                idx, actnode->Id());
          result = (*acc_)[lid];
        }
      }

      // catch position std::strings, which are not handled by structure result test
      if (unknownpos)
        FOUR_C_THROW("Quantity '%s' not supported in structure testing", position.c_str());

      // compare values
      const int err = compare_values(result, "NODE", res);
      nerr += err;
      test_count++;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StruResultTest::TestSpecial(INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // extract name of quantity to be tested
  std::string quantity;
  res.ExtractString("QUANTITY", quantity);

  // get result to be tested on all processors
  const double result = get_special_result_for_testing(quantity);

  // compare values on one processor only, as they are the same everywhere
  if (strudisc_->Comm().MyPID() == 0)
  {
    const int err = compare_values(result, "SPECIAL", res);
    nerr += err;
    test_count++;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double StruResultTest::get_special_result_for_testing(const std::string& quantity)
{
  // initialize variable for result
  double result(0.);

  if (quantity == "lin_iters")  // number of iterations in solid linear solver
    result = static_cast<double>(timeintegrator_->Solver()->getNumIters());
  else if (quantity == "lin_iters_contact")  // number of iterations in contact linear solver
    result = static_cast<double>(timeintegrator_->ContactSolver()->getNumIters());
  else  // Catch unknown quantity strings
    FOUR_C_THROW("Quantity '%s' not supported in structure result test!", quantity.c_str());

  return result;
}

FOUR_C_NAMESPACE_CLOSE
