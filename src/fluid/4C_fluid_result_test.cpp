/*-----------------------------------------------------------*/
/*! \file

\brief testing of fluid calculation results


\level 1

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_result_test.hpp"

#include "4C_fluid_implicit_integration.hpp"
#include "4C_fluid_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"

#include <string>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::FluidResultTest::FluidResultTest(FluidImplicitTimeInt& fluid)
    : Core::UTILS::ResultTest("FLUID")
{
  fluiddis_ = fluid.discret_;
  myerror_ = fluid.evaluate_error_compared_to_analytical_sol();
  mysol_ = fluid.velnp_;

  // quantities not implemented in the HDG formulation
  if (Global::Problem::instance()->spatial_approximation_type() != Core::FE::ShapeFunctionType::hdg)
  {
    mytraction_ = fluid.stressmanager_->get_pre_calc_stresses(fluid.trueresidual_);
    mywss_ = fluid.stressmanager_->get_pre_calc_wall_shear_stresses(fluid.trueresidual_);
    mydivu_ = fluid.evaluate_div_u();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::FluidResultTest::test_node(Input::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.extract_string("DIS", dis);
  if (dis != fluiddis_->name()) return;

  int node;
  res.extract_int("NODE", node);
  node -= 1;

  int havenode(fluiddis_->have_global_node(node));
  int isnodeofanybody(0);
  fluiddis_->get_comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    FOUR_C_THROW(
        "Node %d does not belong to discretization %s", node + 1, fluiddis_->name().c_str());
  }
  else
  {
    if (fluiddis_->have_global_node(node))
    {
      const Core::Nodes::Node* actnode = fluiddis_->g_node(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->owner() != fluiddis_->get_comm().MyPID()) return;

      double result = 0.;

      const Epetra_BlockMap& velnpmap = mysol_->Map();

      const int numdim = Global::Problem::instance()->n_dim();

      std::string position;
      res.extract_string("QUANTITY", position);
      if (position == "velx")
        result = (*mysol_)[velnpmap.LID(fluiddis_->dof(0, actnode, 0))];
      else if (position == "vely")
        result = (*mysol_)[velnpmap.LID(fluiddis_->dof(0, actnode, 1))];
      else if (position == "velz")
      {
        if (numdim == 2) FOUR_C_THROW("Cannot test result for velz in 2D case.");
        result = (*mysol_)[velnpmap.LID(fluiddis_->dof(0, actnode, 2))];
      }
      else if (position == "pressure")
      {
        if (fluiddis_->num_dof(0, actnode) < (numdim + 1))
          FOUR_C_THROW("too few dofs at node %d for pressure testing", actnode->id());
        result = (*mysol_)[velnpmap.LID(fluiddis_->dof(0, actnode, numdim))];
      }
      else if (position == "tractionx")
        result = (*mytraction_)[(mytraction_->Map()).LID(fluiddis_->dof(0, actnode, 0))];
      else if (position == "tractiony")
        result = (*mytraction_)[(mytraction_->Map()).LID(fluiddis_->dof(0, actnode, 1))];
      else if (position == "tractionz")
      {
        if (numdim == 2) FOUR_C_THROW("Cannot test result for tractionz in 2D case.");
        result = (*mytraction_)[(mytraction_->Map()).LID(fluiddis_->dof(0, actnode, 2))];
      }
      else if (position == "wssx")
        result = (*mywss_)[(mywss_->Map()).LID(fluiddis_->dof(0, actnode, 0))];
      else if (position == "wssy")
        result = (*mywss_)[(mywss_->Map()).LID(fluiddis_->dof(0, actnode, 1))];
      else if (position == "wssz")
      {
        if (numdim == 2) FOUR_C_THROW("Cannot test result for wssz in 2D case.");
        result = (*mywss_)[(mywss_->Map()).LID(fluiddis_->dof(0, actnode, 2))];
      }
      else if (position == "L2errvel")
        result = (*myerror_)[0];
      else if (position == "L2errpre")
        result = (*myerror_)[1];
      else if (position == "H1errvel")
        result = (*myerror_)[2];
      else if (position == "divu")
        result = (*mydivu_);
      else if (position == "L2errmix")
        result = (*myerror_)[0];
      else if (position == "L2errden")
        result = (*myerror_)[1];
      else if (position == "L2errmom")
        result = (*myerror_)[2];
      else
        FOUR_C_THROW("Quantity '%s' not supported in fluid testing", position.c_str());

      nerr += compare_values(result, "NODE", res);
      test_count++;
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
