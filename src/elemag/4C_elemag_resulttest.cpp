/*----------------------------------------------------------------------*/
/*! \file

\brief testing of electromagnetic calculation results

\level 2

*/
/*----------------------------------------------------------------------*/


#include "4C_elemag_resulttest.hpp"

#include "4C_elemag_timeint.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                     berardocco 10/18 |
 *----------------------------------------------------------------------*/
EleMag::ElemagResultTest::ElemagResultTest(ElemagTimeInt& elemagalgo)
    : Core::UTILS::ResultTest("ELECTROMAGNETIC")
{
  dis_ = elemagalgo.discretization();
  // mysol_ = Core::LinAlg::CreateVector(*(dis_->NodeRowMap()), true);
  error_ = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(4));
  error_ = elemagalgo.compute_error();
  // elemagalgo.NodalPressureField(mysol_);
}

/*----------------------------------------------------------------------*
 |                                                     berardocco 10/18 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagResultTest::test_node(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis = container.get<std::string>("DIS");
  if (dis != dis_->name()) return;

  int node = container.get<int>("NODE");
  node -= 1;

  int havenode(dis_->have_global_node(node));
  int isnodeofanybody(0);
  dis_->get_comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    FOUR_C_THROW("Node %d does not belong to discretization %s", node + 1, dis_->name().c_str());
  }
  else
  {
    if (dis_->have_global_node(node))
    {
      Core::Nodes::Node* actnode = dis_->g_node(node);

      // Here, we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->owner() != dis_->get_comm().MyPID()) return;

      double result = 0.;
      // const Epetra_BlockMap& map = mysol_->Map();
      std::string variable = container.get<std::string>("QUANTITY");

      if (variable == "electric")
      {
        std::cout << "This error function has not yet been implemented." << std::endl;
        // result = (*mysol_)[map.LID(actnode->Id())];  // Here, we got a node map, not a dof map
        // with dofs based on nodes!
      }
      else if (variable == "L2electric")
      {
        result = std::sqrt((*error_)[0]);
      }
      else if (variable == "L2magnetic")
      {
        result = std::sqrt((*error_)[2]);
      }
      else if (variable == "L2electric-rel")
      {
        if ((*error_)[1] > 0.0)
          result = std::sqrt((*error_)[0] / (*error_)[1]);
        else
          FOUR_C_THROW(
              "Impossible to compute the electric relative error. The L2-norm of the "
              "analytical solution is zero, resulting in a division by zero.");
      }
      else if (variable == "L2magnetic-rel")
      {
        if ((*error_)[3] > 0.0)
          result = std::sqrt((*error_)[2] / (*error_)[3]);
        else
          FOUR_C_THROW(
              "Impossible to compute the magnetic relative error. The L2-norm of the "
              "analytical solution is zero, resulting in a division by zero.");
      }
      else if (variable == "Hdiv-electric")
      {
        result = std::sqrt((*error_)[4]);
      }
      else if (variable == "Hdiv-magnetic")
      {
        result = std::sqrt((*error_)[5]);
      }
      else if (variable == "Hcurl-electric")
      {
        result = std::sqrt((*error_)[6]);
      }
      else if (variable == "Hcurl-magnetic")
      {
        result = std::sqrt((*error_)[7]);
      }
      else if (variable == "L2electric-post")
      {
        result = std::sqrt((*error_)[8]);
      }
      else if (variable == "L2electric-rel-post")
      {
        if ((*error_)[1] > 0.0)
          result = std::sqrt((*error_)[8] / (*error_)[1]);
        else
          FOUR_C_THROW(
              "Impossible to compute the post-processed electric relative error. The L2-norm of "
              "the analytical solution is zero, resulting in a division by zero.");
      }
      else if (variable == "Hdiv-electric-post")
      {
        result = std::sqrt((*error_)[9]);
      }
      else if (variable == "Hcurl-electric-post")
      {
        result = std::sqrt((*error_)[10]);
      }
      else
      {
        FOUR_C_THROW("Quantity '%s' not supported in result-test of elemagstic transport problems",
            variable.c_str());
      }

      nerr += compare_values(result, "NODE", container);
      test_count++;
    }
  }
}
FOUR_C_NAMESPACE_CLOSE
