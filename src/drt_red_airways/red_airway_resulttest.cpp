/*---------------------------------------------------------------------*/
/*! \file

\brief Testing of calculation results for reduced elements


\level 3

*/
/*---------------------------------------------------------------------*/


#include "airwayimplicitintegration.H"
#include "red_airway_resulttest.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
AIRWAY::RedAirwayResultTest::RedAirwayResultTest(RedAirwayImplicitTimeInt& airways)
    : DRT::ResultTest("RED_AIRWAY")
{
  dis_ = airways.Discretization();
  mynodesol_pressure_ = airways.Pnp();
  mynodesol_flow_in_ = airways.Qin_np();
  mynodesol_flow_out_ = airways.Qout_np();

  myelemsol_pressure_external_ = airways.Pext_np();
  myelemsol_acinivol_ = airways.AciniVolume();
  myelemsol_airwayvol_ = airways.AirwayVolume();
  myelemsol_open_ = airways.Open();
  myelemsol_openingTrajectory_ = airways.OpeningTrajectory();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AIRWAY::RedAirwayResultTest::TestNode(
    DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != dis_->Name()) return;

  int node;
  res.ExtractInt("NODE", node);
  node -= 1;

  int havenode(dis_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  dis_->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    dserror("Node %d does not belong to discretisation %s", node + 1, dis_->Name().c_str());
  }
  else
  {
    if (dis_->HaveGlobalNode(node))
    {
      DRT::Node* actnode = dis_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != dis_->Comm().MyPID()) return;

      double result = 0.;
      const Epetra_BlockMap& nodemap = mynodesol_pressure_->Map();
      std::string position;
      res.ExtractString("QUANTITY", position);

      // test result value of single scalar field
      if (position == "pressure")
      {
        result = (*mynodesol_pressure_)[nodemap.LID(dis_->Dof(actnode, 0))];
      }
      else if (position == "flow_in")
      {
        result = (*mynodesol_flow_in_)[nodemap.LID(dis_->Dof(actnode, 0))];
      }
      else if (position == "flow_out")
      {
        result = (*mynodesol_flow_out_)[nodemap.LID(dis_->Dof(actnode, 0))];
      }
      // test result values for a system of scalars
      else
      {
        dserror(
            "Quantity '%s' not supported in result-test of red_airway problems", position.c_str());
      }

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
}


/*----------------------------------------------------------------------*
 * Element based result test for red-airway problems. Tests the results *
 * of acini_volume.                                                     *
 *                                                         roth 11/2014 *
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayResultTest::TestElement(
    DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != dis_->Name()) return;

  int element;
  res.ExtractInt("ELEMENT", element);
  element -= 1;

  int haveelement(dis_->HaveGlobalElement(element));
  int iselementofanybody(0);
  dis_->Comm().SumAll(&haveelement, &iselementofanybody, 1);

  if (iselementofanybody == 0)
  {
    dserror("Node %d does not belong to discretisation %s", element + 1, dis_->Name().c_str());
  }
  else
  {
    if (dis_->HaveGlobalElement(element))
    {
      const DRT::Element* actelement = dis_->gElement(element);

      // Here we are just interested in the elements that we own (i.e. a row element)!
      if (actelement->Owner() != dis_->Comm().MyPID()) return;

      double result = 0.;
      const Epetra_BlockMap& elementmap = myelemsol_acinivol_->Map();
      std::string position;
      res.ExtractString("QUANTITY", position);
      if (position == "pressure_external")
      {
        result = (*myelemsol_pressure_external_)[elementmap.LID(actelement->Id())];
      }
      else if (position == "acini_volume")
      {
        result = (*myelemsol_acinivol_)[elementmap.LID(actelement->Id())];
      }
      else if (position == "airway_volume")
      {
        result = (*myelemsol_airwayvol_)[elementmap.LID(actelement->Id())];
      }
      else if (position == "opening_status")
      {
        result = (*myelemsol_open_)[elementmap.LID(actelement->Id())];
      }
      else if (position == "opening_trajectory")
      {
        result = (*myelemsol_openingTrajectory_)[elementmap.LID(actelement->Id())];
      }
      else
      {
        dserror(
            "Quantity '%s' not supported in result-test of red_airway problems.", position.c_str());
      }

      nerr += CompareValues(result, "ELEMENT", res);
      test_count++;
    }
  }
}
