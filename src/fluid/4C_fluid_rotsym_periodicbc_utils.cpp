/*-----------------------------------------------------------*/
/*! \file

\brief Methods needed to apply rotationally symmetric periodic boundary
       conditions for fluid problems


\level 1

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_rotsym_periodicbc_utils.hpp"

#include "4C_lib_condition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_node.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::GetComponentOfRotatedVectorField(const int idf,
    const Teuchos::RCP<const Epetra_Vector> proc0data, const int lid, const double rotangle)
{
  switch (idf)
  {
    case 0:
    {
      // we assume that local dof id of y-component is lid+1
      double xvalue = (*proc0data)[lid];
      double yvalue = (*proc0data)[lid + 1];
      return (xvalue * cos(rotangle) - yvalue * sin(rotangle));
      break;
    }
    case 1:
    {
      // we assume that local dof id of x-component is lid-1
      double xvalue = (*proc0data)[lid - 1];
      double yvalue = (*proc0data)[lid];
      return (xvalue * sin(rotangle) + yvalue * (cos(rotangle)));
      break;
    }
    default:  // we only allow rotation around z-axis!
      break;
  }

  return (*proc0data)[lid];  // case > 1: return unchanged value
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FLD::IsSlaveNodeOfRotSymPBC(const DRT::Node* node, double& rotangle)
{
  // get periodic surface/line boundary conditions
  std::vector<DRT::Condition*> pbc;
  node->GetCondition("SurfacePeriodic", pbc);
  if (pbc.empty()) node->GetCondition("LinePeriodic", pbc);

  bool isrotsymslave(false);
  for (unsigned int j = 0; j < pbc.size(); ++j)
  {
    const std::string& isslave =
        pbc[j]->parameters().Get<std::string>("Is slave periodic boundary condition");
    if (isslave == "Slave")
    {
      rotangle = GetRotAngleFromCondition(pbc[j]);
      if (abs(rotangle) > 1e-13)  // angle is not zero
      {
        if (isrotsymslave) FOUR_C_THROW("Node is slave of more than one rot.sym. periodic bc");
        isrotsymslave = true;
      }
    }
  }

  return isrotsymslave;  // yes, it is a slave node with non-zero angle
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::GetRotAngleFromCondition(const DRT::Condition* cond)
{
  const double rotangle_deg = cond->parameters().Get<double>("Angle of rotation");

  return rotangle_deg * M_PI / 180.0;  // angle of rotation (RAD);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::GetRelevantSlaveNodesOfRotSymPBC(
    std::map<int, double>& pbcslavenodemap, Teuchos::RCP<DRT::Discretization> dis)
{
  // get all periodic boundary conditions
  std::vector<DRT::Condition*> mypbccond;
  dis->GetCondition("SurfacePeriodic", mypbccond);
  if (mypbccond.empty())
  {
    dis->GetCondition("LinePeriodic", mypbccond);
  }

  // loop the periodic boundary conditions
  for (unsigned numcond = 0; numcond < mypbccond.size(); ++numcond)
  {
    const std::string& mymasterslavetoggle =
        mypbccond[numcond]->parameters().Get<std::string>("Is slave periodic boundary condition");
    const double rotangle = FLD::GetRotAngleFromCondition(mypbccond[numcond]);

    // only slave nodes with non-zero angle of rotation require rotation
    // of vector results!
    if ((mymasterslavetoggle == "Slave") && (abs(rotangle) > 1e-13))
    {
      const std::vector<int>* nodes = mypbccond[numcond]->GetNodes();
      for (unsigned int inode = 0; inode < nodes->size(); inode++)
      {
        const int nodegid = nodes->at(inode);
        pbcslavenodemap[nodegid] = rotangle;
      }
    }  // end is slave condition?
  }    // end loop periodic boundary conditions
}

FOUR_C_NAMESPACE_CLOSE
