/*!----------------------------------------------------------------------
\file fluid_rotsym_periodicbc_utils.cpp

\brief Methods to apply rotationally symmetric periodic boundary conditions
       for fluid problems on element level

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "fluid_rotsym_periodicbc_utils.H"
#include "../drt_lib/standardtypes_cpp.H"



using namespace Teuchos;

/*!
\brief  This class manages local transformations(rotation) of velocity fields

        o used for rotationally symmetric boundary conditions
        o used for rotationally symmetric boundary conditions

\author gjb

 */

/*----------------------------------------------------------------------*/
// this function is operating on global level and used
/*----------------------------------------------------------------------*/
double FLD::GetComponentOfRotatedVectorField
(const int idf,
    const RCP<const Epetra_Vector> proc0data,
    const int lid,
    const double rotangle
)
{
  switch (idf)
  {
  case 0:
  {
    // we assume that local dof id of y-component is lid+1
    double xvalue = (*proc0data)[lid];
    double yvalue = (*proc0data)[lid+1];
    return (xvalue*cos(rotangle) - yvalue*sin(rotangle));
    break;
  }
  case 1:
  {
    // we assume that local dof id of x-component is lid-1
    double xvalue = (*proc0data)[lid-1];
    double yvalue = (*proc0data)[lid];
    return (xvalue*sin(rotangle) + yvalue*(cos(rotangle)));
    break;
  }
  default: // we only allow rotation around z-axis!
    break;
  }

  return (*proc0data)[lid]; // case > 1: return unchanged value
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FLD::IsSlaveNodeOfRotSymPBC(
    const DRT::Node* node,
    double&          rotangle
)
{
  // get periodic surface/line boundary conditions
  vector<DRT::Condition*> pbc;
  node->GetCondition("SurfacePeriodic", pbc);
  if (pbc.empty())
    node->GetCondition("LinePeriodic", pbc);

  bool isrotsymslave(false);
  for (unsigned int j = 0; j < pbc.size();++j)
  {
    const std::string* isslave =
      pbc[j]->Get<std::string>("Is slave periodic boundary condition");
      if (*isslave == "Slave")
      {
        const double rotangle_deg = pbc[j]->GetDouble("Angle of rotation");
        if (abs(rotangle_deg) > EPS13) // angle is not zero
        {
          if (isrotsymslave)
            dserror("Node is slave of more than one rot.sym. periodic bc");
          isrotsymslave=true;
        }
        rotangle = rotangle_deg*PI/180.0; // angle of rotation (RAD)
      }
  }

      return isrotsymslave; // slave node with non-zero angle
    }


#endif  // #ifdef CCADISCRET
