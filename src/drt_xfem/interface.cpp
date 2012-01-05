/*!
\file interface.cpp

\brief provides a class that represents an enriched physical scalar field

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
#ifdef CCADISCRET

#include "interface.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_gmsh_xfem_extension.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_geometry/position_array.H"



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandle::InterfaceHandle(
    const Teuchos::RCP<DRT::Discretization>&  xfemdis
    ) :
      xfemdis_(xfemdis),
      octTreenp_(rcp( new GEO::SearchTree(20))),
      octTreen_(rcp( new GEO::SearchTree(1))) //TODO: searchtree does not return nearest object for max-depth=20. Find out. This is needed for time step update, where I extrapolate from the nearest object...
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandle::~InterfaceHandle()
{
    return;
}


/*----------------------------------------------------------------------*
 * implement this member function in derived classes!
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandle::toGmsh(const int step) const
{
  dserror ("not implemented for the InterfaceHandle base class");
  return;
}


/*----------------------------------------------------------------------*
 * implement this member function in derived classes!
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandle::PositionWithinConditionNP(const LINALG::Matrix<3,1>& x_in) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}


/*----------------------------------------------------------------------*
 * implement this member function in derived classes!
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandle::PositionWithinConditionN(const LINALG::Matrix<3,1>& x_in) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}


/*----------------------------------------------------------------------*
 * implement this member function in derived classes!
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandle::PositionWithRespectToInterfaceNP(
    const LINALG::Matrix<3,1>& x_in,
    const int label) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}


/*----------------------------------------------------------------------*
 * implement this member function in derived classes!
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandle::PositionWithRespectToInterfaceN(
    const LINALG::Matrix<3,1>& x_in,
    const int label) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}


/*----------------------------------------------------------------------*
 * implement this member function in derived classes!
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandle::PositionWithinConditionNP(
    const LINALG::Matrix<3,1>&     x_in,
    GEO::NearestObject&  nearestobject) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}



/*----------------------------------------------------------------------*
 * implement this member function in derived classes!
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandle::PositionWithinConditionN(
    const LINALG::Matrix<3,1>&     x_in,
    GEO::NearestObject&  nearestobject) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}


/*----------------------------------------------------------------------*
 * implement this member function in derived classes!
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandle::PositionWithRespectToInterfaceNP(
    const LINALG::Matrix<3,1>&     x_in,
    const int label,
    GEO::NearestObject&  nearestobject) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}



/*----------------------------------------------------------------------*
 * implement this member function in derived classes!
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandle::PositionWithRespectToInterfaceN(
    const LINALG::Matrix<3,1>&     x_in,
    const int label,
    GEO::NearestObject&  nearestobject) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}


/*----------------------------------------------------------------------*
 * to string
 *----------------------------------------------------------------------*/
std::string XFEM::InterfaceHandle::toString() const
{
  std::stringstream s(" ");
  return s.str();
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::DomainIntCells XFEM::InterfaceHandle::GetDomainIntCells(
    const DRT::Element* xfemElement) const
{
  std::map<int,GEO::DomainIntCells>::const_iterator tmp = elementalDomainIntCells_.find(xfemElement->Id());
  if (tmp == elementalDomainIntCells_.end())
  {
    // create default set with one dummy DomainIntCell of proper size
    GEO::DomainIntCells cells;
    cells.push_back(GEO::DomainIntCell(xfemElement->Shape(), GEO::InitialPositionArray(xfemElement)));
    return cells;
  }
  return tmp->second;
}


/*------------------------------------------------------------------------------------------------*
 | return number of domain integration cells for a given element                      henke 07/09 |
 *------------------------------------------------------------------------------------------------*/
std::size_t XFEM::InterfaceHandle::GetNumDomainIntCells(
    const DRT::Element* xfemElement) const
{
  std::map<int,GEO::DomainIntCells>::const_iterator tmp = elementalDomainIntCells_.find(xfemElement->Id());
  if (tmp == elementalDomainIntCells_.end())
  {
    return 0;
  }
  return (tmp->second).size();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCells XFEM::InterfaceHandle::GetBoundaryIntCells(
    const int gid) const
{
  std::map<int,GEO::BoundaryIntCells>::const_iterator tmp = elementalBoundaryIntCells_.find(gid);
  if (tmp == elementalBoundaryIntCells_.end())
  {
    // return empty list
    return GEO::BoundaryIntCells();
  }
  return tmp->second;
}


/*------------------------------------------------------------------------------------------------*
 | return number of boundary integration cells for a given element                     henke 06/10 |
 *------------------------------------------------------------------------------------------------*/
std::size_t XFEM::InterfaceHandle::GetNumBoundaryIntCells(
    const DRT::Element* xfemElement) const
{
  std::map<int,GEO::BoundaryIntCells>::const_iterator tmp = elementalBoundaryIntCells_.find(xfemElement->Id());
  if (tmp == elementalBoundaryIntCells_.end())
  {
    return 0;
  }
  return (tmp->second).size();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::set<int> XFEM::InterfaceHandle::GetAvailableBoundaryLabels() const
{
  std::set<int> labels;
  for(std::map<int,std::set<int> >::const_iterator conditer = boundaryElementsByLabel_.begin();
      conditer!=boundaryElementsByLabel_.end();
      ++conditer)
  {
    labels.insert(conditer->first);
  }
  return labels;
}


/*------------------------------------------------------------------------------------------------*
 | decide if this element is intersected (touched elements are intersected)                       |
 *------------------------------------------------------------------------------------------------*/
bool XFEM::InterfaceHandle::ElementIntersected(
    const int element_gid) const
{
  if (elementalBoundaryIntCells_.find(element_gid) == elementalBoundaryIntCells_.end())
    return false;
  else
    return true;
}


/*------------------------------------------------------------------------------------------------*
 | decide if this element is bisected (truely split into two parts, not intersected)              |
 *------------------------------------------------------------------------------------------------*/
bool XFEM::InterfaceHandle::ElementSplit(
    const DRT::Element* xfemElement) const
{
  bool bisected = false;

  std::size_t numcells = this->GetNumDomainIntCells(xfemElement);

  if (numcells > 1) // more than one domain integration cell -> element bisected or numerically touched
  {
    if (TouchedType(xfemElement)==0) // not touched
      bisected = true;
    else
      bisected = false;
  }
  else
    bisected = false;

  return bisected;
}


/*------------------------------------------------------------------------------------------------*
 | decide if this element is touched on a whole 2D face                                           |
 *------------------------------------------------------------------------------------------------*/
bool XFEM::InterfaceHandle::ElementTouched(
    const DRT::Element* xfemElement) const
{
  if (TouchedType(xfemElement)!=0) // not untouched = touched
    return true;
  else
    return false;
}


/*------------------------------------------------------------------------------------------------*
 | decide if this element is touched on a whole 2D face and lies in the plus domain               |
 *------------------------------------------------------------------------------------------------*/
bool XFEM::InterfaceHandle::ElementTouchedPlus(
    const DRT::Element* xfemElement) const
{
  if (TouchedType(xfemElement)==+1)
    return true;
  else
    return false;
}


/*------------------------------------------------------------------------------------------------*
 | decide if this element is touched on a whole 2D face and lies in the minus domain              |
 *------------------------------------------------------------------------------------------------*/
bool XFEM::InterfaceHandle::ElementTouchedMinus(
    const DRT::Element* xfemElement) const
{
  if (TouchedType(xfemElement)==-1)
    return true;
  else
    return false;
}


/*------------------------------------------------------------------------------------------------*
 | decide if this element is numerically touched                                                  |
 *------------------------------------------------------------------------------------------------*/
int XFEM::InterfaceHandle::TouchedType(
    const DRT::Element* xfemElement) const
{
  int touchedType = 0; // untouched

  std::size_t numDomainCells = this->GetNumDomainIntCells(xfemElement);
  std::size_t numBoundaryCells = this->GetNumBoundaryIntCells(xfemElement);

  if (numDomainCells == 1 and numBoundaryCells == 1) // real touched element, touchedminus or touchedplus possible
  {
    // reinitialization (just for sure)
    touchedType = 0;

    // get the boundary integration cell (there is exactly one)
    GEO::BoundaryIntCells boundarycells = elementalBoundaryIntCells_.find(xfemElement->Id())->second;
    if(boundarycells.size() != 1) dserror("there must be exactly one boundary cell");

    for(GEO::BoundaryIntCells::const_iterator iter=boundarycells.begin(); iter!=boundarycells.end(); iter++)
    {
      if (iter->getDomainPlus())
        touchedType = 1;
      else
        touchedType = -1;
    }
  }
  else if (numDomainCells > 1) // potential numerically touched
  {
    // reinitialization (just for sure)
    touchedType = 0;

    GEO::DomainIntCells domaincells = elementalDomainIntCells_.find(xfemElement->Id())->second;
    if (domaincells.size()<1) dserror("there must be at least one domain cell");

    for(GEO::DomainIntCells::const_iterator iter=domaincells.begin(); iter!=domaincells.end(); iter++)
    {
      if (iter->getDomainPlus())
      {
        if (touchedType==0) // first time modified
          touchedType = +1; // potentially touched plus
        else if (touchedType==-1) // domain cell in omega+ and omega- ->not touched
        {
          touchedType = 0;
          break;
        } // if touched = 1, the element was and is again potentially touched plus
      }
      else
      {
        if (touchedType==0) // first time modified
          touchedType = -1; // potentially touched minus
        else if (touchedType==+1) // domain cell in omega+ and omega- ->not touched
        {
          touchedType = 0;
          break;
        }
        // if touched = -1, the element was and is again potentially touched minus
      }
    }
  }
  else
    touchedType = 0; // untouched

  return touchedType;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool XFEM::InterfaceHandle::ElementHasLabel(
    const int element_gid,
    const int label) const
{
  if(elementalBoundaryIntCells_.empty())
    dserror("boundary intcells are empty");

  const GEO::BoundaryIntCells& bcells = elementalBoundaryIntCells_.find(element_gid)->second;
  bool has_label = false;
  for (GEO::BoundaryIntCells::const_iterator bcell = bcells.begin(); bcell != bcells.end(); ++bcell)
  {
    const int surface_ele_gid = bcell->GetSurfaceEleGid();
    const int label_for_current_bele = labelPerBoundaryElementId_.find(surface_ele_gid)->second;
    if (label == label_for_current_bele)
    {
      has_label = true;
      break;
    }
  }
  return has_label;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::set<int> XFEM::InterfaceHandle::GetIntersectingBoundaryElementsGID(
    const int element_gid
    ) const
{
  std::set<int> begids;

  if(elementalBoundaryIntCells_.empty())
    dserror("boundary intcells are empty");

  const GEO::BoundaryIntCells& bcells = elementalBoundaryIntCells_.find(element_gid)->second;
  for (GEO::BoundaryIntCells::const_iterator bcell = bcells.begin(); bcell != bcells.end(); ++bcell)
  {
    begids.insert(bcell->GetSurfaceEleGid());
  }
  return begids;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::set<int> XFEM::InterfaceHandle::LabelsPerElement(
    const int element_gid) const
{
  std::set<int> labelset;
  if(elementalBoundaryIntCells_.empty())
    return labelset;

  const GEO::BoundaryIntCells& bcells = elementalBoundaryIntCells_.find(element_gid)->second;
  for (GEO::BoundaryIntCells::const_iterator bcell = bcells.begin(); bcell != bcells.end(); ++bcell)
  {
    labelset.insert(bcell->GetSurfaceEleGid());
  }
  return labelset;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandle::InvertElementsPerLabel()
{
  labelPerBoundaryElementId_.clear();
  for(std::map<int,std::set<int> >::const_iterator conditer = boundaryElementsByLabel_.begin();
      conditer!=boundaryElementsByLabel_.end();
      ++conditer)
  {
    const int xfemlabel = conditer->first;
    for(std::set<int>::const_iterator eleiditer = conditer->second.begin(); eleiditer!=conditer->second.end(); ++eleiditer)
    {
      const int eleid = *eleiditer;
      if (labelPerBoundaryElementId_.count(eleid) == 1)
        dserror("Assumption violation: there should be exactly ONE xfem condition per boundary element id!");
      labelPerBoundaryElementId_[eleid] = xfemlabel;
    }
  }
  return;
}



#endif  // #ifdef CCADISCRET
