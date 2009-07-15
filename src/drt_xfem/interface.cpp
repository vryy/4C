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
      octTreen_(rcp( new GEO::SearchTree(20)))
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
int XFEM::InterfaceHandle::PositionWithinConditionNP(const LINALG::Matrix<3,1>&     x_in,
                                                     GEO::NearestObject&  nearestobject) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}



/*----------------------------------------------------------------------*
 * implement this member function in derived classes!
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandle::PositionWithinConditionN(const LINALG::Matrix<3,1>&     x_in,
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


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool XFEM::InterfaceHandle::ElementIntersected(
    const int element_gid) const
{
  if (elementalDomainIntCells_.find(element_gid) == elementalDomainIntCells_.end())
    return false;
  else
    return true;
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
