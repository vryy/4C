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

#include "xfem_condition.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_gmsh_xfem_extension.H"
#include "../drt_geometry/integrationcell.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandle::InterfaceHandle(
    const Teuchos::RCP<DRT::Discretization>  xfemdis 
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

//! implement this member function in derived classes!
void XFEM::InterfaceHandle::toGmsh(const int step) const
{
  dserror ("not implemented for the InterfaceHandle base class");
  return;
}

//! implement this member function in derived classes!
int XFEM::InterfaceHandle::PositionWithinConditionNP(const LINALG::Matrix<3,1>& x_in) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}

//! implement this member function in derived classes!
int XFEM::InterfaceHandle::PositionWithinConditionN(const LINALG::Matrix<3,1>& x_in) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}


//! implement this member function in derived classes!
int XFEM::InterfaceHandle::PositionWithinConditionNP(const LINALG::Matrix<3,1>&     x_in,
                                                     GEO::NearestObject&  nearestobject) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}

//! implement this member function in derived classes!
int XFEM::InterfaceHandle::PositionWithinConditionN(const LINALG::Matrix<3,1>&     x_in,
                                                    GEO::NearestObject&  nearestobject) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::string XFEM::InterfaceHandle::toString() const
{
  std::stringstream s(" ");
  return s.str();
}

void XFEM::InterfaceHandle::SanityChecks() const
{
  //  std::cout << "numcuttedelements (elementalDomainIntCells_)   = " << elementalDomainIntCells_.size() << endl;
  //  std::cout << "numcuttedelements (elementalBoundaryIntCells_) = " << elementalBoundaryIntCells_.size() << endl;
  if (elementalDomainIntCells_.size() != elementalBoundaryIntCells_.size())
  {
    dserror("mismatch in cutted elements maps!");  
  }
  
  // sanity check, whether, we really have integration cells in the map
  for (std::map<int,GEO::DomainIntCells>::const_iterator 
      tmp = elementalDomainIntCells_.begin();
      tmp != elementalDomainIntCells_.end();
      ++tmp)
  {
    dsassert(tmp->second.empty() == false, "this is a bug!");
  }
  
  // sanity check, whether, we really have integration cells in the map
  for (std::map<int,GEO::BoundaryIntCells>::const_iterator 
      tmp = elementalBoundaryIntCells_.begin();
      tmp != elementalBoundaryIntCells_.end();
      ++tmp)
  {
    dsassert(tmp->second.empty() == false, "this is a bug!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::DomainIntCells XFEM::InterfaceHandle::GetDomainIntCells(
    const int gid,
    const DRT::Element::DiscretizationType distype
) const
{
  std::map<int,GEO::DomainIntCells>::const_iterator tmp = elementalDomainIntCells_.find(gid);
  if (tmp == elementalDomainIntCells_.end())
  {   
    // create default set with one dummy DomainIntCell of proper size
    GEO::DomainIntCells cells;
    cells.push_back(GEO::DomainIntCell(distype));
    return cells;
  }
  return tmp->second;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCells XFEM::InterfaceHandle::GetBoundaryIntCells(
    const int gid
) const
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
bool XFEM::InterfaceHandle::ElementIntersected(
    const int element_gid
) const
{
  if (elementalDomainIntCells_.find(element_gid) == elementalDomainIntCells_.end())
  {   
    return false;
  }
  else
  {
    return true;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool XFEM::InterfaceHandle::ElementHasLabel(
    const int element_gid,
    const int label
) const
{
  const GEO::BoundaryIntCells& bcells = elementalBoundaryIntCells().find(element_gid)->second;
  bool has_label = false;
  for (GEO::BoundaryIntCells::const_iterator bcell = bcells.begin(); bcell != bcells.end(); ++bcell)
  {
    const int surface_ele_gid = bcell->GetSurfaceEleGid();
    const int label_for_current_bele = labelPerElementId_.find(surface_ele_gid)->second;
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
    const int element_gid
) const
{
  std::set<int> labelset;
  const GEO::BoundaryIntCells& bcells = elementalBoundaryIntCells().find(element_gid)->second;
  for (GEO::BoundaryIntCells::const_iterator bcell = bcells.begin(); bcell != bcells.end(); ++bcell)
  {
    labelset.insert(bcell->GetSurfaceEleGid());
  }
  return labelset;
}


#endif  // #ifdef CCADISCRET
