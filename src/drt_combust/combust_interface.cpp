/*!-----------------------------------------------------------------------------------------------*
 \file combust_interface.cpp

 \brief interface handle that transports the intersection related things around for combustion problems

  detailed description in header file combust_interface.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/

#include "combust_interface.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_geometry/position_array.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"


/*------------------------------------------------------------------------------------------------*
 | constructor                                                                         henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::InterfaceHandleCombust::InterfaceHandleCombust(
    const Teuchos::RCP<DRT::Discretization> fluiddis,
    const Teuchos::RCP<const DRT::Discretization> gfuncdis,
    const Teuchos::RCP<const COMBUST::FlameFront> flamefront
    ) : xfemdis_(fluiddis),
        gfuncdis_(gfuncdis),
        flamefront_(flamefront)
{
  elementcutstatus_.clear();
}
/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::InterfaceHandleCombust::~InterfaceHandleCombust()
{
  return;
}

//! implement this function if needed for combustion!
void COMBUST::InterfaceHandleCombust::toGmsh(const int step) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = DRT::INPUT::IntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");

  const bool screen_out = true;

  //const bool gmsh_tree_output = false;

  const int myrank = xfemdis_->Comm().MyPID();

  //if (gmshdebugout)
  //{
  //  // debug: write both meshes to file in Gmsh format
  //  std::stringstream filename;
  //  std::stringstream filenamedel;
  //  filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".elements_coupled_system_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
  //  filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".elements_coupled_system_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
  //  std::remove(filenamedel.str().c_str());
  //  if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
  //  std::ofstream f_system(filename.str().c_str());
  //  IO::GMSH::XdisToStream("Fluid", 0.0, xfemdis_, elementalDomainIntCells_, elementalBoundaryIntCells_, f_system);
  //  //f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis_, cutterposnp_);
  //  f_system.close();
  //  if (screen_out) cout << " done" << endl;
  //}

  if (gmshdebugout)
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".domains_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".domains_" << std::setw(5) << setfill('0') << step-500 << ".p" << myrank << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";

    std::ofstream f_system(filename.str().c_str());
    {
      // stringstream for domains
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "Domains using CellCenter of Elements and Integration Cells \" {" << endl;

      for (int i=0; i<xfemdis_->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = xfemdis_->lRowElement(i);
        const GEO::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele);
        GEO::DomainIntCells::const_iterator cell;
        for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
        {
          const LINALG::SerialDenseMatrix& cellpos = cell->CellNodalPosXYZ();
          //const LINALG::Matrix<3,1> cellcenterpos(cell->GetPhysicalCenterPosition());
          int domain_id = 0;
          if (cell->getDomainPlus())
            domain_id = 1;
          //const double color = domain_id*100000+(closestElementId);
          const double color = domain_id;
          gmshfilecontent << IO::GMSH::cellWithScalarToString(cell->Shape(), color, cellpos) << endl;
        };
      };
      gmshfilecontent << "};" << endl;

      // the cut is done in element coordinates; for distorted elements, the sence of orientation can be reversed in physical coordinates
      gmshfilecontent << "View \" " << "Volume negative physical coordinates \" {" << endl;
      for (int i=0; i<xfemdis_->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = xfemdis_->lRowElement(i);
        const GEO::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele);
        GEO::DomainIntCells::const_iterator cell;
        for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
        {
          const LINALG::SerialDenseMatrix& cellpos = cell->CellNodalPosXYZ();
          const double volphys = cell->VolumeInPhysicalDomain();
          int vol_id = 3;
          if (volphys<=0)
            vol_id = 2;
          //const double color = domain_id*100000+(closestElementId);
          const double color = vol_id;
          gmshfilecontent << IO::GMSH::cellWithScalarToString(cell->Shape(), color, cellpos) << endl;
        };
      };
      gmshfilecontent << "};" << endl;

      f_system << gmshfilecontent.str();
    }
    f_system.close();
    if (screen_out) cout << " done" << endl;
  }

//  if (gmshdebugout) // print space time layer
//  {
//    std::stringstream filename;
//    std::stringstream filenamedel;
//    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".spacetime_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
//    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".spacetime_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
//    std::remove(filenamedel.str().c_str());
//    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
//
//    std::ofstream f_system(filename.str().c_str());
//    {
//      // stringstream for domains
//      stringstream gmshfilecontent;
//      gmshfilecontent << "View \" " << "SpaceTime cells \" {" << endl;
//      LINALG::SerialDenseVector vals(8);
//      vals(0) = 0.0;vals(1) = 0.0;vals(2) = 0.0;vals(3) = 0.0;
//      vals(4) = 1.0;vals(5) = 1.0;vals(6) = 1.0;vals(7) = 1.0;
//      for (std::map<int,XFEM::SpaceTimeBoundaryCell>::const_iterator slabiter = stlayer_.begin(); slabiter != stlayer_.end(); ++slabiter)
//      {
//        const XFEM::SpaceTimeBoundaryCell& slabitem = slabiter->second;
//
//        gmshfilecontent << IO::GMSH::cellWithScalarFieldToString(DRT::Element::hex8, vals, slabitem.get_xyzt()) << endl;
//      }
//      gmshfilecontent << "};" << endl;
//      f_system << gmshfilecontent.str();
//    }
//    f_system.close();
//    if (screen_out) cout << " done" << endl;
//  }
//
//
//  if (gmsh_tree_output)
//  {
//    // debug: write information about which structure we are in
//    std::stringstream filenameP;
//    std::stringstream filenamePdel;
//    filenameP    << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_points_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
//    filenamePdel << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_points_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
//    std::remove(filenamePdel.str().c_str());
//
//    std::cout << "writing " << left << std::setw(50) <<filenameP.str()<<"...";
//    std::ofstream f_systemP(filenameP.str().c_str());
//    {
//      // stringstream for cellcenter points
//      stringstream gmshfilecontentP;
//      gmshfilecontentP << "View \" " << "CellCenter of Elements and Integration Cells \" {" << endl;
//
//      for (int i=0; i<xfemdis_->NumMyRowElements(); ++i)
//      {
//        const DRT::Element* actele = xfemdis_->lRowElement(i);
//        const GEO::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele);
//        GEO::DomainIntCells::const_iterator cell;
//        for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
//        {
//          const LINALG::Matrix<3,1> cellcenterpos(cell->GetPhysicalCenterPosition());
//
//          //const int domain_id = PositionWithinConditionNP(cellcenterpos);
//
//          LINALG::SerialDenseMatrix point(3,1);
//          point(0,0)=cellcenterpos(0);
//          point(1,0)=cellcenterpos(1);
//          point(2,0)=cellcenterpos(2);
//
//          gmshfilecontentP << IO::GMSH::cellWithScalarToString(DRT::Element::point1, (actele->Id()), point) << endl;
//        };
//      };
//      gmshfilecontentP << "};" << endl;
//      f_systemP << gmshfilecontentP.str();
//    }
//    f_systemP.close();
//    cout << " done" << endl;
//
//    octTreenp_->printTree(DRT::Problem::Instance()->OutputControlFile()->FileName(), step);
//    octTreenp_->evaluateTreeMetrics(step);
//  }
  return;
}

/*------------------------------------------------------------------------------------------------*
 | fill integration cells according to current flame front                            henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::InterfaceHandleCombust::UpdateInterfaceHandle()
{
  elementalDomainIntCells_.clear();
  elementalBoundaryIntCells_.clear();
  elementcutstatus_.clear();

  elementalDomainIntCells_ = flamefront_->DomainIntCells();
  elementalBoundaryIntCells_ = flamefront_->BoundaryIntCells();
  elementcutstatus_ = flamefront_->ElementCutStatus();
  return;
}

/*----------------------------------------------------------------------*
 * to string
 *----------------------------------------------------------------------*/
std::string COMBUST::InterfaceHandleCombust::toString() const
{
  std::stringstream s(" ");
  return s.str();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::DomainIntCells COMBUST::InterfaceHandleCombust::GetDomainIntCells(
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
std::size_t COMBUST::InterfaceHandleCombust::GetNumDomainIntCells(
    const DRT::Element* xfemElement) const
{
  std::map<int,GEO::DomainIntCells>::const_iterator tmp = elementalDomainIntCells_.find(xfemElement->Id());
  if (tmp == elementalDomainIntCells_.end())
  {
    return 0;
  }
  return (tmp->second).size();
}


std::map<int, GEO::BoundaryIntCells> COMBUST::InterfaceHandleCombust::GetElementalBoundaryIntCells (
) const
{
  return elementalBoundaryIntCells_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCells COMBUST::InterfaceHandleCombust::GetBoundaryIntCells(
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
std::size_t COMBUST::InterfaceHandleCombust::GetNumBoundaryIntCells(
    const DRT::Element* xfemElement) const
{
  std::map<int,GEO::BoundaryIntCells>::const_iterator tmp = elementalBoundaryIntCells_.find(xfemElement->Id());
  if (tmp == elementalBoundaryIntCells_.end())
  {
    return 0;
  }
  return (tmp->second).size();
}
/*------------------------------------------------------------------------------------------------*
 | compute the volume of the minus domain (mass conservation check)               rasthofer 06/09 |
 *------------------------------------------------------------------------------------------------*/
const double COMBUST::InterfaceHandleCombust::ComputeVolumeMinus()
{
  //dserror("not conservative in parallel: loop over row elements -> get cells belonging to this element -> loop over cells");

  double volume = 0.0;

  for (int iele=0; iele<xfemdis_->NumMyRowElements(); ++iele)
  {
    const DRT::Element* ele = xfemdis_->lRowElement(iele);
    const GEO::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(ele);
    GEO::DomainIntCells::const_iterator itercell;
    for(itercell = elementDomainIntCells.begin(); itercell != elementDomainIntCells.end(); ++itercell )
    {
      // pick cells belonging to minus domain
      if (itercell->getDomainPlus() == false)
        volume += itercell->VolumeInPhysicalDomain();
    }
  }
  return volume;
}


/*------------------------------------------------------------------------------------------------*
 | compute surface of interface                                                   rasthofer 07/11 |
 |                                                                                    DA wichmann |
 *------------------------------------------------------------------------------------------------*/
const double COMBUST::InterfaceHandleCombust::ComputeSurface()
{
  double area = 0.0;

  for (int iele=0; iele<xfemdis_->NumMyRowElements(); ++iele)
  {
    const DRT::Element* ele = xfemdis_->lRowElement(iele);
    const GEO::BoundaryIntCells& elementBoundaryIntCells = this->GetBoundaryIntCells(ele->Id());
    GEO::BoundaryIntCells::const_iterator itercell;
    for(itercell = elementBoundaryIntCells.begin(); itercell != elementBoundaryIntCells.end(); ++itercell )
    {
      if (!(itercell->Shape() == DRT::Element::tri3 or itercell->Shape() == DRT::Element::quad4))
        dserror("invalid type of boundary integration cell for surface area calculation");

      // get coordinates of vertices defining flame front patch
      const LINALG::SerialDenseMatrix& coords = itercell->CellNodalPosXYZ();

      // first point of flame front patch
      LINALG::Matrix<3,1> pointA;
      pointA(0) = coords(0,0);
      pointA(1) = coords(1,0);
      pointA(2) = coords(2,0);

      // second point of flame front patch
      LINALG::Matrix<3,1> pointB;
      pointB(0) = coords(0,1);
      pointB(1) = coords(1,1);
      pointB(2) = coords(2,1);

      // first edge of flame front patch
      LINALG::Matrix<3,1> edgeBA;
      edgeBA.Update(1.0, pointA, -1.0, pointB);

      // third point of flame front patch
      LINALG::Matrix<3,1> pointC;
      pointC(0) = coords(0,2);
      pointC(1) = coords(1,2);
      pointC(2) = coords(2,2);

      // second edge of flame front patch
      LINALG::Matrix<3,1> edgeBC;
      edgeBC.Update(1.0, pointC, -1.0, pointB);

      LINALG::Matrix<3,1> crossP;
      crossP(0) = edgeBA(1)*edgeBC(2) - edgeBA(2)*edgeBC(1);
      crossP(1) = edgeBA(2)*edgeBC(0) - edgeBA(0)*edgeBC(2);
      crossP(2) = edgeBA(0)*edgeBC(1) - edgeBA(1)*edgeBC(0);

      if (itercell->Shape() == DRT::Element::tri3)
        area += crossP.Norm2() / 2.0;
      else
        area += crossP.Norm2();
    }
  }

  return area;
}


// return whether the element has a whole touched face and lies in the plus domain or not
bool COMBUST::InterfaceHandleCombust::ElementTouched(const int xfemeleid) const
{
  std::map<int,COMBUST::FlameFront::CutStatus>::const_iterator iter = elementcutstatus_.find(xfemeleid);
  if (iter != elementcutstatus_.end())
  {
    if (iter->second == COMBUST::FlameFront::touched)
      return true;
  }
  else
    dserror("cut status not available for ths element");

  return false;
}

// return whether the element is bisected or not
bool COMBUST::InterfaceHandleCombust::ElementBisected(const int xfemeleid) const
{
  std::map<int,COMBUST::FlameFront::CutStatus>::const_iterator iter = elementcutstatus_.find(xfemeleid);
  if (iter != elementcutstatus_.end())
  {
    if (iter->second == COMBUST::FlameFront::bisected)
      return true;
  }
  else
    dserror("cut status not available for this element");

  return false;
}

// return whether the element is trisected or not
bool COMBUST::InterfaceHandleCombust::ElementTrisected(const int xfemeleid) const
{
  std::map<int,COMBUST::FlameFront::CutStatus>::const_iterator iter = elementcutstatus_.find(xfemeleid);
  if (iter != elementcutstatus_.end())
  {
    if (iter->second == COMBUST::FlameFront::trisected)
      return true;
  }
  else
    dserror("cut status not available for ths element");

  return false;
}

// return whether the element is trisected or not
COMBUST::FlameFront::CutStatus COMBUST::InterfaceHandleCombust::ElementCutStatus(const int xfemeleid) const
{
  COMBUST::FlameFront::CutStatus cutstat = COMBUST::FlameFront::undefined;

  std::map<int,COMBUST::FlameFront::CutStatus>::const_iterator iter = elementcutstatus_.find(xfemeleid);
  if (iter != elementcutstatus_.end())
    cutstat = iter->second;

  return cutstat;
}

/*------------------------------------------------------------------------------------------------*
 | decide if this element is intersected (touched elements are intersected)                       |
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::InterfaceHandleCombust::ElementIntersected(
    const int element_gid) const
{
  if (elementalBoundaryIntCells_.find(element_gid) == elementalBoundaryIntCells_.end())
    return false;
  else
    return true;
}


/*------------------------------------------------------------------------------------------------*
 | decide if this element is split (= bi- or trisected or numerically bi- or trisected            |
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::InterfaceHandleCombust::ElementSplit(
    const DRT::Element* xfemElement) const
{
  bool bisected = false;

  std::size_t numcells = this->GetNumDomainIntCells(xfemElement);

  if (numcells > 1) // more than one domain integration cell -> element bisected or numerically touched
    bisected = true;
  else
    bisected = false;

  return bisected;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool COMBUST::InterfaceHandleCombust::ElementHasLabel(
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
 * implement this member function in derived classes!
 *----------------------------------------------------------------------*/
int COMBUST::InterfaceHandleCombust::PositionWithinConditionNP(const LINALG::Matrix<3,1>& x_in) const
{
  dserror("not implemented for the InterfaceHandle base class");
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::set<int> COMBUST::InterfaceHandleCombust::GetIntersectingBoundaryElementsGID(
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
std::set<int> COMBUST::InterfaceHandleCombust::LabelsPerElement(
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
void COMBUST::InterfaceHandleCombust::InvertElementsPerLabel()
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
