/*!-----------------------------------------------------------------------------------------------*
 \file combust_interface.cpp

 \brief interface handle that transports the intersection related things around for combustion problems

  detailed description in header file combust_interface.H

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#include "combust_interface.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"
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
    const Teuchos::RCP<DRT::Discretization> gfuncdis
    ) : fluiddis_(fluiddis),
        gfuncdis_(gfuncdis)
{
  elementalDomainIntCells_.clear();
  elementalBoundaryIntCells_.clear();
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

  const int myrank = fluiddis_->Comm().MyPID();

  //if (gmshdebugout)
  //{
  //  // debug: write both meshes to file in Gmsh format
  //  std::stringstream filename;
  //  std::stringstream filenamedel;
  //  filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".elements_coupled_system_" << std::setw(5) << std::setfill('0') << step   << ".p" << myrank << ".pos";
  //  filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".elements_coupled_system_" << std::setw(5) << std::setfill('0') << step-5 << ".p" << myrank << ".pos";
  //  std::remove(filenamedel.str().c_str());
  //  if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
  //  std::ofstream f_system(filename.str().c_str());
  //  IO::GMSH::XdisToStream("Fluid", 0.0, fluiddis_, elementalDomainIntCells_, elementalBoundaryIntCells_, f_system);
  //  //f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis_, cutterposnp_);
  //  f_system.close();
  //  if (screen_out) std::cout << " done" << std::endl;
  //}

  if (gmshdebugout)
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".domains_" << std::setw(5) << std::setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".domains_" << std::setw(5) << std::setfill('0') << step-500 << ".p" << myrank << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << std::left << std::setw(50) <<filename.str()<<"...";

    std::ofstream f_system(filename.str().c_str());
    {
      // std::stringstream for domains
      std::stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "Domains using CellCenter of Elements and Integration Cells \" {" << std::endl;

      for (int i=0; i<fluiddis_->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = fluiddis_->lRowElement(i);
        const GEO::DomainIntCells& elementDomainIntCells = this->ElementDomainIntCells(actele->Id());
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
          gmshfilecontent << IO::GMSH::cellWithScalarToString(cell->Shape(), color, cellpos) << std::endl;
        };
      };
      gmshfilecontent << "};" << std::endl;

      // the cut is done in element coordinates; for distorted elements, the sence of orientation can be reversed in physical coordinates
      gmshfilecontent << "View \" " << "Volume negative physical coordinates \" {" << std::endl;
      for (int i=0; i<fluiddis_->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = fluiddis_->lRowElement(i);
        const GEO::DomainIntCells& elementDomainIntCells = this->ElementDomainIntCells(actele->Id());
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
          gmshfilecontent << IO::GMSH::cellWithScalarToString(cell->Shape(), color, cellpos) << std::endl;
        };
      };
      gmshfilecontent << "};" << std::endl;

      f_system << gmshfilecontent.str();
    }
    f_system.close();
    if (screen_out) std::cout << " done" << std::endl;
  }

//  if (gmshdebugout) // print space time layer
//  {
//    std::stringstream filename;
//    std::stringstream filenamedel;
//    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".spacetime_" << std::setw(5) << std::setfill('0') << step   << ".p" << myrank << ".pos";
//    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".spacetime_" << std::setw(5) << std::setfill('0') << step-5 << ".p" << myrank << ".pos";
//    std::remove(filenamedel.str().c_str());
//    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
//
//    std::ofstream f_system(filename.str().c_str());
//    {
//      // std::stringstream for domains
//      std::stringstream gmshfilecontent;
//      gmshfilecontent << "View \" " << "SpaceTime cells \" {" << std::endl;
//      LINALG::SerialDenseVector vals(8);
//      vals(0) = 0.0;vals(1) = 0.0;vals(2) = 0.0;vals(3) = 0.0;
//      vals(4) = 1.0;vals(5) = 1.0;vals(6) = 1.0;vals(7) = 1.0;
//      for (std::map<int,XFEM::SpaceTimeBoundaryCell>::const_iterator slabiter = stlayer_.begin(); slabiter != stlayer_.end(); ++slabiter)
//      {
//        const XFEM::SpaceTimeBoundaryCell& slabitem = slabiter->second;
//
//        gmshfilecontent << IO::GMSH::cellWithScalarFieldToString(DRT::Element::hex8, vals, slabitem.get_xyzt()) << std::endl;
//      }
//      gmshfilecontent << "};" << std::endl;
//      f_system << gmshfilecontent.str();
//    }
//    f_system.close();
//    if (screen_out) std::cout << " done" << std::endl;
//  }
//
//
//  if (gmsh_tree_output)
//  {
//    // debug: write information about which structure we are in
//    std::stringstream filenameP;
//    std::stringstream filenamePdel;
//    filenameP    << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_points_" << std::setw(5) << std::setfill('0') << step   << ".p" << myrank << ".pos";
//    filenamePdel << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_points_" << std::setw(5) << std::setfill('0') << step-5 << ".p" << myrank << ".pos";
//    std::remove(filenamePdel.str().c_str());
//
//    std::cout << "writing " << left << std::setw(50) <<filenameP.str()<<"...";
//    std::ofstream f_systemP(filenameP.str().c_str());
//    {
//      // std::stringstream for cellcenter points
//      std::stringstream gmshfilecontentP;
//      gmshfilecontentP << "View \" " << "CellCenter of Elements and Integration Cells \" {" << std::endl;
//
//      for (int i=0; i<fluiddis_->NumMyRowElements(); ++i)
//      {
//        const DRT::Element* actele = fluiddis_->lRowElement(i);
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
//          gmshfilecontentP << IO::GMSH::cellWithScalarToString(DRT::Element::point1, (actele->Id()), point) << std::endl;
//        };
//      };
//      gmshfilecontentP << "};" << std::endl;
//      f_systemP << gmshfilecontentP.str();
//    }
//    f_systemP.close();
//    std::cout << " done" << std::endl;
//
//    octTreenp_->printTree(DRT::Problem::Instance()->OutputControlFile()->FileName(), step);
//    octTreenp_->evaluateTreeMetrics(step);
//  }
  return;
}

/*------------------------------------------------------------------------------------------------*
 | fill integration cells according to current flame front                            henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::InterfaceHandleCombust::UpdateInterfaceHandle(
    std::map<int,GEO::DomainIntCells >&                         elementalDomainIntCells,
    std::map<int,GEO::BoundaryIntCells >&                       elementalBoundaryIntCells,
    std::map<int,COMBUST::InterfaceHandleCombust::CutStatus >&  elementcutstatus
)
{
  elementalDomainIntCells_.clear();
  elementalBoundaryIntCells_.clear();
  elementcutstatus_.clear();

  elementalDomainIntCells_ = elementalDomainIntCells;
  elementalBoundaryIntCells_ = elementalBoundaryIntCells;
  elementcutstatus_ = elementcutstatus;
  return;
}

/*----------------------------------------------------------------------*
 * to std::string
 *----------------------------------------------------------------------*/
std::string COMBUST::InterfaceHandleCombust::toString() const
{
  std::stringstream s(" ");
  return s.str();
}

//! return map of elemental domain integration cells
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------*
 | return number of domain integration cells for a given element                      henke 07/09 |
 *------------------------------------------------------------------------------------------------*/
std::map<int, GEO::DomainIntCells> COMBUST::InterfaceHandleCombust::DomainIntCells (
) const
{
  return elementalDomainIntCells_;
}


std::map<int, GEO::BoundaryIntCells> COMBUST::InterfaceHandleCombust::BoundaryIntCells (
) const
{
  return elementalBoundaryIntCells_;
}


std::map<int, COMBUST::InterfaceHandleCombust::CutStatus> COMBUST::InterfaceHandleCombust::CutState (
) const
{
  return elementcutstatus_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::DomainIntCells COMBUST::InterfaceHandleCombust::ElementDomainIntCells(
    const int gid) const
{
  std::map<int,GEO::DomainIntCells>::const_iterator tmp = elementalDomainIntCells_.find(gid);
  if (tmp == elementalDomainIntCells_.end())
  {
    // create default set with one dummy DomainIntCell of proper size
//    GEO::DomainIntCells cells;
//    cells.push_back(GEO::DomainIntCell(xfemElement->Shape(), GEO::InitialPositionArray(xfemElement)));
//    return cells;
    return GEO::DomainIntCells();
  }
  return tmp->second;
}


GEO::BoundaryIntCells COMBUST::InterfaceHandleCombust::ElementBoundaryIntCells(
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
COMBUST::InterfaceHandleCombust::CutStatus COMBUST::InterfaceHandleCombust::ElementCutStatus(const int xfemeleid) const
{
  COMBUST::InterfaceHandleCombust::CutStatus cutstat = COMBUST::InterfaceHandleCombust::undefined;

  std::map<int,COMBUST::InterfaceHandleCombust::CutStatus>::const_iterator iter = elementcutstatus_.find(xfemeleid);
  if (iter != elementcutstatus_.end())
    cutstat = iter->second;

  return cutstat;
}


std::size_t COMBUST::InterfaceHandleCombust::NumDomainIntCells(
    const int gid) const
{
  std::map<int,GEO::DomainIntCells>::const_iterator tmp = elementalDomainIntCells_.find(gid);
  if (tmp == elementalDomainIntCells_.end())
  {
    return 0;
  }
  return (tmp->second).size();
}


std::size_t COMBUST::InterfaceHandleCombust::NumBoundaryIntCells(
    const int gid) const
{
  std::map<int,GEO::BoundaryIntCells>::const_iterator tmp = elementalBoundaryIntCells_.find(gid);
  if (tmp == elementalBoundaryIntCells_.end())
  {
    return 0;
  }
  return (tmp->second).size();
}
/*------------------------------------------------------------------------------------------------*
 | compute the volume of the minus domain (mass conservation check)               rasthofer 06/09 |
 *------------------------------------------------------------------------------------------------*/
double COMBUST::InterfaceHandleCombust::ComputeVolumeMinus()
{
  double myvolume = 0.0;

  for (int iele=0; iele<fluiddis_->NumMyRowElements(); ++iele)
  {
    const DRT::Element* ele = fluiddis_->lRowElement(iele);
    const GEO::DomainIntCells& elementDomainIntCells = this->ElementDomainIntCells(ele->Id());
    GEO::DomainIntCells::const_iterator itercell;
    for(itercell = elementDomainIntCells.begin(); itercell != elementDomainIntCells.end(); ++itercell )
    {
      // pick cells belonging to minus domain
      if (itercell->getDomainPlus() == false)
        myvolume += itercell->VolumeInPhysicalDomain();
    }
  }

  double volume = 0.0;

  fluiddis_->Comm().SumAll(&myvolume,&volume,1);

  return volume;
}


/*------------------------------------------------------------------------------------------------*
 | compute surface of interface                                                   rasthofer 07/11 |
 |                                                                                    DA wichmann |
 *------------------------------------------------------------------------------------------------*/
double COMBUST::InterfaceHandleCombust::ComputeSurface()
{
  double myarea = 0.0;

  for (int iele=0; iele<fluiddis_->NumMyRowElements(); ++iele)
  {
    const DRT::Element* ele = fluiddis_->lRowElement(iele);
    const GEO::BoundaryIntCells& elementBoundaryIntCells = this->ElementBoundaryIntCells(ele->Id());
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
        myarea += crossP.Norm2() / 2.0;
      else
        myarea += crossP.Norm2();
    }
  }

  double area = 0.0;

  fluiddis_->Comm().SumAll(&myarea,&area,1);

  return area;
}


// return whether the element has a whole touched face and lies in the plus domain or not
bool COMBUST::InterfaceHandleCombust::ElementTouched(const int xfemeleid) const
{
  std::map<int,COMBUST::InterfaceHandleCombust::CutStatus>::const_iterator iter = elementcutstatus_.find(xfemeleid);
  if (iter != elementcutstatus_.end())
  {
    if (iter->second == COMBUST::InterfaceHandleCombust::touched)
      return true;
  }
  else
    dserror("cut status not available for ths element");

  return false;
}

// return whether the element is bisected or not
bool COMBUST::InterfaceHandleCombust::ElementBisected(const int xfemeleid) const
{
  std::map<int,COMBUST::InterfaceHandleCombust::CutStatus>::const_iterator iter = elementcutstatus_.find(xfemeleid);
  if (iter != elementcutstatus_.end())
  {
    if (iter->second == COMBUST::InterfaceHandleCombust::bisected)
      return true;
  }
  else
    dserror("cut status not available for this element");

  return false;
}

// return whether the element is trisected or not
bool COMBUST::InterfaceHandleCombust::ElementTrisected(const int xfemeleid) const
{
  std::map<int,COMBUST::InterfaceHandleCombust::CutStatus>::const_iterator iter = elementcutstatus_.find(xfemeleid);
  if (iter != elementcutstatus_.end())
  {
    if (iter->second == COMBUST::InterfaceHandleCombust::trisected)
      return true;
  }
  else
    dserror("cut status not available for ths element");

  return false;
}

// return whether the element is trisected or not
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
    const int gid) const
{
  bool bisected = false;

  std::size_t numcells = this->NumDomainIntCells(gid);

  if (numcells > 1) // more than one domain integration cell -> element bisected or numerically touched
    bisected = true;
  else
    bisected = false;

  return bisected;
}

