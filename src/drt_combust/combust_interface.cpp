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
#ifdef CCADISCRET

#include "combust_interface.H"

#include "../drt_lib/standardtypes_cpp.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_globalproblem.H"
 #include "../drt_lib/drt_utils.H"

// #include "../drt_io/io_gmsh.H"
// #include "../drt_io/io_gmsh_xfem_extension.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_gmsh_xfem_extension.H"


/*------------------------------------------------------------------------------------------------*
 | constructor                                                                         henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::InterfaceHandleCombust::InterfaceHandleCombust(
    const Teuchos::RCP<DRT::Discretization> fluiddis,
    const Teuchos::RCP<const DRT::Discretization> gfuncdis,
    const Teuchos::RCP<const COMBUST::FlameFront> flamefront
    ) : InterfaceHandle(fluiddis),
        gfuncdis_(gfuncdis),
        flamefront_(flamefront)
{
  //if (fluiddis->Comm().MyPID() == 0)
  //  std::cout << "Construct InterfaceHandleCombust" << std::endl;

  //URSULA
  /*
   * die DomainIntCells für alle Element werden jetzt in der FlameFront berechnet
   * das InterfaceHandle wird im Moment in der Sysmat verwendet um die Intergrationszellen
   * zu bekommen (aus ih->elementalDomainIntCells_) und um über die FlameFront die phi-Werte
   * zur Berechnung der Enrichmentfunction zu erhalten
   * das heißt: InterfaceHandle ist nur Verbindung zwischen FlameFront() und Sysmat und daher
   * eigentlich nicht zwingend notwendig
   *
   * dennoch könnte man das InterfaceHandle weiter verwenden
   * um bei der Sysmat nichts ändern zu müssen
   * dazu müssen die Integrationszellen von der FlameFront an das InterfaceHandle übergeben werden
   * elementalDomainIntCells_ = flamefront_->DomainIntCells();
   * und zwar immer dann nachdem ProcessFlameFront() aufgerufen wurde
   * daher wäre es auch besser ProcessFlameFront nur über das InterfaceHandle aufzurufen
   * und nicht wie bisher direkt im Combust-Algorithm
   * für den Konstruktor heißt das dann
   * elementalDomainIntCells_.clear();
   * flamefront_->ProcessFlameFront();
   * elementalDomainIntCells_ = flamefront_->DomainIntCells();
   * und die Funktion UpdateInterfaceHandle() erledigt das auch mit
   * anstelle das extra Aufrufs (ebenfalls in Combust-Algorithm)
   * elementalDomainIntCells_.clear();
   * flamefront_->ProcessFlameFront();
   * elementalDomainIntCells_ = flamefront_->DomainIntCells();
   */
  //NEIN, da in Constructor Fluidimpl... als Dummy verwendet
//  elementalDomainIntCells_.clear();
//  //flamefront_->ProcessFlameFront();
//  elementalDomainIntCells_ = flamefront_->DomainIntCells();
  //URSULA

  // Dinge, die hier passieren müssen, sind in diesen Funktionen zu finden:
  // computeIntersection
  // computePLC
  // computeCDT
  // Aufpassen! Die FlameFront enthält Infos über alle Fluid Col Elemente auf einem Proc!
  // Die Triangulierung sollte aber Row-mässig, d.h. eindeutig durchgeführt werden!
  // computeAverageGfuncValuePerIntCell   so wird bestimmt auf welcher Seite die Zelle liegt

  //if (fluiddis->Comm().MyPID() == 0)
  //  std::cout << "Construct InterfaceHandleCombust done" << std::endl;
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
  const bool gmshdebugout = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;

  const bool screen_out = true;

  //const bool gmsh_tree_output = false;

  const int myrank = xfemdis_->Comm().MyPID();

  if (gmshdebugout)
  {
    // debug: write both meshes to file in Gmsh format
    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".elements_coupled_system_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".elements_coupled_system_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());
    f_system << IO::GMSH::XdisToString("Fluid", 0.0, xfemdis_, elementalDomainIntCells_, elementalBoundaryIntCells_);
    //f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis_, cutterposnp_);
    f_system.close();
    if (screen_out) cout << " done" << endl;
  }

  if (gmshdebugout)
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".domains_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".domains_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
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
          const LINALG::Matrix<3,1> cellcenterpos(cell->GetPhysicalCenterPosition());
          int domain_id = 0;
          if (cell->getDomainPlus())
        	  domain_id = 1;
          //const double color = domain_id*100000+(closestElementId);
          const double color = domain_id;
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
  //uebergebe combustdyn und phinp an UpdateInterfaceHandle
  //rufe von dort ProcessFlameFront() und lasse die maps
  //elementalDomainIntCells_ und elementalBoundaryIntCells füllen
  //flamefront_->ProcessFlameFront(combustdyn, phinp, elementalDomainIntCells_, elementalBoundaryIntCells_);
  //flamefront_->ProcessFlameFront();
  elementalDomainIntCells_ = flamefront_->DomainIntCells();
  elementalBoundaryIntCells_ = flamefront_->BoundaryIntCells();
  return;
}

/*------------------------------------------------------------------------------------------------*
 | compute the volume of the minus domain (mass conservation check)               rasthofer 06/09 |
 *------------------------------------------------------------------------------------------------*/
const double COMBUST::InterfaceHandleCombust::ComputeVolumeMinus()
{
  double volume = 0.0;
  // loop over map entries for elements
  for(std::map<int,GEO::DomainIntCells>::iterator itermap=elementalDomainIntCells_.begin(); itermap!=elementalDomainIntCells_.end(); ++itermap)
  {
    // loop over integration cells for this element
    for(GEO::DomainIntCells::iterator itercell=itermap->second.begin(); itercell!=itermap->second.end(); itercell++)
    {
      // pick cells belonging to minus domain
      if (itercell->getDomainPlus() == false)
        volume += itercell->VolumeInPhysicalDomain();
    }
  }
  return volume;
}


//! implement this function if needed for combustion!
int COMBUST::InterfaceHandleCombust::PositionWithinConditionNP(const LINALG::Matrix<3,1>& x_in) const
{
  dserror("not implemented");
  return 0;
}

//! implement this function if needed for combustion!
int COMBUST::InterfaceHandleCombust::PositionWithinConditionN(const LINALG::Matrix<3,1>& x_in) const
{
  dserror("not implemented");
  return 0;
}

//! implement this function if needed for combustion!
int COMBUST::InterfaceHandleCombust::PositionWithinConditionNP(const LINALG::Matrix<3,1>&     x_in,
                              GEO::NearestObject&  nearestobject) const
{
  dserror("not implemented");
  return 0;
}

//! implement this function if needed for combustion!
int COMBUST::InterfaceHandleCombust::PositionWithinConditionN(const LINALG::Matrix<3,1>&     x_in,
                             GEO::NearestObject&  nearestobject) const
{
  dserror("not implemented");
  return 0;
}

#endif // #ifdef CCADISCRET
