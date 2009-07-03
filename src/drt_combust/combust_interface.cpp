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
// #include "../drt_geometry/integrationcell.H"


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
  if (fluiddis->Comm().MyPID() == 0)
    std::cout << "Construct InterfaceHandleCombust" << std::endl;

/* Ich muss erstmal schauen, ob die DomainIntCell für meine Zwecke sinnvoll ist. Falls nicht, steht
 * das ganze InterfaceHandle in Frage. Es könnte auch alles in die FlameFront integriert werden.
 * 
 * henke 03/09 */
  
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
  
  std::cout << "Proc " << fluiddis->Comm().MyPID() << ": Hier passiert absolut nichts" << std::endl;
  // Dinge, die hier passieren müssen, sind in diesen Funktionen zu finden:
  // computeIntersection
  // computePLC
  // computeCDT
  // Aufpassen! Die FlameFront enthält Infos über alle Fluid Col Elemente auf einem Proc!
  // Die Triangulierung sollte aber Row-mässig, d.h. eindeutig durchgeführt werden!
  // computeAverageGfuncValuePerIntCell   so wird bestimmt auf welcher Seite die Zelle liegt

  if (fluiddis->Comm().MyPID() == 0)
    std::cout << "Construct InterfaceHandleCombust done" << std::endl;

  // Kläre was mit diesen Bäumen passieren soll!
  //  octTreenp_ = rcp( new GEO::SearchTree(5));
  //  octTreenp_->initializeTree(AABB, elementsByLabel_, GEO::TreeType(GEO::OCTTREE));
  //  octTreen_ = rcp( new GEO::SearchTree(5));
  //  octTreen_->initializeTree(AABB, elementsByLabel_, GEO::TreeType(GEO::OCTTREE));
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
  dserror ("not implemented");
  return;
}

/*------------------------------------------------------------------------------------------------*
 | fill integration cells according to current flame front                            henke 10/08 | 
 *------------------------------------------------------------------------------------------------*/
void COMBUST::InterfaceHandleCombust::UpdateInterfaceHandle()
{
  elementalDomainIntCells_.clear();
  elementalBoundaryIntCells_.clear();
  //flamefront_->ProcessFlameFront();
  elementalDomainIntCells_ = flamefront_->DomainIntCells();
  elementalBoundaryIntCells_ = flamefront_->BoundaryIntCells();
  return;
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

//URSULA
/*------------------------------------------------------------------------------------------------*
 | determine which side of the interface an integration cell belongs to; the domain is devided    |
 | into a "+"-domain and a "-"-domain, separated by the level set function        rasthofer 07/09 |
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::InterfaceHandleCombust::GetIntCellPosition(
		const GEO::DomainIntCell& intcell,
		const DRT::Element*       xfemele
		)
{
	bool inphiplus = false;
	
	//berechnen zunächst phi an jedem Knoten der DomainIntCell
	//anschließen wird der Mittelwert gebildet
	//falls Mittelwert
	// >0 true
	// <0 false
	
	const LINALG::SerialDenseMatrix IntCellCoord = intcell.CellNodalPosXiDomain();
	const DRT::Element::DiscretizationType cell_distype = intcell.Shape();
	int numcellnodes = 0;
	switch(cell_distype){
	case DRT::Element::tet4:
	{
		numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;
		break;
	}
	case DRT::Element::hex8:
	{
		numcellnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
		break;
	}
	default:
		dserror("Discretization Type (intcell) not supported yet!");
	}
	
	//const DRT::Element::DiscretizationType xfem_distype = xfemele->Shape();
	const int numelenodes = xfemele->NumNode();
	
	//get phi at the nodes of the XFEM-Element
    const Teuchos::RCP<Epetra_Vector>  phinp = flamefront_->Phinp();
    // remark: vector "lm" is neccessary, because ExtractMyValues() only accepts "vector<int>"
    // arguments, but ele->NodeIds delivers an "int*" argument
    vector<int> lm(numelenodes);
    // get vector of node GIDs for this element
    const int* nodeids = xfemele->NodeIds();
    for (unsigned inode=0; inode < lm.size(); inode++)
      lm[inode] = nodeids[inode];
    // create vector "phi" holding phi values for this element
    vector<double> phi(numelenodes,1000.0);
    // get entries in "phinp" corresponding to node GIDs "lm" and store them in "phi"
    DRT::UTILS::ExtractMyValues(*phinp,phi,lm);
    
    //calculate phi at the nodes of the integrationcell
	std::vector<double> phicellnodes (numcellnodes);
	for (int icellnode=0; icellnode<numcellnodes; icellnode++)
		phicellnodes[icellnode] = 0.0;
	
	for (int icellnode=0; icellnode<numcellnodes; icellnode++)
	{
		//calculate shape function at IntCell node
		Epetra_SerialDenseVector  funct(numelenodes);
	    DRT::UTILS::shape_function_3D(funct,IntCellCoord(0,icellnode),IntCellCoord(1,icellnode),IntCellCoord(2,icellnode),xfemele->Shape());
	    /*
	     *         ___
	     * gval(x)=\   N(x)*gvali
	     *         /
	     *         ___
	     */
		for (int ielenode=0; ielenode<numelenodes; ielenode++)
		{
			phicellnodes[icellnode] += funct(ielenode) * phi[ielenode];
		}
	}
	
	//calculate average Gfunc value
	double averagePhiValue = 0.0;
	for (int icellnode=0; icellnode<numcellnodes; icellnode++)
		averagePhiValue += phicellnodes[icellnode];
	averagePhiValue /= numcellnodes;
	
	//determine DomainIntCell position
	if(averagePhiValue>0)
		inphiplus = true;
	if(averagePhiValue==0)
		dserror("can't determine DomainIntCell position");
	
	return inphiplus;
}

#endif // #ifdef CCADISCRET
