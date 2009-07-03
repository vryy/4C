/*!-----------------------------------------------------------------------------------------------*
 \file combust_flamefront.cpp

 \brief

  detailed description in header file combust_interface.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "combust_flamefront.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/linalg_utils.H"  // LINALG::Export
#include <Teuchos_StandardParameterEntryValidators.hpp>
// #include "../drt_geometry/integrationcell.H"

//URSULA
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
//URSULA


// extern struct _FILES  allfiles;

/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::FlameFront::FlameFront(
    const Teuchos::RCP<const DRT::Discretization> fluiddis,
    const Teuchos::RCP<DRT::Discretization> gfuncdis
    ) :
    fluiddis_(fluiddis),
    gfuncdis_(gfuncdis),
    phinp_(Teuchos::null)
{
  if (fluiddis->Comm().MyPID() == 0)
    std::cout << "Constructing FlameFront done" << std::endl;
}
/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::FlameFront::~FlameFront()
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 | public: flame front control routine                                                henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ProcessFlameFront(
       const Teuchos::ParameterList& combustdyn,
       const Teuchos::RCP<Epetra_Vector> phinp)
{
  /* This function is accessible from outside the FlameFront class. It can be called e.g. by the
   * interface handle constructor. If the FlameFront class will always do the same thing, this
   * function could be spared. Then the content of this function could be added to the constructor
   * of this class.
   *
   * henke 10/08
   */

  /* In the following processing of the flame front the fluid discretization is cut by the level set
   * function (G-function). This is the reason why fluid elements need to be able to access the
   * corresponding scalar G-function values for all their nodes. Therefore, the ScaTra dof-based
   * vector phinp has to be rearranged in a parallel environment to represent a Fluid node-based
   * vector. This involves two steps:
   * 1. ScaTra DofRowMap -> Fluid NodeRowMap
   * 2. Fluid NodeRowMap -> Fluid NodeColMap
   *
   * henke 02/09
   */
	
  //URSULA
  //diese Maps müssen bei jeden Durchlauf von ProcessFlameFront
  //neu gefüllt werden
  elementintcells_.clear();
  boundaryintcells_.clear();
  //URSULA

  //------------------------------------------------------------------------------------------------
  // Rearranging vector phinp from gfuncdis DofRowMap to fluiddis NodeRowMap
  //------------------------------------------------------------------------------------------------
  const Teuchos::RCP<Epetra_Vector> phinprow = rcp(new Epetra_Vector(*fluiddis_->NodeRowMap()));

  if (phinprow->MyLength() != phinp->MyLength())
    dserror("vectors phinprow and phinp must have the same length");

  //----------------------------------------------------------------------------------------------
  // congruent (matching) discretizations (Fluid == G-function)
  //----------------------------------------------------------------------------------------------
  if(true)
  {
    /* The premise of geometrically identical discretizations leads to a significant simplification
     * of the rearrangment process. If nodes of both discretizations are distributed in the same way
     * (which they are) and the dofs belonging to nodes lie on the same processors, then the
     * following holds:
     * G-function NodeMap   identical to   Fluid NodeMap
     * G-function NodeMap     similar to   G-function DofMap   and therefore
     * G-function DofMap      similar to   Fluid NodeMap
     *
     * This means that vector phinp living on the G-function DofRowMap can be directly copied to the
     * Fluid NodeRowMap without any rearrangement for congruent discretizations.
     *
     * henke 02/09
     */

    // get the G-function values corresponding to fluid nodes
    *phinprow = *phinp;
    /* ausgeschrieben, was *phinprow = *phinp in einer Zeile macht:
       for(int lnodeid=0;lnodeid<fluiddis_->NumMyRowNodes();lnodeid++)
       {
         // get the G-function value corresponding to this fluid node
         double value = (*phinp)[lnodeid];
         phinprow->ReplaceMyValue(lnodeid,value);
       }
     */
  }
  //----------------------------------------------------------------------------------------------
  // no congruent (matching) discretizations (Fluid != G-function)
  //----------------------------------------------------------------------------------------------
  else
  {
    /* In general this operation is very complex. It involves parallel communication and
     * interpolation of the G-function at fluid nodes, if non-matching discretizations for fluid
     * and G-function are used.
     *
     * henke 02/09
     */
    dserror("non-congruent discretizations cannot be handled yet!");
  }

  //------------------------------------------------------------------------------------------------
  // Exporting vector phinp from fluiddis NodeRowMap to fluiddis NodeColMap for parallel
  // accessibility.
  // remark: SetState() can not be used here, because it is designed for dof-based vectors only.
  //------------------------------------------------------------------------------------------------
  const Teuchos::RCP<Epetra_Vector> phinpcol = rcp(new Epetra_Vector(*fluiddis_->NodeColMap()));
  LINALG::Export(*phinprow,*phinpcol);
  
  // Jetzt hab ich einen Vektor auf der ColMap. Was mach ich jetzt damit? Muss an Elemente/Dis gehen!
  // eleparams.set("velocity field",tmp);
  // store vector on fluiddis NodeColMap holding G-function values in member variable
  phinp_ = phinpcol;

  // loop over fluid (combustion) column elements
  // remark: loop over row elements would be sufficient, but enrichment is done in column loop
  for (int iele=0; iele<fluiddis_->NumMyColElements(); ++iele)
  {
    // get element from discretization
    const DRT::Element *ele = fluiddis_->lColElement(iele);

#ifdef DEBUG
    if(ele->Type() != DRT::Element::element_combust3)
      // this is not compulsory, but combust3 elements are expected here!
      dserror("unexpected element type: this should be of combust3 type!");
#endif

    // create refinement cell from a fluid element -> cell will have same geometry as element!
    const Teuchos::RCP<COMBUST::RefinementCell> rootcell = rcp(new COMBUST::RefinementCell(ele));

    // refinement strategy is turned on
    if (Teuchos::getIntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REFINEMENT") == true)
      RefineFlameFront(rootcell);

    else // refinement strategy is turned off
      FindFlameFront(phinp_,rootcell);

    /* jetzt habe ich für jedes Element eine "rootcell", an der entweder (refinement on) ein ganzer
     * Baum von Zellen mit Interfaceinformationen hängt, oder (refinement off) eine einzige Zelle
     * hängt. Aus dieser Information muss jetzt das Oberfläche der Flammenfront berechnet werden. */

    // generate interface (flame front) surface
    CaptureFlameFront(rootcell);
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*
 * class section: refinement                                                                      *
 *------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------*
 | refine the region around the flame front                                           henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::RefineFlameFront(const Teuchos::RCP<const COMBUST::RefinementCell> cell)
{
/*
  hier muss der gesamte rekursive Verfeinerungsprozess gestartet werden, d.h. dass hier der Suchbaum
  gestartet wird bzw. ein Ersatzalgorithmus gerufen wird, bis der Suchbaum läuft.

  rufe FindFlameFront() für eine Verfeinerungszelle
  falls Rückgabe = true (Zelle wird geschnitten)
    teile die Verfeinerungszelle
    -> RefineCell() (input: Zelle; output: geteilte Zellen)
    für alle neuen Verfeinerungszellen
      falls maximale Anzahl von Verfeinerungen noch nicht erreicht
        rekursiver Aufruf: RefineFlamefront() (RefineCell())
*/
//  FindFlameFront(refinementcell); // Der phinpcol Vektor muss hier auch noch reingegeben werden!

  dserror("Die Verfeinerung funktioniert noch nicht!");
  return;
}

/*------------------------------------------------------------------------------------------------*
 | split given refinement into 8 refinement cells                                     henke 12/08 |
 *------------------------------------------------------------------------------------------------*/
/*std::vector<Teuchos::RCP<COMBUST::RefinementCell> > SplitRefinementCell(Teuchos::RCP<COMBUST::RefinementCell> cell)
{
  dserror("Thou shalt not call this function!");
  finercells = null;
  return finercells;
}*/

/*------------------------------------------------------------------------------------------------*
 | find the flame front within a refinement cell according to G-function field        henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::FindFlameFront(
       const Teuchos::RCP<Epetra_Vector>           gfuncvalues, // vector holding G-function field
       const Teuchos::RCP<COMBUST::RefinementCell> cell)
{
  //URSULA
	/* gfuncvalues muss doch eigentlich nicht mehr übergeben werden
	 * habe doch dafür nun Varible phinp_ 
	 */
  //URSULA
	
  // get the element which this cell belongs to
  const DRT::Element* ele = cell->Ele();

  //------------------------------------------------------------------------------------------------
  // congruent (matching) discretizations (Fluid == G-function)
  //------------------------------------------------------------------------------------------------
  if (true)
  {
    //----------------------------------------------------------------------------------------------
    // (refinement level == 0) also means: refinement turned off (fluid element == cell!)
    //----------------------------------------------------------------------------------------------
    if (cell->RefinementLevel() == 0)
    {
      // remark: vector "lm" is neccessary, because ExtractMyValues() only accepts "vector<int>"
      // arguments, but ele->NodeIds delivers an "int*" argument
      vector<int> lm(ele->NumNode());
      // get vector of node GIDs for this element
      const int* nodeids = ele->NodeIds();
      for (unsigned inode=0; inode < lm.size(); inode++)
        lm[inode] = nodeids[inode];
/*
  weitere Implementierungsmöglichkeiten
  2.Möglichkeit, führt Schleife für NodeIds() von Hand aus und spart so Umschreiben vin nodeids auf lm
      // get pointer to nodes of this element
      const DRT::Node* const* nodes = ele->Nodes();
      // loop all nodes of this element
      for (int inode = 0; inode < ele->NumNode(); inode++)
      {
        // get GID of this node (fluid node!)
        const int nodeid = nodes[inode]->Id();
        // build vector of all node GIDs of this element
        lm[inode] = nodeid;
      }
  3.Möglichkeit, die funktioniert, wenn phi auf der gfuncdis DofColMap vorliegen würde.
        // get vector of all dof GIDs of this element
        const vector<int> lm = gfuncdis_->Dof(ele);
        // create local vector "myphinp"
        vector<double> myphinp(lm.size()); // vector<double> degeneriert hier zum scalaren double!
*/
      //--------------------------------------------------------------------------------------------
      // extract G-function values for all nodes belonging to this fluid element ( == cell!)
      //--------------------------------------------------------------------------------------------
      // create vector "mygfuncvalues" holding G-function values for this element
      vector<double> mygfuncvalues(ele->NumNode(),1000.0);
      // get entries in "gfuncvalues" corresponding to node GIDs "lm" and store them in "mygfuncvalues"
      DRT::UTILS::ExtractMyValues(*gfuncvalues,mygfuncvalues,lm);
      //TEST
//      if (ele->Id()==0)
//      {
//        std::cout<< "Gfunc " << ele->Id() << std::endl;
//        for(std::size_t ig=0; ig<mygfuncvalues.size(); ig++)
//    	    std::cout << mygfuncvalues[ig] << std::endl;
//      }

#ifdef DEBUG
        if (lm.size() != mygfuncvalues.size() )
          dserror("unexpected number of nodes: there should be 8 nodes per element!");
#endif

      //--------------------------------------------------------------------------------------------
      // store G-function values in refinement cell and decide whether this cell is intersected
      //--------------------------------------------------------------------------------------------
      // node numbering in element has to match data structure in cell
      cell->SetGfuncValues(mygfuncvalues);

//URSULA
      // TEST
//        vector<double> phis (8);
//        phis[0]=-1;
//        phis[1]=-1;
//        phis[2]=1;
//        phis[3]=1;
//        phis[4]=1;
//        phis[5]=-1;
//        phis[6]=1;
//        phis[7]=2;
        //
//        phis[0]=0;
//        phis[1]=-1;
//        phis[2]=0;
//        phis[3]=1;
//        phis[4]=0;
//        phis[5]=-1;
//        phis[6]=0;
//        phis[7]=1;
        //
//        phis[0]=-1.5;
//        phis[1]=0.5;
//        phis[2]=-1.5;
//        phis[3]=0.5;
//        phis[4]=-1.5;
//        phis[5]=0.5;
//        phis[6]=-1.5;
//        phis[7]=0.5;
        //
//        phis[0]=0;
//        phis[1]=0;
//        phis[2]=0;
//        phis[3]=0;
//        phis[4]=0;
//        phis[5]=0;
//        phis[6]=0;
//        phis[7]=1;
        //
//       phis[0]=1.5;
//       phis[1]=-0.5;
//       phis[2]=1.5;
//       phis[3]=2;
//       phis[4]=1.5;
//       phis[5]=-0.5;
//       phis[6]=1.5;
//       phis[7]=2;
//        cell->SetGfuncValues(phis);
//URSULA

//      std::cout << "G-Funktionswerte in Zelle verpackt für Element: " << cell->Ele()->Id() << " - geschnitten?: " << cell->Intersected() << endl;
//      for (unsigned i=0; i < mygfuncvalues.size(); i++)
//        std::cout << "Wert am Knoten " << i << ": " << mygfuncvalues[i] << endl;
    }
    //----------------------------------------------------------------------------------------------
    // higher level of refinement
    //----------------------------------------------------------------------------------------------
    else
    {
      dserror("refinement not working yet!");
      /* hole GID dieses Fluid Elementes
           // const int gid = ele->Id();
         für alle Eckpunkte der Zelle
           werte G-Funktion an lokalen Koordinaten des Elementes GID aus
             // GID der Elemente beider Dis identisch!
             // lokale Fluid und G-Funktion Koordinaten identisch!
           speichere den Wert am Eckpunkt */
    }
  }
  //------------------------------------------------------------------------------------------------
  // no congruent (matching) discretizations (Fluid != G-function)
  //------------------------------------------------------------------------------------------------
  else
  {
    dserror("non-congruent discretizations cannot be handled yet!");
    if (true) // (refinement level == 0) also means: refinement turned off (fluid element == cell!)
    {
      dserror("non-congruent discretizations not yet working!");
      /* für jeden Knoten des Fluid Elementes
           hole die globalen Koordinaten dieses Knotens
           finde GID des zugehörigen LevelSet Knotens
             // gleiche globale Koordinaten!
           hole Wert der G-Funktion für diese GID aus Lösungsvektor
           speichere den Wert am Eckpunkt
             // hier muss Knoten mit Eckpunkt richtig verbunden werden! */
    }
    else // (refinement level != 0) higher level of refinement
    {
      dserror("refinement for non-congruent discretizations not yet working!");
      /* für alle Eckpunkte der Zelle
           berechne globale Koordinaten des Eckpunktes
           finde GID des zugehörigen LevelSet Elements
             // Eckpunkt wird oft auf einer Elementkante liegen!
           transformiere die Eckpunktkoord. in lokale LevelSet Koord.
           werte G-Funktion an lokalen Koordinaten des Elementes GID aus
           speichere den Wert am Eckpunkt */
    }
  }

  //------------------------------------------------------------------------------------------------
  // find intersection points of G-function (zero level set) with refinement cell edges
  //------------------------------------------------------------------------------------------------
  if(cell->Intersected() == true)
  {
    FindIntersectionPoints(cell);
  }
  else
  {
    // done (-> next fluid element)
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | find intersection points of G-function (level set zero iso-surface) with refinement cell edges |
 |
 | this function does default bullshit for the moment                                 henke 03/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::FindIntersectionPoints(const Teuchos::RCP<COMBUST::RefinementCell> cell)
{
  // get G-function values from refinement cell
  std::vector<double> gfuncvalues = cell->GetGfuncValues(); //-> an den Ecken
  //TEST
//  std::cout<< "Gfunc " << std::endl;
//  for(std::size_t ig=0; ig<gfuncvalues.size(); ig++)
//	  std::cout << gfuncvalues[ig] << std::endl;
  // TEST mit vorgegebenen phi-Werten
  //std::vector<double> gfuncvalues (8);
//  gfuncvalues[0]=1.5;
//  gfuncvalues[1]=1.5;
//  gfuncvalues[2]=1.5;
//  gfuncvalues[3]=1.5;
//  gfuncvalues[4]=-0.5;
//  gfuncvalues[5]=-0.5;
//  gfuncvalues[6]=-0.5;
//  gfuncvalues[7]=-0.5;
//  gfuncvalues[0]=-1.5;
//  gfuncvalues[1]=0.5;
//  gfuncvalues[2]=-1.5;
//  gfuncvalues[3]=-3;
//  gfuncvalues[4]=-1.5;
//  gfuncvalues[5]=0.5;
//  gfuncvalues[6]=-1.5;
//  gfuncvalues[7]=-3;

  
  // ich gehe im Moment davon aus, dass GetVertexCoord vertexcoordinaten liefert,
  // die im Elementkoordinatensystem des zugehörigen Elements dargestellt sind
  // Nummerierung entsprechend Konvention fuer Elemente -> hier hex8
  std::vector<std::vector<double> > vertexcoord = cell->GetVertexCoord();
  
  std::map<int,std::vector<double> > intersectionpoints;
  
  // mögliche Vorgehensweise:
  // bekannt: phi an den Ecken der Zelle
  //          falls hex8 maximal ein Schnittpunkt je Kante (linearer Verlauf entlang Kante)
  // Schnittpunkt kann direkt durch lin Interpolation oder mit Hilfe
  // der Ansatzfunktionen bestimmt werden
  // außerdem ist es möglicherweise ganz gut zu wissen zu welcher Kante der
  // Schnittpunkt gehört, daher
  // große Schleife über alle Kanten (Reihenfolge entsprechend Festlegung in globalreport)
  // bestimme Schnittpunkt in lokalen (xi) Koordinaten (Element)
  // muss dann eventuell in globale Koordinaten umgerechnet werden?
  // Umrechnung dann aber auch für die Ecken der Zelle
  // ich glaube in lokalen Koordinaten des zugrundeliegenden Elements reicht aus (siehe Scilab Testprogramm)
  // falls Schnittpunkt vorhanden wird dieser zur Map intersectionpoints dazugefügt
  // dabei wird der int-Wert mit der Kantennummer gefüllt und der vector
  // erhält die Koordinaten (global oder lokal)
  
  // Lines and corresponding nodes (hex8)
  /* L1:  0 1
   * L2:  1 2
   * L3:  2 3
   * L4:  0 3
   * L5:  0 4
   * L6:  1 5
   * L7:  2 6
   * L8:  3 7
   * L9:  4 5
   * L10: 5 6
   * L11: 6 7
   * L12: 4 7
   */
  if (cell->Ele()->Shape()!=DRT::Element::hex8)
	  dserror("FindIntersectionPoints not supported for this element shape!");
	  
  std::vector<std::vector<int> > lines = DRT::UTILS::getEleNodeNumberingLines(DRT::Element::hex8);
//  switch(cell->Ele()->Shape())
//  {
//  case DRT::Element::hex8:
//	  {
//		  lines = DRT::UTILS::getEleNodeNumberingLines(DRT::Element::hex8);
//		  break;
//	  }
//  default:
//	  dserror("FindIntersectionPoints not supported for this element shape!");
//  }
  
  // hex8 only
  // hex20, hex27 schwieriger, da an den Kanten quad Verlaeufe
  for(std::size_t ilines=0; ilines<lines.size(); ilines++)
  {
	 double gfuncval1 = gfuncvalues[lines[ilines][0]];
	 double gfuncval2 = gfuncvalues[lines[ilines][1]];
	 
	 std::vector<double> coordinates(3);
	 // test intersection: change of sign
	 if (gfuncval1*gfuncval2<0) //real intersection point
	 {
	      
		  for (int dim = 0; dim < 3; dim++)
		  {
			  if(vertexcoord[lines[ilines][0]][dim]==vertexcoord[lines[ilines][1]][dim])
			  {
				  coordinates[dim] = vertexcoord[lines[ilines][0]][dim];
			  }
			  else // calculate intersectionpoint
			  {
				  // in case of hex8, edge order = 1 -> lin interpolation
				  /*
				   * x = x1 + (phi(=0) - phi1)/(phi2 - phi1)*(x2 - x1)
				   */
				  
				  coordinates[dim] = vertexcoord[lines[ilines][0]][dim] - gfuncval1 / (gfuncval2 - gfuncval1) * (vertexcoord[lines[ilines][1]][dim] - vertexcoord[lines[ilines][0]][dim]);
			  }
		      // now vector "cooordinates" contains intersection point in local element coordinates
		  }

			 /*
			 // umrechnen in globale (physikalische) Koordinaten 
			 // x = Nxi
			 const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
			 static LINALG::Matrix<numnode,1> funct;
			 DRT::UTILS::shape_function_3D(funct,coordinates[0],coordinates[1],coordinates[2],DRT::Element::hex8);
			 //jetzt braucht man noch die Knotenkoordinaten des zugehörigen Elements (global)
			 std::vector<double> global_coordinates(3);
			 for(int dim=0; dim<3; dim++)
			 {
				 global_coordinates[dim] = 0.0;
			 }
			 for(int inode=0; inode<numnode; inode++)
			 {
				 for(int dim=0; dim<3; dim++)
				 {
					 global_coordinates[dim] += funct(inode) * cell->Ele()->Nodes()[inode]->X()[dim];//ele.Nodes()[inode]->X()[dim];//Knotenkoordinaten des zugehörigen Elements (global)[inode][dim];
				 }
			 }
		     */

			 intersectionpoints[ilines] = coordinates;
	 }
	 else
	 {
		 //do nothing and go to the next line
	 }

  }
  
  //TEST Ausgabe
//  for (std::map<int,std::vector<double> >::const_iterator iter = intersectionpoints.begin(); iter != intersectionpoints.end(); ++iter)
//  {
//	  std::cout<< iter->first << std::endl;
//	  std::vector<double> coord = iter->second;
//	  for (std::size_t isd=0; isd<3; isd++)
//	  {
//		  std::cout<< coord[isd] << std::endl;
//	  }
//  }
  
  // store intersection points in refinement cell
  cell->intersectionpoints_ = intersectionpoints;
  return;
}
//URSULA

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*
 * class section: capture of flame front                                                          *
 *------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------*
 | capture flame front within one refinement cell to provide integration cells        henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::CaptureFlameFront(const Teuchos::RCP<const COMBUST::RefinementCell> cell)
{
  /*
     triangulate flame front surface (spanned by the intersection points) in every refinement cell
     create a piecewise linear complex (PLC) in Tetgen format from every refinement cell
     call the Constraint Delaunay Tetrahedrization (CDT) to produce burnt and unburnt integration cells
     transform integration cells from Tetgen format to baci format

     order function calls:

     TriangulateFlameFront()
     buildPLC()
  */
	
  // URSULA 10.6.09
  // hier geht die rootcell rein, also das Element
  // falls das Element geschnitten ist, hängen weitere Refinementcells dran
  // diese können geschnitten sein oder nicht
  // ich brauche außen eine Schleife über alle Refinementcells der rootcell
  // darin kommt
  // - frage cell, ob geschnitten
  //    falls nicht, ist diese cell bereits eine zum Element gehörende Intergrationszelle
  //	und man kann zur nächsten cell gehen
  //	falls geschnitten, muss Interface trianguliert werden um mit TetGen die
  //	Intergrationszellen zu bestimmen
  //    - TriangulateFlameFront()
  //    - buildPLC() mit CreateIntegrationCells() und TransformIntegrationCells()
  //    - besser ist:
  //      rufe buildPLC() und darin dann erst TriangulateFlameFront()
  // - letztendlich müssen in der Sysmat über das InterfaceHandle Integrationszellen ankommen

  // das Folgende ist bisher nur vorläufig und funktioniert daher nur, falls noch kein
  // Refinement vorgenommen wurde
	
//  // das sind Vektoren mit Integrationzellen bzw Integrationsflächen
//  GEO::DomainIntCells listDomainIntCellsperEle;
//  GEO::BoundaryIntCells listBoundaryIntCellsperEle;
  if (cell->Intersected())
  {	  
	  // das sind Vektoren mit Integrationzellen bzw Integrationsflächen
	  GEO::DomainIntCells listDomainIntCellsperEle;
	  GEO::BoundaryIntCells listBoundaryIntCellsperEle;
	  buildPLC(cell, listDomainIntCellsperEle, listBoundaryIntCellsperEle);
	  elementintcells_[cell->Ele()->Id()] = listDomainIntCellsperEle;
	  if(listBoundaryIntCellsperEle.size()>0)
	    boundaryintcells_[cell->Ele()->Id()] = listBoundaryIntCellsperEle;
	  
	  std::cout << cell->Ele()->Id() << "DomainIntCell " << listDomainIntCellsperEle.size() << std::endl;
	  std::cout << cell->Ele()->Id() << "BoundaryIntCell " << listBoundaryIntCellsperEle.size() << std::endl;
  }
  //falls nicht geschnitten
  // - cell == ele: -> fertig, da GetDomainIntCell diesen Fall berücksichtigt
  //   aber GetDomainIntCell weiß nicht, wo das Element liegt, + oder -
  //   also doch in GetDomainIntCell schreiben 
  //   macht jetzt InterfaceHandle
  // - cell == refinement of ele: -> celle ist DomainIntCell von ele
  else
  {
//	  std::vector<std::vector<double> > vertexcoord = cell->GetVertexCoord();
//	  std::vector<double> gfuncvalue = cell->GetGfuncValues();
//	  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
//	  LINALG::SerialDenseMatrix  xyze(3,numnode);
//	  for(int inode=0;inode<numnode;inode++)
//	  {
//	    xyze(0,inode) = cell->Ele()->Nodes()[inode]->X()[0];//ele.Nodes()[inode]->X()[0];
//	    xyze(1,inode) = cell->Ele()->Nodes()[inode]->X()[1];//ele.Nodes()[inode]->X()[1];
//	    xyze(2,inode) = cell->Ele()->Nodes()[inode]->X()[2];//ele.Nodes()[inode]->X()[2];
//	  }
//	  //Zelle ist nicht geschnitten-> Zelle ist Integrationszelle
//	  std::cout<< "füge DomIntCell hinzu" <<std::endl;
//	  //entspricht vertexcoord
//	  LINALG::SerialDenseMatrix CellCoord(3, numnode);
//	  //falls Element, dann ist das xyze
//	  LINALG::SerialDenseMatrix physCellCoord(3, numnode);
//	  for(int inode=0; inode<numnode; inode++)
//	  {
//		  static LINALG::Matrix<3,1> ccoord;
//		  for(int dim=0; dim<3; dim++)
//		  {
//			  CellCoord(dim,inode) = vertexcoord[inode][dim];
//			  ccoord(dim) = CellCoord(dim,inode);
//		  }
//	      GEO::elementToCurrentCoordinatesInPlace(DRT::Element::hex8, xyze, ccoord);  
//	      for(int  dim=0; dim<3; dim++)
//	        physCellCoord(dim,inode) = ccoord(dim);
//	  }
//	  bool inGplus = GetIntCellDomain(CellCoord, gfuncvalue, DRT::Element::hex8, DRT::Element::hex8);
//	  //TEST
////	  std::cout << "physCellCoord " << physCellCoord(0,3) << physCellCoord(1,3) << std::endl;
////	    if(inGplus)
////	    {
////	        std::cout << "In G plus" << std::endl;
////	    }
////	    else
////	    {
////	    	std::cout << "In G minus" << std::endl;
////	    }
//	  listDomainIntCellsperEle.push_back(GEO::DomainIntCell(DRT::Element::hex8, CellCoord, physCellCoord, inGplus));
  }
//  std::cout << cell->Ele()->Id() << "DomainIntCell " << listDomainIntCellsperEle.size() << std::endl;
  //TEST
//  GEO::DomainIntCell intcelltest=listDomainIntCellsperEle[0];
//  //std::cout << intcelltest.getLabel() << std::endl;
//  if (intcelltest.getDomainPlus())
//   std::cout << "domain +" << std::endl;
//  else
//   std::cout << "domain -" << std::endl;
  
  //elementintcells_[cell->Ele()->Id()] = listDomainIntCellsperEle;
  //dserror("Ende Schnitt");
  // URSULA

  // hier muss jetzt die map "flamefrontpatches_" mit Oberflächenstückchen gefüllt werden

  // eventuell macht es keinen Sinn die übergebene Zelle "cell" als doppelt "const" zu deklarieren

  flamefrontpatches_.insert(make_pair(cell->Ele()->Id(),cell));

  return;
}

//URSULA
/*------------------------------------------------------------------------------------------------*
 | triangulate the intersection surface of a refinement cell                          henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::TriangulateFlameFront(
		std::vector<std::vector<int> >&        trianglelist,
		std::multimap<int,std::vector<int> >&  segmentlist,
		std::vector<std::vector<double> >&     pointlist,
		std::map<int,int>&          intersectionpoints_id, //
		const std::vector<double>&  gfuncvalue //
		) 
{
  //std::cout << "Triangulierung" << std::endl;
  //check -> is cell really intersected?
 if(segmentlist.count(-1)<segmentlist.size()) //oder !=, segmentlist darf nicht nur Würfelkanten enthalten (key=-1)
 {
  std::map<int,std::vector<int> > polygonpoints;
  int actpoint;

  //Suche Einstieg, falls eine Fläche zwei Interfacestuecke besitzt, beginne ich dort
  //warum: -damit ich das zweite Polygon nicht übersehe
  //       -damit ich nicht mit Kante beginne, da es möglich ist, das diese die Interfacefläche nur berührt
  int numpolygon = 1;
  int beginpolygon=-1;
  for (std::multimap<int,std::vector<int> >::iterator seglistiter = segmentlist.begin(); seglistiter != segmentlist.end(); seglistiter++)
  {
	  if(seglistiter->first!=-1)
	  {
		  if (segmentlist.count(seglistiter->first)==2)
		  {
			  if (numpolygon==1)
			  {
				  numpolygon = 2;
				  beginpolygon = seglistiter->first;
			  }
		  }
		  else if (segmentlist.count(seglistiter->first)==1)
		  {
			  if (beginpolygon<0)
			  {
				  beginpolygon = seglistiter->first;
			  }
		  }
		  else
			  dserror("can't build polygon");
	  }
  }

  //hier werden die Polygone gebildet
  //also die zugehoerigen Punkte in die richtige Reihenfolge gebracht
  //einschließlich des Umlaufsinns
  int j=0;
  //loop over all polygons
  for (std::multimap<int,std::vector<int> >::iterator segmentiter=segmentlist.equal_range(beginpolygon).first; segmentiter!=segmentlist.equal_range(beginpolygon).second; segmentiter++) 
  {
//	  int actkey = segmentiter->first; 
//    std::vector<int> actsegmentpoints = segmentiter->second;
	  // es reicht aus dem ersten Segment den richtigen Umlaufsinn zugeben
	  // der Rest ergibt sich im Folgenden automatisch
      IdentifyPolygonOrientation(segmentiter->second, segmentiter->first, intersectionpoints_id, gfuncvalue); //
//	  IdentifyPolygonOrientation(actsegmentpoints, actkey, intersectionpoints_id, gfuncvalue);
      std::vector<int> actsegmentpoints = segmentiter->second;
      //TEST
//      std::cout << "Startsegment " << segmentiter->first << std::endl;
//      std::cout << actsegmentpoints[0] << std::endl;
//      std::cout << actsegmentpoints[1] << std::endl;

      //contains the vertices of the polygon
	  std::vector<int> polypoints;
	  //store points of start segment
      polypoints.push_back(actsegmentpoints[0]);
      actpoint = actsegmentpoints[1];
      polypoints.push_back(actpoint);

      // jetzt wird Segment gesucht, das daran anschließt
      // also als Endpunkt auch den gerade aktuellen Punkt besitzt
      // dann wird für dieses Segment der Anschluss gesucht
      // das geht solange bis man wieder beim ersten Punkt von polypoints ankommt
	  while(actpoint!=polypoints[0])
      {
	      //muss abgebrochen werden sobald naechster Punkt gefunden ist
          for (std::multimap<int,std::vector<int> >::const_iterator iter = segmentlist.begin(); iter != segmentlist.end(); ++iter)
          {
			  std::vector<int> it_points = iter->second;
    	      if (it_points!=actsegmentpoints) //damit man nicht mit dem im letzten Schritt gefunden Segemnt vergleicht
    	      {
				  if (it_points[0]==actpoint)
				  {
					  actpoint = it_points[1];
					  polypoints.push_back(actpoint);
					  actsegmentpoints = it_points;
					  break; //finished: find next segment
				  }
				  else if (it_points[1]==actpoint)
				  {
					  actpoint = it_points[0];
					  polypoints.push_back(actpoint);
					  actsegmentpoints = it_points;
					  break; //finished: find next segment
				  }
				  else
				  {
				  }
    	    	  
			  }
			  else
			  {
			  }
		  }
	  }
	  //store polypoints in a map containing all polygons
	  polygonpoints[j] = polypoints;
	  j++;
	  //TEST
	    std::cout<<"polygonecken"<< std::endl;
	    for (std::size_t ipoly=0; ipoly<polypoints.size(); ipoly++)
	    {
	  	std::cout<< polypoints[ipoly] << std::endl;
	    }
  }
    
  //jetzt müssen Dreiecke gebildet werden
  //falls polygon nur 3 Punkte verschiedene Punkte enthält ist man fertig
  //bei mehr als 3 Punkten muss "Mittelpunkt" bestimmt werden
  //dieser bildet zusammen mit zwei aufeinanderfolgenden Punkten von polypoints
  //ein Dreieck, das dann bereits den richtigen Umlaufsinn besitzt, so dass Normalenvektor
  // von + nach - zeigt
  for (std::size_t ipolygons=0; ipolygons<polygonpoints.size(); ipolygons++)
  {
	  std::vector<int> polypoints = polygonpoints[ipolygons];
      if (polypoints.size()<4)
	      dserror("TriangulateFlameFront needs at least 3 intersectionpoints");  
  
      if (polypoints.size()==4) //erster und letzter Punkt identisch
      {
	      std::vector<int> trianglepoints (3);
	      for (int i=0; i<3; i++)
	      {
		      trianglepoints[i] = polypoints[i];
	      }
	      trianglelist.push_back(trianglepoints);
      }
      else
      {
	  
	      //calculate midpoint of interface first
	      std::size_t numpoints = polypoints.size() - 1;
	      std::vector<double> midpoint (3);
	      std::vector<double> point1 (3);
	      std::vector<double> point2 (3);
	      //even
	      if (numpoints%2==0)
	      {
		      point1 = pointlist[polypoints[0]];
		      point2 = pointlist[polypoints[numpoints/2]];
	      }
	      //odd
	      else
	      {
		      point1 = pointlist[polypoints[0]];
		      point2 = pointlist[polypoints[(numpoints+1)/2]]; 
	      }
	      for (int dim=0; dim<3; dim++)
	      {
		      midpoint[dim] = (point2[dim] + point1[dim]) * 0.5;
	      } 
	      std::size_t midpoint_id = pointlist.size();//id beginnen bei 0
	      //std::cout<< "pointlistsize " << pointlist.size() << "midpointid " << midpoint_id<<std::endl;
		  pointlist.push_back(midpoint);
		  
          //build triangles
	      for (std::size_t j=0; j<polypoints.size()-1; j++)
	      {
		      std::vector<int> trianglepoints (3);
              trianglepoints[0] = polypoints[j];
              trianglepoints[1] = polypoints[j+1];
              trianglepoints[2] = midpoint_id;
		      trianglelist.push_back(trianglepoints); 
	      }
	  
      }
  }

  //TEST
  std::cout<<"dreiecke"<< std::endl;
  for (std::size_t itriangle=0; itriangle<trianglelist.size(); itriangle++)
  {
	  std::cout<< "dreieck " << itriangle<< std::endl;
	  for (int i=0; i<3; i++)
		  std::cout<< trianglelist[itriangle][i]<<std::endl;
  }
  
 }
  
  return;
}
//URSULA

//URSULA
/*------------------------------------------------------------------------------------------------*
 | indentify the orientation of the interface polygon, normal vector + -> -                       |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::IdentifyPolygonOrientation(
		std::vector<int>&           segment,
		const int                   surf_id,
		std::map<int,int>&          intersectionpoints_id,
		const std::vector<double>&  gfuncvalue
		)
{
	// array containing the line_id's for each surface
	int surfacelines[6][4] = {{3, 2, 1, 0},
			                  {0, 5, 8, 4},
			                  {1, 6, 9, 5},
			                  {2, 7,10, 6},
			                  {4,11, 7, 3},
			                  {8, 9,10,11}};
	// array containing test vertex for each surface
	// that means : test of sign of g-func 
	int testnode[6][4] = {{0, 3, 2, 1},
					      {0, 1, 5, 4},
					      {1, 2, 6, 5},
					      {2, 3, 7, 6},
					      {0, 4, 7, 3},
					      {4, 5, 6, 7}};
	int line_id = -1;
	int testpoint = -1;
	
	const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
	//segment contains intersection point
	if(segment[0]>numnode-1)
	{
		int	interpoint = segment[0];

		for(std::map<int,int>::const_iterator iterinterpoint=intersectionpoints_id.begin(); iterinterpoint!=intersectionpoints_id.end(); iterinterpoint++)
		{
			if(iterinterpoint->second==interpoint)
			{
				line_id = iterinterpoint->first;
			}
		}
		//get point to test, depends on the line_id
		for(int k=0; k<4; k++)
		{
			if(line_id==surfacelines[surf_id][k])
			{
				testpoint = k;
				break;
			}
		}
	}
	else //Diagonalenschnitt
	{
		for(int k=0; k<4; k++)
		{
			if(testnode[surf_id][k]==segment[0])
			{
				if(k==0)
					testpoint = 3;
				else
					testpoint = k-1;
				break;
			}
		}
	}
	
	if(gfuncvalue[testnode[surf_id][testpoint]]>0) //tausche, (in plus <0)
	{
		int temp = segment[0];
		segment[0] = segment[1];
		segment[1] = temp;
	}
	
	return;
}
//URSULA

//URSULA
/*------------------------------------------------------------------------------------------------*
 | build polygon segments of intersection surface of a refinement cell                            |
 | hex8 only                                    |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::buildFlameFrontSegments(
        std::map<int,int>&                    intersectionpoints_id,
        std::multimap<int,std::vector<int> >& segmentlist,
		const std::vector<double>&            gfuncvalue,
		const std::vector<std::vector<double> >& pointlist
		)
{
	  // array containing the line_id's for each surface
	  int surface[6][4] = {{3, 2, 1, 0},
			               {0, 5, 8, 4},
			               {1, 6, 9, 5},
			               {2, 7,10, 6},
			               {4,11, 7, 3},
			               {8, 9,10,11}};
	
	  //int numsurf = DRT::UTILS::getNumberOfElementSurfaces(DRT::Element::hex8);
	  //das ist hex8 spezifisch
	  for (int i=0; i<6; i++)//loop over all surfaces
	  {
		  // gets end points of the segment
		  std::vector<int> segmentpoints;

		  for (int j=0; j<4; j++)//loop over all lines of the surface
		  {
              // find line
			  std::map<int,int>::const_iterator it_intersectionpoint = intersectionpoints_id.find(surface[i][j]);
			  if (it_intersectionpoint!=intersectionpoints_id.end()) //check: is line intersected?
			  {
				  //get id of the intersectionpoint and store it in segmentpoints vector
				  segmentpoints.push_back(it_intersectionpoint->second);
			  }
		  }
		  
		  //std::cout << "segmentpoints größe" << segmentpoints.size() << std::endl;
		  
		  switch (segmentpoints.size()) {
		  case 0:
		  {
			  // check, if gfunc==0 at the vertices
			  // gets nodes with gfunc==0
			  std::vector<int> zeropoints;
			  std::vector<std::vector<int> > surfacepointlist = DRT::UTILS::getEleNodeNumberingSurfaces(DRT::Element::hex8);
			  for(int k=0; k<4; k++)//loop over nodes
			  {
				  if (gfuncvalue[surfacepointlist[i][k]]==0)
					  zeropoints.push_back(surfacepointlist[i][k]);
			  }
			  //std::cout << "zeropoints größe" << zeropoints.size() << std::endl;
			  switch (zeropoints.size()){
			  case 0: //surface not intersected
			  case 1: //surface not intersected, but one vertex touchs the interface
			  {
				  break;
			  }
			  case 2: //zwei gegenueberliegende Ecken sind Null
				      //   - echter Schnitt, entspricht dann der Diagonalen, falls die anderen Ecken unterschiedliches Vorzeichen haben
				      //   - sonst nur Beruehrung
				      //zwei benachbarte Ecken haben phi==0
				      //   - Kante ist "Tangente" (Teil des Interfaces), aber Zelle selbst wird nicht geschnitten
				      //   - Kante ist interfacebegrenzendes Segement -> man braucht Segment in TriangulateFlameFront
			  {
				  // Schnitt nur bei zwei gegenüberliegenden Nullen möglich
				  if(((zeropoints[0]==surfacepointlist[i][0])and(zeropoints[1]==surfacepointlist[i][2])))
				  {
					  //std::cout<< "Diagonalenschnitt"<<std::endl;
					  //dann braucht man noch unterschiedliches Vorzeichen bei den verbleibenden Knoten	
					  if (gfuncvalue[surfacepointlist[i][1]]*gfuncvalue[surfacepointlist[i][3]]<0)
					  {
						  //store in segmentlist
						  segmentlist.insert(pair<int,std::vector<int> >(i,zeropoints));
					  }
				  }
				  else if (((zeropoints[0]==surfacepointlist[i][1])and(zeropoints[1]==surfacepointlist[i][3])))
				  {
					  //std::cout<< "Diagonalenschnitt"<<std::endl;
					  if (gfuncvalue[surfacepointlist[i][0]]*gfuncvalue[surfacepointlist[i][2]]<0)
					  {
						  segmentlist.insert(pair<int,std::vector<int> >(i,zeropoints));
					  }					  
				  }
				  else
				  {
					  //no intersection
					  //jedoch wird eventuell segment benoetigt um ein vollstaendiges Polygon zu haben
					  //ich ordnen die Punkte in zeropoints nach ihrer ID
					  //damit ich dann vergleichen kann, ob das segment schon in der segmentlist
					  //vorhanden ist, nicht das ich es doppelt habe
					  //Key=-1, da es sich nicht eindeutig einer Fläche zuordnen lässt
					  if (zeropoints[0]>zeropoints[1])
					  {
						  int temp =zeropoints[1];
	                      zeropoints[1] = zeropoints[0];
						  zeropoints[0] = temp;
					  }
					  //schon in segmentlist?
					  bool not_in_segmentlist = true;
					  for (std::multimap<int,std::vector<int> >::iterator iter=segmentlist.equal_range(-1).first; iter!=segmentlist.equal_range(-1).second; iter++)
					  {
						  if (iter->second == zeropoints)
							  not_in_segmentlist = false;
					  }
	                  if(not_in_segmentlist)
					     //store in segmentlist
	                     segmentlist.insert(pair<int,std::vector<int> >(-1,zeropoints));
				  }
				  break;
			  }
			  case 3: //Flaeche ist nicht geschnitten, aber zwei aneinander anstoßende Kanten liegen auf dem Interface
			  {
				  //ich füge die Segmente noch hinzu, key=-1
				  //aber nur falls noch nicht vorhanden
				  //am besten nach Punkt_id ordnen, wie bei case2
				  
				  //erst die zeropoints den segmenten zu ordnen
				  //suche dazu welcher Knoten fehlt
				  if(zeropoints[0]==surfacepointlist[i][0])
				  {
					  if(zeropoints[1]==surfacepointlist[i][1])
					  {
						 if(zeropoints[2]==surfacepointlist[i][2])
						 {
							 //der vierte Knoten fehlt -> kein Problem
						 }
						 else
						 {
							  //der dritte Knoten fehlt -> zeropoints neu ordnen
							  int temp = zeropoints[0];
							  zeropoints[0] = zeropoints[2];
							  zeropoints[2] = zeropoints[1];
							  zeropoints[1] = temp;
						 }
					  }
					  else
					  {
						  //der zweite Knoten fehlt -> zeropoints neu ordnen
						  int temp = zeropoints[0];
						  zeropoints[0] = zeropoints[1];
						  zeropoints[1] = zeropoints[2];
						  zeropoints[2] = temp;
					  }
				  }
				  else
				  {
					  //der erste Knoten fehlt -> kein Problem
				  }
				  
				  for (std::size_t k=0; k<zeropoints.size()-1; k++)
				  {
				      std::vector<int> segment (2);
				      segment[0] = zeropoints[k];
				      segment[1] = zeropoints[k+1];
					  if (segment[0]>segment[1])
					  {
						  int temp =segment[1];
						  segment[1] = segment[0];
						  segment[0] = temp;
					  }
				      bool not_in_segmentlist = true;
				      for (std::multimap<int,std::vector<int> >::iterator iter=segmentlist.equal_range(-1).first; iter!=segmentlist.equal_range(-1).second; iter++)
				      {
					       if (iter->second == segment)
						       not_in_segmentlist = false;
				      }
			          if(not_in_segmentlist)
					     segmentlist.insert(pair<int,std::vector<int> >(-1,segment));
				  }
			       
			       
	              break;
			  }
			  case 4: //Seitenflaeche ist Interfacestueck -> Zelle ist nicht geschnitten
				      //aber weitere Seitenflaechen koennen Interfacestuecke sein (ebenfalls 4 Nullwerte)
				      //diese Seitenflaechen sind BoundaryIntCells und muessen da hinzugefuegt werden
			  {
				  //zu beachten ist noch, dass an der Interfaceseitenflache immer 2 Zellen zusammen stoßen
				  //die BoundaryIntCell darf man aber nur einmal haben
				  //Auswahlkriterium: Normalenvektor zeigt von - nach +
				  //dserror("Seite gleich Interface geht noch nicht");

				  //mit 2 und 3 findet man bereits alle Segmente die das Interface begrenzen
				  //und man muss hier nichts mehr tun

				  break;
			  }
			  default:
				  dserror("impossible number of zero values");
			  }
			  break; 
		  }
		  case 1:
		  {
			  //std::cout << "ein Schnittpunkt" << std::endl;
			  //find zeropoint to build segment
			  std::vector<int> zeropoints;
			  std::vector<std::vector<int> > surfacepointlist = DRT::UTILS::getEleNodeNumberingSurfaces(DRT::Element::hex8);
			  for(int k=0; k<4; k++)//loop over nodes
			  {
				  if (gfuncvalue[surfacepointlist[i][k]]==0)
					  zeropoints.push_back(surfacepointlist[i][k]);
			  }
			  if (zeropoints.size()!=1)
				  dserror("can't build intersection segment");
			  
			  std::vector<int> segment (2);
			  segment[0] = segmentpoints[0];
			  segment[1] = zeropoints[0];
			  segmentlist.insert(pair<int,std::vector<int> >(i,segment));
			  
			  break;
		  }
		  case 2:
		  {		  
			  //store segment in segmentlist
			  segmentlist.insert(pair<int,std::vector<int> >(i,segmentpoints));
			  break;
		  }
		  case 3:
		  {
			  dserror("impossible number of intersectionpoints for hex8 element surface");
			  break;
		  }
		  case 4:
		  {
			  // zwei Segmente pro Fläche sind bei hex8-Elementen theoretisch möglich
			  // das entspricht vier Intersectionpoints
			  
			  //suche Segment mit maximaler Länege
			  //das ist nicht mögliche Kombination der Endpunkte
			  //und damit kenne ich den Rest
			  //std::cout << "4 Schnittpunkte " << std::endl;
			  double distance;
			  double maxdist = 0.0;
			  int maxdistcounter = 0;
			  for(int k=0; k<4; k++)//loop over segmentpoints
			  {
			      std::vector<double> point1 = pointlist[segmentpoints[k]];
			      std::vector<double> point2 (3);
				  if(k<3)
				  {
					  std::vector<double> point2 = pointlist[segmentpoints[k+1]]; 
				  }
				  else
				  {
					  std::vector<double> point2 = pointlist[segmentpoints[0]];
				  }
				  distance = sqrt((point1[0]-point2[0])*(point1[0]-point2[0])+(point1[1]-point2[1])*(point1[1]-point2[1])+(point1[2]-point2[2])*(point1[2]-point2[2]));
				  if (distance>maxdist)
					  maxdistcounter = k;
			  }
			  
			  //bilde Segmente
			  //std::cout << "Maxdist" << maxdistcounter << std::endl;
			  if (maxdistcounter==0 or maxdistcounter==2)
			  {
				  std::vector<int> segment1 (2);
				  segment1[0] = segmentpoints[3];
				  segment1[1] = segmentpoints[0];
				  std::vector<int> segment2 (2);
				  segment2[0] = segmentpoints[1];
				  segment2[1] = segmentpoints[2];
				  segmentlist.insert(pair<int,std::vector<int> >(i,segment1));
				  segmentlist.insert(pair<int,std::vector<int> >(i,segment2));
			  }
			  else
			  {
				  std::vector<int> segment1 (2);
				  segment1[0] = segmentpoints[0];
				  segment1[1] = segmentpoints[1];
				  std::vector<int> segment2 (2);
				  segment2[0] = segmentpoints[2];
				  segment2[1] = segmentpoints[3];
				  segmentlist.insert(pair<int,std::vector<int> >(i,segment1));
				  segmentlist.insert(pair<int,std::vector<int> >(i,segment2));
			  }
			  break;
		  }
		  default:
		      dserror("unexpected number of intersectionpoints!");
	      }
	  }
	
	  return;
}
//URSULA

//URSULA
/*------------------------------------------------------------------------------------------------*
 | build piecewise linear complex (PLC) in Tetgen format,based on Ursulas preparePLC() henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::buildPLC(
		const Teuchos::RCP<const COMBUST::RefinementCell> cell,
		GEO::DomainIntCells& domainintcelllist,
		GEO::BoundaryIntCells& boundaryintcelllist)
{
  // übergeben wird RefinementCell
  // für PLC braucht man (siehe auch Ursulas Intersection-Klasse)
  // - pointlist-> Ecken, Schnittpunkte, Mittelpunkt des Interfaces(Triangulierung)
  // - XFEMSurfaces-> Elementflächen bzw Seitenfläche der Refinementcell
  // - trianglelist-> Dreiecke, die Interface darstellen, werden in TriangulateFlameFront() gebaut
  // - segmentlist -> Schnittkurve Interface mit Elementflächen, werden in buildFlameFrontSegments() gebaut
  // das muss schließlich an CreateIntegrationCells() übergeben werden
	
  std::vector<std::vector<double> > pointlist;
  std::multimap<int,std::vector<int> > segmentlist;
  std::vector<std::vector<int> > trianglelist;
  //sollte hier eigentlich nicht mehr nötig sein
  //da ich nur die Nummerierung der Ecken möchte
  //ich lasse es trotzdem mal drin
  if (cell->Ele()->Shape()!=DRT::Element::hex8)
	  dserror("FindIntersectionPoints not supported for this element shape!");
  //für jede Seitenfläche zugehörige Knoten-ID
  //diese wird dann als Punkt-ID verwendet
  std::vector<std::vector<int> > XFEMsurfacelist = DRT::UTILS::getEleNodeNumberingSurfaces(DRT::Element::hex8);
  //entspricht intersectionpoints, enthält jedoch anstatt der Punktkoordinaten Punkt-ID
  std::map<int,int> intersectionpoints_id;
   
  std::vector<std::vector<double> > vertexcoord = cell->GetVertexCoord();
  std::map<int,std::vector<double> > intersectionpoints = cell->intersectionpoints_;
  std::vector<double> gfuncvalue = cell->GetGfuncValues();

  //brauche ich später für die Integrationszellen
  //globalen Koordinaten des zugehörigen Elements
  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
  LINALG::SerialDenseMatrix  xyze(3,numnode);
  for(int inode=0;inode<numnode;inode++)
  {
    xyze(0,inode) = cell->Ele()->Nodes()[inode]->X()[0];//ele.Nodes()[inode]->X()[0];
    xyze(1,inode) = cell->Ele()->Nodes()[inode]->X()[1];//ele.Nodes()[inode]->X()[1];
    xyze(2,inode) = cell->Ele()->Nodes()[inode]->X()[2];//ele.Nodes()[inode]->X()[2];
  }
  //std::cout << "xyze " << xyze(0,3) << xyze(1,3) << std::endl;
  
  // jetzt fülle ich die pointlist
  // das ist allgemein gültig, egal welcher Elementtyp
  // als erstes die Ecken des Würfels
  int numofpoints = 0;
  int numvertex = DRT::UTILS::getNumberOfElementCornerNodes(DRT::Element::hex8);  	
  for (int ivertex=0; ivertex<numvertex; ivertex++) //Hex hat 8 Ecken, weitere Koord in vertexcoord bei hex20-27 sind innere Knoten und keine Ecken
  {
	  pointlist.push_back(vertexcoord[ivertex]);
	  numofpoints++;
  }
  // und dann die Intersectionpoints
  // gleichzeitig wird intersectionpoints_id mit den IDs gefüllt
  for (std::map<int,std::vector<double> >::const_iterator iter = intersectionpoints.begin(); iter != intersectionpoints.end(); ++iter)
  {
	  pointlist.push_back(iter->second);
	  intersectionpoints_id[iter->first] = numofpoints;
	  numofpoints++;
  }
  //TEST Ausgabe
//    for (std::size_t iter=0; iter<pointlist.size(); ++iter)
//    {
//  	  std::cout<< iter << std::endl;
//  	  std::vector<double> coord = pointlist[iter];
//  	  for (std::size_t isd=0; isd<3; isd++)
//  	  {
//  		  std::cout<< coord[isd] << std::endl;
//  	  }
//    }
  
  //nun bilde ich die Segmente, die die Interfacefläche im Element begrenzen
  //diese werden dann an TriangulateFlameFront() übergeben,
  //so das darin Dreiecke gebildet werden
  //das geschieht durch Aufruf einer extra Funktion, da einige Sonderfaelle unterschieden
  //werden muessen
  //das ist hex8 spezifisch
  buildFlameFrontSegments(intersectionpoints_id, segmentlist, gfuncvalue, pointlist);
  // TEST
  std::cout<<"Segmente"<< std::endl;
  for (std::map<int,std::vector<int> >::const_iterator iter = segmentlist.begin(); iter != segmentlist.end(); ++iter)
  {
	std::cout<<"Segment "<< iter->first << std::endl;
  	std::vector<int> point = iter->second;
  	for (std::size_t isd=0; isd<2; isd++)
  	{
  		std::cout<< point[isd] << std::endl;
  	}
  }
  
  TriangulateFlameFront(trianglelist, segmentlist, pointlist, intersectionpoints_id, gfuncvalue);  
  //TEST
  for (std::size_t iter=0; iter<pointlist.size(); ++iter)
  {
	  std::cout<< iter << std::endl;
	  std::vector<double> coord = pointlist[iter];
	  for (std::size_t isd=0; isd<3; isd++)
	  {
		  std::cout<< coord[isd] << std::endl;
	  }
  }

  // falls nun Dreiecke vorhanden sind, dann ist das Element irgendwie geschnitten und
  // die Bildung von Tetraeder zur Integration ist notwendig
  // obwohl man aus Tetgen auch Dreiecke als Interfacestücke bekommt, schreibe ich diese
  // nicht heraus, da der Umlaufsinn nicht eindeutig ist
  // stattdessen werden die Dreicke der trianglelist direkt verwendet
  if(trianglelist.size()!=0)
  {
	//determine integration cells
    CreateIntegrationCells(pointlist, segmentlist, XFEMsurfacelist, trianglelist, domainintcelllist, xyze);
    //store boundary integration cells in boundaryintcelllist
    for (std::size_t itriangle=0; itriangle<trianglelist.size(); itriangle++)
    {
    	LINALG::SerialDenseMatrix TriangleCoord(3,3); //3Richtungen, 3Knoten
    	LINALG::SerialDenseMatrix physTriangleCoord(3,3);
    	for (int inode=0; inode<3; inode++)
    	{
    		static LINALG::Matrix<3,1> tcoord;
    		for (int dim=0; dim<3; dim++)
    		{
    			TriangleCoord(dim,inode) = pointlist[trianglelist[itriangle][inode]][dim];
    			tcoord(dim) = TriangleCoord(dim,inode);
    		}
  	      GEO::elementToCurrentCoordinatesInPlace(DRT::Element::hex8, xyze, tcoord);  
  	      for(int  dim=0; dim<3; dim++)
  	        physTriangleCoord(dim,inode) = tcoord(dim);
    	}
    	boundaryintcelllist.push_back(GEO::BoundaryIntCell(DRT::Element::tri3, -1, TriangleCoord, Teuchos::null, physTriangleCoord));
    }
  }
  else
  {
	  //Zelle ist nicht geschnitten-> Zelle ist Intergrationszelle
	  std::cout<< "Zelle nicht geschnitten" <<std::endl;
	  std::cout<< "füge DomIntCell hinzu" <<std::endl;
	  // was jetzt kommt berauche ich für DomainIntCell
	  //entspricht vertexcoord
	  LINALG::SerialDenseMatrix CellCoord(3, numnode);
	  //falls Element, dann ist das xyze
	  LINALG::SerialDenseMatrix physCellCoord(3, numnode);
	  for(int inode=0; inode<numnode; inode++)
	  {
		  static LINALG::Matrix<3,1> ccoord;
		  for(int dim=0; dim<3; dim++)
		  {
			  CellCoord(dim,inode) = vertexcoord[inode][dim];
			  ccoord(dim) = CellCoord(dim,inode);
		  }
	      GEO::elementToCurrentCoordinatesInPlace(DRT::Element::hex8, xyze, ccoord);  
	      for(int  dim=0; dim<3; dim++)
	        physCellCoord(dim,inode) = ccoord(dim);
	  }
	  
	  //das macht jetzt das InterfaceHandle
	  // es reicht daher 
	  domainintcelllist.push_back(GEO::DomainIntCell(DRT::Element::hex8, CellCoord, physCellCoord));
//	  // nun noch das Fluidgebiet bestimmen, in dem die Zelle liegt
//	  bool inGplus = GetIntCellDomain(CellCoord, gfuncvalue, DRT::Element::hex8, DRT::Element::hex8);
//	  //TEST
////	  std::cout << "physCellCoord " << physCellCoord(0,3) << physCellCoord(1,3) << std::endl;
//	    if(inGplus)
//	    {
//	        std::cout << "In G plus" << std::endl;
//	    }
//	    else
//	    {
//	    	std::cout << "In G minus" << std::endl;
//	    }
//	  domainintcelllist.push_back(GEO::DomainIntCell(DRT::Element::hex8, CellCoord, physCellCoord, inGplus));
	  
	  if(segmentlist.size()>=4) //ab hier nur hex8
	  {
		  //möglicherweise müssen noch Randintergationszellen berücksichtigt werden
		  //Normalenvektor soll von + nach - zeigen
		  //Seitenfläche ist nur zu berücksichtigen, falls sich zugehörige Zelle in + befindet
		  //macht man das nicht, hat man die Fläche doppelt
		  bool in_fluid_plus = true;
		  std::vector<int> boundarysurfaces;
		  
		  int isurf = 0;
		  int numsurf = DRT::UTILS::getNumberOfElementSurfaces(DRT::Element::hex8);
		  while((isurf<numsurf) and (in_fluid_plus==true))
		  {
			  int numzeronodes = 0;
			  //std::vector<std::vector<int> > surfacepointlist = DRT::UTILS::getEleNodeNumberingSurfaces(DRT::Element::hex8);
			  for(int k=0; k<4; k++)//loop over nodes
			  {
				  if (gfuncvalue[XFEMsurfacelist[isurf][k]]==0)
					  numzeronodes++;
				  if (gfuncvalue[XFEMsurfacelist[isurf][k]]<0)
					  in_fluid_plus = false;
			  }
			  if (numzeronodes==4)
				  boundarysurfaces.push_back(isurf);
			  isurf++;
		  }
		  
		  if(in_fluid_plus)
		  {
			  //store boundary integration cells in boundaryintcelllist
			  // Umlaufsinn entspricht den Knoten
			  std::cout<< "füge BoundIntCell hinzu" <<std::endl;
			  for(std::size_t ibound=0; ibound<boundarysurfaces.size(); ibound++)
			  {
				  //std::cout<< boundarysurfaces[ibound] <<std::endl;
			      LINALG::SerialDenseMatrix QuadCoord(3,4); //3Richtungen, 4Knoten
			      LINALG::SerialDenseMatrix physQuadCoord(3,4);
				  for (int inode=0; inode<4; inode++)
				  {
					  static LINALG::Matrix<3,1> qcoord;
				  	  for (int dim=0; dim<3; dim++)
				      {
				    	  QuadCoord(dim,inode) = pointlist[XFEMsurfacelist[boundarysurfaces[ibound]][inode]][dim];
				    	  qcoord(dim) = QuadCoord(dim,inode);
				      }
				  	  GEO::elementToCurrentCoordinatesInPlace(DRT::Element::hex8, xyze, qcoord);  
				  	  for(int  dim=0; dim<3; dim++)
				  	     physQuadCoord(dim,inode) = qcoord(dim);
				  }
				  boundaryintcelllist.push_back(GEO::BoundaryIntCell(DRT::Element::quad4, -1, QuadCoord, Teuchos::null, physQuadCoord)); 
			  }
		  }
	  }
	  
  }

  return;
}
//URSULA

//URSULA
/*------------------------------------------------------------------------------------------------*
 | calls the CDT to create burnt and unburnt integration cells                        henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::CreateIntegrationCells(
		const std::vector<std::vector<double> >&  pointlist,
		std::multimap<int,std::vector<int> >&     segmentlist,
		const std::vector<std::vector<int> >&     XFEMsurfacelist,
		const std::vector<std::vector<int> >&     trianglelist,
		GEO::DomainIntCells&                      domainintcelllist,
		//GEO::BoundaryIntCells&                    boundaryintcelllist,
		const LINALG::SerialDenseMatrix    xyze
		//const std::vector<double>& gfuncvalue
        ) // output: integration cells in Tetgen format
{
  // URSULA
  // übergeben werden Referenzen auf die in buildPLC genannten Listen
  // damit muss nun die Variable tetgenio in gefüllt werden
  // zusätzlich braucht man noch tetgenio out, wird von TetGen gefüllt
  // damit dann in TransformIntegrationCells()

  const int dim = 3; 
  tetgenio in;
  tetgenio out;
  char switches[] = "pQ";    //- p     tetrahedralizes a PLC 
                             //-Q      no terminal output except errors
  tetgenio::facet *f;
  tetgenio::polygon *p;
  
  //brauche ich auch scalefactor (siehe GEO::Intersection)?
  
  //allocate point list
  in.numberofpoints = pointlist.size();
  in.pointlist = new REAL[in.numberofpoints * dim];

  // fill point list
  int fill = 0;
  for(int i = 0; i <  in.numberofpoints; i++)
  {
      for(int j = 0; j < dim; j++)
      {
          double coord = pointlist[i][j];
          in.pointlist[fill] = (REAL) coord; //(REAL)?
          fill++;
      }
  }
  
  in.numberoffacets = XFEMsurfacelist.size() + trianglelist.size();
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  
  in.facetmarkerlist = new int[in.numberoffacets];
  
  // loop over all element surfaces
  for(std::size_t i=0; i<XFEMsurfacelist.size(); i++)
  {
	  f = &in.facetlist[i];
	  //map
//	  int numsegments = 0;
//	  std::map<int,std::vector<int> >::const_iterator it_segmentlist = segmentlist.find(i);
//	  if (it_segmentlist!=segmentlist.end()) //check: is surface intersected?
//	  {
//	  	 numsegments = 1;
//	  }
//	 // falls zweimal geschnitten : numpolygons = 3
	 //multimap
	  int numsegments = (int) segmentlist.count(i);
	 f->numberofpolygons = 1 + numsegments;
	 f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
	 f->numberofholes = 0;
	 f->holelist = NULL;
	 // face of hex
	 p = &f->polygonlist[0];
	 p->numberofvertices = 4; //Hexaederseitenfläche hat vier Ecken, evtl allgemeiner XFEMsurfacelist[i].size()
	 p->vertexlist = new int[p->numberofvertices];
	 for(int ivertex = 0; ivertex < p->numberofvertices; ivertex++)
	 {
	     p->vertexlist[ivertex] = XFEMsurfacelist[i][ivertex]; //hier braucht man jetzt die richtigen Id's
	 }
	 
	 //store segments if present
	 if (numsegments)
	 {
		 //map
//		 for(int j=1; j<(numsegments+1); j++) //vorsicht bei Obergrenze
//		 {
//			 p = &f->polygonlist[j];
//			 p->numberofvertices = 2;
//			 p->vertexlist = new int[p->numberofvertices];
//			 for(int k=0; k<2; k++)
//			 {
//				 p->vertexlist[k] = segmentlist[i][k]; //surface i, point k
//				 //TEST
////				 std::cout << "Surface " << i << std::endl;
////				 std::cout << p->vertexlist[k] << std::endl;
//			 }
//		 }
		 //multimap
         int j=1;
		 for (std::multimap<int,std::vector<int> >::iterator iter=segmentlist.equal_range(i).first; iter!=segmentlist.equal_range(i).second; iter++)
		 {
			 p = &f->polygonlist[j];
			 p->numberofvertices = 2;
			 p->vertexlist = new int[p->numberofvertices];
			 std::vector<int> segmentpoints = iter->second;
			 for(int k=0; k<2; k++)
			 {
				 p->vertexlist[k] = segmentpoints[k];
			 }
			 j++;
		 }
	 }
	 
	 //in.facetmarkerlist[i] = 0; nur für BoundaryIntCells
  }
  // store triangles
  std::size_t k = XFEMsurfacelist.size();
  for(std::size_t i=0; i<trianglelist.size(); i++)
  {
	    f = &in.facetlist[k];
	    f->numberofpolygons = 1;
	    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
	    f->numberofholes = 0;
	    f->holelist = NULL;
	    p = &f->polygonlist[0];
	    p->numberofvertices = 3; //Dreieck hat 3 Ecken
	    p->vertexlist = new int[p->numberofvertices];
	    for(int j=0; j<3; j++)
	    {
	    	p->vertexlist[j] = trianglelist[i][j];
	    	//TEST
//			std::cout << "Triangle " << i << std::endl;
//			std::cout << k << std::endl;
//			std::cout << p->vertexlist[j] << std::endl;
	    }
	    
	    //in.facetmarkerlist[k] = 1; nur für BoundaryIntCells
	    
	    k++;
  }
  
  std::cout << "-----Start Tetgen-----" << std::endl;
  //in.save_nodes("tetin");
  //in.save_poly("tetin");
  tetrahedralize(switches, &in,  &out);
  //out.save_nodes("tetout");
  //out.save_elements("tetout");
  //out.save_faces("tetout");
  std::cout << "----- End Tetgen -----" << std::endl;
  
  TransformIntegrationCells(out, domainintcelllist, xyze);
  
  // URSULA
  return;
}
//URSULA

//URSULA
/*------------------------------------------------------------------------------------------------*
 | this could easlily be integrated into CreateIntegrationCells()                     henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::TransformIntegrationCells(
		tetgenio&                out,
		GEO::DomainIntCells&     domainintcelllist,
		//GEO::BoundaryIntCells&   boundaryintcelllist,
		const LINALG::SerialDenseMatrix xyze
		//const std::vector<double>& gfuncvalue
		) // output: integration cells in baci format
{
  DRT::Element::DiscretizationType distype = DRT::Element::tet4;
  const int numTetNodes = DRT::UTILS::getNumberOfElementNodes(distype);
  
  for(int i=0; i<out.numberoftetrahedra; i++ )
  {
    LINALG::SerialDenseMatrix tetrahedronCoord(3, numTetNodes);
    LINALG::SerialDenseMatrix physTetrahedronCoord(3, numTetNodes);
    for(int j = 0; j < numTetNodes; j++)
    {
      static LINALG::Matrix<3,1> tetCoord;
      for(int k = 0; k < 3; k++)
      {
        tetrahedronCoord(k,j) = out.pointlist[out.tetrahedronlist[i*out.numberofcorners+j]*3+k];
        // k: drei aufeinanderfolgende Einträge in der pointlist, entspricht den 3 Richtungen
        // *3: Richtungen je Knoten
        // tetraherdronlist[i*out.numberofcorners+j]: KontenID für i-tes Element j-ter Knoten
        tetCoord(k) = tetrahedronCoord(k,j);
        //TEST
        //std::cout << "Tetraeder " << i << "Knoten " << j << "Koord " << k << "ist " << tetCoord(k) << std::endl;
      }
      // compute physical coordinates
      GEO::elementToCurrentCoordinatesInPlace(DRT::Element::hex8, xyze, tetCoord);  
      for(int k = 0; k < 3; k++)
        physTetrahedronCoord(k,j) = tetCoord(k);
    }
    
    // if degenerate don t store
    if(!GEO::checkDegenerateTet(numTetNodes, tetrahedronCoord, physTetrahedronCoord))
    {
    	//Das übernimmt das InterfaceHandle
    	domainintcelllist.push_back(GEO::DomainIntCell(distype, tetrahedronCoord, physTetrahedronCoord));
//       //ich habe einen weiteren Constructor für DomainIntCell gebaut,
//       //da ich label gleich 0 setzen möchte und das ist das wichtige der Intergrationszelle gleich
//       //mitgeben möchte, ob sie in + oder minus liegt
//       //das heißt dann aber auch, dass ich auch nicht geschnittenen Elemente der in CaptureFlameFront
//       //definierten map hinzufügen muss, da nun der Defaultfall nicht mehr zuverwenden ist,
//       //da nichts mehr über das Fluidgebiet bekannt ist
//       bool inGplus = GetIntCellDomain(tetrahedronCoord, gfuncvalue, DRT::Element::hex8, distype);
//       if(inGplus)
//       {
//           std::cout << "In G plus" << std::endl;
//       }
//       else
//       {
//    	   std::cout << "In G minus" << std::endl;
//       }
//       domainintcelllist.push_back(GEO::DomainIntCell(distype, tetrahedronCoord, physTetrahedronCoord, inGplus));
     }
  }

  //das mache ich nicht!!!!!!!
  //dann kommen noch die BoundaryIntCells dran
  //wenn man die BoundaryIntCells aus Tetgen übernehmen möchte, muss facetmarker verwenden
  //Interfaceflachen haben falschen Umlaufsinn 
  //aber eigentlich muss man die Interfacedreiecke nicht rausschrieben, da man die Dreiecke bereits
  //hat und im Moment auch keine Verscheibung der Knoten auf das Interface vorgesehen ist
////  int anzahl = 0;
//  for(int i=0; i<out.numberoftrifaces; i++)
//  {
//	  if(out.trifacemarkerlist[i]==1)
//	  {
//		  LINALG::SerialDenseMatrix triangleCoord(3, 3);
//		  for(int j = 0; j < 3; j++)
//		  {
//		    //static LINALG::Matrix<3,1> triCoord;
//		    for(int k = 0; k < 3; k++)
//		    {
//		    	triangleCoord(k,j) = out.pointlist[out.trifacelist[i*3+j]*3+k];
////		    	std::cout << "Tri " << i << std::endl;
////		    	std::cout << "k " << k << " j " << j << std::endl;
////		    	std::cout << out.trifacelist[i*3+j] << std::endl;
////		    	std::cout << out.pointlist[out.trifacelist[i*3+j]*3+k] << std::endl;
//		    }		  
//		  }
////		  anzahl++;
//	  }
//  }
////  std::cout << anzahl << std::endl;
  
  return;
}
//URSULA

//URSULA
/*------------------------------------------------------------------------------------------------*
 | compute average GfuncValue for integration cell                                                |
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::FlameFront::GetIntCellDomain(
		const LINALG::SerialDenseMatrix         IntCellCoord,
		const std::vector<double>&              gfuncvalue,
		const DRT::Element::DiscretizationType  xfem_distype,
		const DRT::Element::DiscretizationType  cell_distype
		)
{
	//berechnen zunächst phi an jedem Knoten der DomainIntCell
	//anschließen wird der Mittelwert gebildet
	//falls Mittelwert
	// >0 true
	// <0 false
	
	bool inGplus = false;
	
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
		dserror("Discretization Type (IntCell) not supported yet!");
	}

	int numelenodes = 0;
	if (xfem_distype==DRT::Element::hex8)
	{
		numelenodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
	}
	else
	{
		dserror("Discretization Type (Ele) not supported yet!");
	}
			
	std::vector<double> gvalcellnodes (numcellnodes);
	for (int icellnode=0; icellnode<numcellnodes; icellnode++)
		gvalcellnodes[icellnode] = 0.0;
	
	for (int icellnode=0; icellnode<numcellnodes; icellnode++)
	{
		//calculate shape function at IntCell node
		Epetra_SerialDenseVector  funct(numelenodes);
	    DRT::UTILS::shape_function_3D(funct,IntCellCoord(0,icellnode),IntCellCoord(1,icellnode),IntCellCoord(2,icellnode),DRT::Element::hex8);
	    /*
	     *         ___
	     * gval(x)=\   N(x)*gvali
	     *         /
	     *         ___
	     */
		for (int ielenode=0; ielenode<numelenodes; ielenode++)
		{
			gvalcellnodes[icellnode] += funct(ielenode) * gfuncvalue[ielenode];
		}
	}
	
	//calculate average Gfunc value
	double averageGvalue = 0.0;
	for (int icellnode=0; icellnode<numcellnodes; icellnode++)
		averageGvalue += gvalcellnodes[icellnode];
	averageGvalue /= numcellnodes;
	
	//determine DomainIntCell position
	if(averageGvalue>0)
		inGplus = true;
	if(averageGvalue==0)
		dserror("can't determine DomainIntCell position");
	
	return inGplus;
}
//URSULA

#endif // #ifdef CCADISCRET
