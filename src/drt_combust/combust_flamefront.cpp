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
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_parobject.H"
#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


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
    phinp_(Teuchos::null),
    maxRefinementLevel_(0)
{
  //if (fluiddis->Comm().MyPID() == 0)
  //  std::cout << "Constructing FlameFront done" << std::endl;
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
  // TODO define communicator
  //if (Comm().MyPID()==0)
    std::cout << "\n---  processing flame front ... " << std::flush;;

  //diese Maps müssen bei jeden Durchlauf von ProcessFlameFront neu gefüllt werden
  myelementintcells_.clear();
  myboundaryintcells_.clear();

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
  // Export vector phinp from fluiddis NodeRowMap to fluiddis NodeColMap for parallel
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
    {
//      std::cout << "starting refinement ..." << std::endl;
      maxRefinementLevel_ = Teuchos::getIntParameter(combustdyn.sublist("COMBUSTION GFUNCTION"),"REFINEMENTLEVEL");
//      std::cout << "maximal refinement level " << maxRefinementLevel_ << std::endl;
      //maxRefinementLevel_ = 2;
      if (maxRefinementLevel_ < 0)
        dserror("maximal refinement level not defined");

      RefineFlameFront(rootcell);
    }

    else // refinement strategy is turned off
      FindFlameFront(rootcell);

    /* jetzt habe ich für jedes Element eine "rootcell", an der entweder (refinement on) ein ganzer
     * Baum von Zellen mit Interfaceinformationen hängt, oder (refinement off) eine einzige Zelle
     * hängt. Aus dieser Information muss jetzt das Oberfläche der Flammenfront berechnet werden. */

    // generate interface (flame front) surface
    CaptureFlameFront(rootcell);
  }

  std::cout << "done" << std::endl;
}


/*------------------------------------------------------------------------------------------------*
 | refine the region around the flame front                                           henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::RefineFlameFront(const Teuchos::RCP<COMBUST::RefinementCell> cell)
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

  // get G-Function values of refinement cell and determine, if it is intersected
  FindFlameFront(cell);
  // is maximal refinement level already reached?
  if (cell->RefinementLevel() < maxRefinementLevel_) //->No
  {
     if (cell->Intersected()) // cell is intersected ....
     {
        cell->RefineCell(); // ... and will be refined
//        std::cout << "RefineCell done" << std::endl;
        // loop over all new refinement cells and call of RefineFlameFront
        for (int icell = 0; icell < cell->NumOfChildren(); icell++)
    	    RefineFlameFront(cell->GetRefinementCell(icell));
     }
  }

  //dserror("Die Verfeinerung funktioniert noch nicht!");
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
       const Teuchos::RCP<COMBUST::RefinementCell> cell)
{
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
      DRT::UTILS::ExtractMyValues(*phinp_,mygfuncvalues,lm);
      //TEST
//      if (ele->Id()==0)
//      {
//        std::cout<< "Gfunc " << ele->Id() << std::endl;
//        for(std::size_t ig=0; ig<mygfuncvalues.size(); ig++)
//          std::cout << mygfuncvalues[ig] << std::endl;
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
//        phis[0]=1;
//        phis[1]=1;
//        phis[2]=1;
//        phis[3]=1;
//        phis[4]=1;
//        phis[5]=1;
//        phis[6]=1;
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
      //
//       phis[0]=1;
//       phis[1]=-1;
//       phis[2]=-1;
//       phis[3]=1;
//       phis[4]=1;
//       phis[5]=-1;
//       phis[6]=-1;
//       phis[7]=1;
      //
//       phis[0]=1.5;
//       phis[1]=1.5;
//       phis[2]=1.5;
//       phis[3]=1.5;
//       phis[4]=-0.5;
//       phis[5]=0.1;
//       phis[6]=-0.5;
//       phis[7]=0.1;
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
      /* hole GID dieses Fluid Elementes
           // const int gid = ele->Id();
         für alle Eckpunkte der Zelle
           werte G-Funktion an lokalen Koordinaten des Elementes GID aus
             // GID der Elemente beider Dis identisch!
             // lokale Fluid und G-Funktion Koordinaten identisch!
           speichere den Wert am Eckpunkt */

      // get the vertex coordinates of the refinement cell
      // given in the coordinate system of the element == rootcell
      const std::vector<std::vector<double> > vertexcoord = cell->GetVertexCoord();
      // to be filled with the G-Function values of the refinement cell
      std::vector<double> gfunctionvalues(8);
      // G-Function values of the corresponding element
      if (cell->ReturnRootCell()==NULL)
    	  dserror("return NULL");
      const std::vector<double> gfuncelement = cell->ReturnRootCell()->GetGfuncValues();
      //TEST
//      std::cout<< "G-Funct of rootcell "<< std::endl;
//      for (std::size_t i=0; i<gfuncelement.size(); i++)
//          std::cout<< gfuncelement[i] << std::endl;

      const int numnode = cell->Ele()->NumNode();
//      if (gfuncelement.size() != numnode)
//        dserror("number of G-Function values does not match number of element nodes");

      // loop over all vertices of the refinement cell
      for (std::size_t ivertex = 0; ivertex < vertexcoord.size(); ivertex++)
      {
        // evaluate shape function at the coordinates of the current refinement cell vertex
        Epetra_SerialDenseVector funct(numnode);
        DRT::UTILS::shape_function_3D(funct,vertexcoord[ivertex][0],vertexcoord[ivertex][1],vertexcoord[ivertex][2],DRT::Element::hex8);
        // compute G-Function value
        gfunctionvalues[ivertex] = 0;
        for (int inode = 0; inode < numnode; inode++)
          gfunctionvalues[ivertex] += gfuncelement[inode] * funct(inode);
      }
//      std::cout << "FindFlameFront finished " << std::endl;
      cell->SetGfuncValues(gfunctionvalues);
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
 |                                                                                rasthofer 06/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::FindIntersectionPoints(const Teuchos::RCP<COMBUST::RefinementCell> cell)
{
  // get G-function values at vertices from refinement cell
  const std::vector<double>& gfuncvalues = cell->GetGfuncValues();
  // get vertex coordinates (local fluid element coordinates) from refinement cell
  const std::vector<std::vector<double> >& vertexcoord = cell->GetVertexCoord();
  // temporary variable to store intersection points
  std::map<int,std::vector<double> > intersectionpoints;

  //-------------------------------------------------
  // get vector of lines with corresponding vertices
  //-------------------------------------------------
  // lines: edge numbers and corresponding vertices (hex8)
  std::vector<std::vector<int> > lines;

  switch(cell->Ele()->Shape())
  {
  case DRT::Element::hex8:
  {
    lines = DRT::UTILS::getEleNodeNumberingLines(DRT::Element::hex8);
    // remark: vertices are assumed to be numbered in the same way the nodes are
    //         convention documented in globalreport.pdf
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
    break;
  }
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  {
    // more difficult, since we have quadratic functions along the edges -> check derivatives, too!
    dserror("FindIntersectionPoints() does not support quadatic elements, yet!");
  }
  default:
    dserror("FindIntersectionPoints() does not support this element shape!");
  }

  //-------------------------------
  // determine intersection points
  //-------------------------------
  // loop edges of refinement cell
  for(std::size_t iline=0; iline<lines.size(); iline++)
  {
    // get G-function value of the two vertices defining an edge
    double gfuncval1 = gfuncvalues[lines[iline][0]];
    double gfuncval2 = gfuncvalues[lines[iline][1]];

    std::vector<double> coordinates(3);

    // check for change of sign along edge
    if (gfuncval1*gfuncval2 < 0.0)
    {
      for (int dim = 0; dim < 3; dim++)
      {
        // vertices have one coordinate in common
        if(vertexcoord[lines[iline][0]][dim] == vertexcoord[lines[iline][1]][dim])
        {
          // intersection point has the same coordinate for that direction
          coordinates[dim] = vertexcoord[lines[iline][0]][dim];
        }
        else // compute intersection point
        {
          // linear interpolation (for hex8)
          // x = x1 + (phi(=0) - phi1)/(phi2 - phi1)*(x2 - x1)
          // store intersection point coordinate (local element coordinates) for every component dim
          coordinates[dim] = vertexcoord[lines[iline][0]][dim] - gfuncval1 / (gfuncval2 - gfuncval1)
          * (vertexcoord[lines[iline][1]][dim] - vertexcoord[lines[iline][0]][dim]);

          // shift intersection point to vertex if it is very close
          if(fabs(vertexcoord[lines[iline][0]][dim]-coordinates[dim]) < 1.0E-8)
          {
            coordinates[dim] = vertexcoord[lines[iline][0]][dim];
            cout << "coordinates shifted to vertex 1" << endl;
          }
          if(fabs(vertexcoord[lines[iline][1]][dim]-coordinates[dim]) < 1.0E-8)
          {
            coordinates[dim] = vertexcoord[lines[iline][1]][dim];
            cout << "coordinates shifted to vertex 2" << endl;
          }
        }
      }

      //--------------------------------------------------------------------
      // transformation to global coordinates via shape functions (x = Nxi)
      //--------------------------------------------------------------------
//      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
//      static LINALG::Matrix<numnode,1> funct;
//      DRT::UTILS::shape_function_3D(funct,coordinates[0],coordinates[1],coordinates[2],cell->Ele()->Shape());
//
//      //jetzt braucht man noch die Knotenkoordinaten des zugehörigen Elements (global)
//      std::vector<double> globalcoordinates(3);
//
//      // initialization loop
//      for(int dim=0; dim<3; dim++)
//      {
//        globalcoordinates[dim] = 0.0;
//      }
//      // coordinate transform (isoparametric appraoch)
//      for(int inode=0; inode<numnode; inode++)
//      {
//        for(int dim=0; dim<3; dim++)
//        {
//          globalcoordinates[dim] += funct(inode) * cell->Ele()->Nodes()[inode]->X()[dim];
//        }
//      }

      // store coordinates of intersection point for each line
      intersectionpoints[iline] = coordinates;
    }
    else
    {
      //do nothing and go to the next line
    }
  }

  //TEST Ausgabe
//  for (std::map<int,std::vector<double> >::const_iterator iter = intersectionpoints.begin(); iter != intersectionpoints.end(); ++iter)
//  {
//    std::cout<< iter->first << std::endl;
//    std::vector<double> coord = iter->second;
//    for (std::size_t isd=0; isd<3; isd++)
//    {
//      std::cout<< coord[isd] << std::endl;
//    }
//  }

  // store intersection points in refinement cell
  cell->intersectionpoints_ = intersectionpoints;
  return;
}


/*------------------------------------------------------------------------------------------------*
 | capture flame front within one refinement cell to provide integration cells        henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::CaptureFlameFront(const Teuchos::RCP<const COMBUST::RefinementCell> cell)
{
  // The rootcell (= element) enters this function. If element is cut, a tree of refinement cells
  // belongs to it. They can be cut, or not.
  //
  // loop refinement cells of rootcell
  // {
  //   if cell is intersected
  //     call buildPLC() (-> Tetgen)
  //   else
  //    refinement cell is also an integration cell

//  std::cout << "CaptureFlameFront " << std::endl;

  std::vector<const COMBUST::RefinementCell* > RefinementCells;
  cell->SearchRefinementCells(RefinementCells);
//  std::cout << "Number of refinement cells " << RefinementCells.size() << std::endl;

  // vectors holding lists if domain integration cells and boundary integration surfaces, respectively
  GEO::DomainIntCells listDomainIntCellsperEle;
  GEO::BoundaryIntCells listBoundaryIntCellsperEle;

  //-------------------------------------------
  // loop over all refinement cells of root cell
  //-------------------------------------------
  for (std::size_t icell = 0; icell < RefinementCells.size(); icell++)
  {
    //---------------------
    // cell is intersected
    //---------------------
    if (RefinementCells[icell]->Intersected())
    {
//   std::cout << "Cell is intersected " << std::endl;
      const std::vector<double>& gfuncvalues = cell->GetGfuncValues();
      buildPLC(RefinementCells[icell], gfuncvalues, listDomainIntCellsperEle, listBoundaryIntCellsperEle);
    }
    //-------------------------
    // cell is not intersected
    //-------------------------
    // remark: if the refinement cell is not cut, it will be an integration cell of either side of
    //         the domain (refinement cell == integration cell)
    //         refinement cell can be identical with element, if refinement strategy is turned off
    else
    {
//    std::cout<< "add DomainIntCell for uncut refinement cell" <<std::endl;

      // get G-function values at vertices from rootcell cell
      const std::vector<double>& gfuncvalues = cell->GetGfuncValues();
      // get vertex coordinates (local fluid element coordinates) from refinement cell
      const std::vector<std::vector<double> >& vertexcoord = RefinementCells[icell]->GetVertexCoord();

      //---------------------------------
      // get global coordinates of nodes
      //---------------------------------
      const int numnode = cell->Ele()->NumNode();
      //const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      LINALG::SerialDenseMatrix xyze(3,numnode);

      for(int inode=0;inode<numnode;inode++)
      {
        xyze(0,inode) = cell->Ele()->Nodes()[inode]->X()[0];
        xyze(1,inode) = cell->Ele()->Nodes()[inode]->X()[1];
        xyze(2,inode) = cell->Ele()->Nodes()[inode]->X()[2];
      }

      //-----------------------------------------
      // get global coordinates of cell vertices
      //-----------------------------------------
      //corresponds to vertexcoord
      LINALG::SerialDenseMatrix cellcoord(3,numnode);
      // if cell == element, globalcellcoord == xyze
      LINALG::SerialDenseMatrix globalcellcoord(3,numnode);

      for(int ivert=0; ivert<numnode; ivert++)
      {
        static LINALG::Matrix<3,1> vertcoord;
        for(int dim=0; dim<3; dim++)
        {
          // cellcoord = transpose of vertexcoord
          cellcoord(dim,ivert) = vertexcoord[ivert][dim];
          // coordinates of a single cell vertex
          vertcoord(dim) = cellcoord(dim,ivert);
        }
        // transform vertex from local (element) coordinates to global (physical) coordinates
        GEO::elementToCurrentCoordinatesInPlace(cell->Ele()->Shape(), xyze, vertcoord);

        // store vertex in array of global cell coordinates
        for(int  dim=0; dim<3; dim++)
          globalcellcoord(dim,ivert) = vertcoord(dim);
        // if cell == element, globalcellcoord == xyze!
      }

      //---------------------------------
      // determine which domain of flame
      //---------------------------------
      // compute average G-function value for this refinement cell (= integration cell)
      bool inGplus = GetIntCellDomain(cellcoord, gfuncvalues, cell->Ele()->Shape(), DRT::Element::hex8);

      //TEST
      //std::cout << "globalcellcoord " << globalcellcoord(0,3) << globalcellcoord(1,3) << std::endl;
//      if(inGplus) {
//        std::cout << "In G plus" << std::endl;
//      }
//      else {
//        std::cout << "In G minus" << std::endl;
//      }
      //-- ----------------------
      // create integration cell
      //--- ---------------------
      // create a hex8 integration cell and add it to the list of integration cells per element
      listDomainIntCellsperEle.push_back(GEO::DomainIntCell(DRT::Element::hex8, cellcoord, globalcellcoord, inGplus));
    }
  }//end loop over all refinement cells
  //TEST
//  std::cout << cell->Ele()->Id() << "size of DomainIntCell: " << listDomainIntCellsperEle.size() << std::endl;
//  std::cout << cell->Ele()->Id() << "size of BoundaryIntCell: " << listBoundaryIntCellsperEle.size() << std::endl;

  //------------------------------------------------------------
  // store list of integration cells per element in flame front
  //------------------------------------------------------------
  myelementintcells_[cell->Ele()->Id()] = listDomainIntCellsperEle;
  // if there exist boundary integration cells
  if(listBoundaryIntCellsperEle.size() > 0)
    myboundaryintcells_[cell->Ele()->Id()] = listBoundaryIntCellsperEle;

  // TODO: temporary stuff -> clarify, remove later
  // hier muss jetzt die map "flamefrontpatches_" mit Oberflächenstückchen gefüllt werden
  // eventuell macht es keinen Sinn die übergebene Zelle "cell" als doppelt "const" zu deklarieren
  flamefrontpatches_.insert(make_pair(cell->Ele()->Id(),cell));

  return;
}


/*------------------------------------------------------------------------------------------------*
 | triangulate the interface (flame front) inside a refinement cell                   henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::TriangulateFlameFront(
       std::vector<std::vector<int> >&       trianglelist,
       std::multimap<int,std::vector<int> >& segmentlist,
       std::vector<std::vector<double> >&    pointlist,
       std::map<int,int>&                    intersectionpointsids,
       const std::vector<double>&            gfuncvalues
       )
{
  //std::cout << "TriangulateFlameFront" << std::endl;
  //check -> is cell really intersected? This means, not only edges of the cell are contained in segmentlist (key=-1)
  if(segmentlist.count(-1)<segmentlist.size()) //or !=
  {
    std::map<int,std::vector<int> > polygonpoints;
    int actpoint;

    // determine start segment for polygon
    // if a surface is intersected twice, the is the start segment
    //  -> so the second segment can not be forgotten
    //  -> and an edge is also not a good beginning
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

    //-----------------------------------------------------------------------
    // build polygon representing the boundary of the interface patch
    //-----------------------------------------------------------------------
    int j=0;
    //loop over all polygons
    for (std::multimap<int,std::vector<int> >::iterator segmentiter=segmentlist.equal_range(beginpolygon).first; segmentiter!=segmentlist.equal_range(beginpolygon).second; segmentiter++)
    {
      // int actkey = segmentiter->first;
      // std::vector<int> actsegmentpoints = segmentiter->second;

      // Its sufficent to determine the orientation of the first segment
      // to get the right direction of rotation (Umlaufsinn) of the polygon.
      IdentifyPolygonOrientation(segmentiter->second, segmentiter->first, intersectionpointsids, gfuncvalues);
      // IdentifyPolygonOrientation(actsegmentpoints, actkey, intersectionpointsids, gfuncvalues);
      std::vector<int> actsegmentpoints = segmentiter->second;
      //TEST
      // std::cout << "Startsegment " << segmentiter->first << std::endl;
      // std::cout << actsegmentpoints[0] << std::endl;
      // std::cout << actsegmentpoints[1] << std::endl;

      //contains the vertices of the polygon
      std::vector<int> polypoints;
      //store points of start segment
      polypoints.push_back(actsegmentpoints[0]);
      actpoint = actsegmentpoints[1];
      polypoints.push_back(actpoint);

      // now, the following segment is searched
      // this means the end point of this segment is equal the end point of the last segment
      // this procedure is repeated until the starting point is reached
      while(actpoint!=polypoints[0])
      {
        for (std::multimap<int,std::vector<int> >::const_iterator iter = segmentlist.begin(); iter != segmentlist.end(); ++iter)
        {
          std::vector<int> it_points = iter->second;
          if (it_points!=actsegmentpoints)
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
//      std::cout<<"polygonvertices"<< std::endl;
//      for (std::size_t ipoly=0; ipoly<polypoints.size(); ipoly++)
//      {
//        std::cout<< polypoints[ipoly] << std::endl;
//      }

      //URSULA 270709 ----- Modification
      //URSULA 310709 ----- Modification
      if (numpolygon==2)
      {
        if (polypoints.size()>(segmentlist.size()-3))
        {
//    		  std::cout<< "G-function values" << std::endl;
//    		  for (std::size_t i=0; i<gfuncvalues.size(); i++)
//    		  {
//    			  std::cout<< gfuncvalues[i] << std::endl;
//    		  }
//    		  dserror("Unexpected intersection");
          // special case: both segments of the double intersected surface belong to the same polygon
          break;
        }
      }
    }

    //--------------------------------------------------------------------------
    // build triangles
    //--------------------------------------------------------------------------
    for (std::size_t ipolygons=0; ipolygons<polygonpoints.size(); ipolygons++)
    {
      std::vector<int> polypoints = polygonpoints[ipolygons];
      if (polypoints.size()<4)
        dserror("TriangulateFlameFront needs at least 3 intersectionpoints");

      if (polypoints.size()==4) //there are only tree different points and tree points form a triangle
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
        if (numpoints%2==0) // even
        {
          point1 = pointlist[polypoints[0]];
          point2 = pointlist[polypoints[numpoints/2]];
        }
        else // odd
        {
          point1 = pointlist[polypoints[0]];
          point2 = pointlist[polypoints[(numpoints+1)/2]];
        }
        for (int dim=0; dim<3; dim++)
        {
          midpoint[dim] = (point2[dim] + point1[dim]) * 0.5;
          //URSULA 310709 ------ Modification
          /*This is to avoid problems with TetGen, when triangles lie competely on a surface
           */
          if (midpoint[dim]==1 or midpoint[dim]==-1) // midpoint lies on surface of the cell
        	  midpoint[dim] = midpoint[dim] * 0.98;  // -> move midpoint into th cell
        }
        std::size_t midpoint_id = pointlist.size();//ids start at 0
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
//    std::cout<<"triangles"<< std::endl;
//    for (std::size_t itriangle=0; itriangle<trianglelist.size(); itriangle++)
//    {
//      std::cout<< "dreieck " << itriangle<< std::endl;
//     for (int i=0; i<3; i++)
//       std::cout<< trianglelist[itriangle][i]<<std::endl;
//    }
  }
  return;
}


/*------------------------------------------------------------------------------------------------*
 | identify the orientation of the interface polygon, normal vector + -> -                        |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::IdentifyPolygonOrientation(
       std::vector<int>&           segment,
       const int                   surf_id,
       std::map<int,int>&          intersectionpointsids,
       const std::vector<double>&  gfuncvalues
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
    int interpoint = segment[0];

    for(std::map<int,int>::const_iterator iterinterpoint=intersectionpointsids.begin(); iterinterpoint!=intersectionpointsids.end(); iterinterpoint++)
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
  else //segment == diagonal
  {
    for(int k=0; k<4; k++)
    {
      if(testnode[surf_id][k]==segment[0])
      {
        if(k==0)
          testpoint = 3;
        else
          testpoint = k-1;
        break; // TODO: Ist das break hier OK?
      }
    }
  }
  if(gfuncvalues[testnode[surf_id][testpoint]]>0) //change! (normal vector in domain plus <0 is needed)
  {
    int temp = segment[0];
    segment[0] = segment[1];
    segment[1] = temp;
  }
  return;
}


/*------------------------------------------------------------------------------------------------*
 | build polygon segments of intersection surface of a refinement cell                            |
 | hex8 only                                                                                      |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::buildFlameFrontSegments(
        std::map<int,int>&                       intersectionpointsids,
        std::multimap<int,std::vector<int> >&    segmentlist,
        const std::vector<double>&               gfuncvalues,
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
  //hex8 only
  for (int i=0; i<6; i++)//loop over all surfaces
  {
    // gets end points of the segment
    std::vector<int> segmentpoints;

    for (int j=0; j<4; j++)//loop over all lines of the surface
    {
      // find line
      std::map<int,int>::const_iterator it_intersectionpoint = intersectionpointsids.find(surface[i][j]);
      if (it_intersectionpoint!=intersectionpointsids.end()) //check: is line intersected?
      {
        //get id of the intersectionpoint and store it in segmentpoints vector
        segmentpoints.push_back(it_intersectionpoint->second);
      }
    }

    //std::cout << "segmentpoints size" << segmentpoints.size() << std::endl;

    switch (segmentpoints.size())
    {
    case 0:
    {
      // check, if gfunc==0 at the vertices
      // gets nodes with gfunc==0
      std::vector<int> zeropoints;
      std::vector<std::vector<int> > surfacepointlist = DRT::UTILS::getEleNodeNumberingSurfaces(DRT::Element::hex8);
      for(int k=0; k<4; k++)//loop over nodes
      {
        if (gfuncvalues[surfacepointlist[i][k]]==0)
          zeropoints.push_back(surfacepointlist[i][k]);
      }
      //std::cout << "zeropoints size" << zeropoints.size() << std::endl;
      switch (zeropoints.size())
      {
      case 0: //surface not intersected
      case 1: //surface not intersected, but one vertex touchs the interface
      {
        break;
      }
      case 2:
        /*
         * 2 opposite vertices with Phi=0
         * - surface is intersected, iff Phi at the remaining vertices has different sign (segment == diagonal)
         * - otherwise vertiches touch interface
         * 2 neighboring vertices with Phi=0
         * - edge touchs interface, but cell is not intersected
         * - edge limits interface patch -> then, it is needed in TriangulateFlameFront()
         */
      {
        // intersection: opposite vertices with Phi=0
        if(((zeropoints[0]==surfacepointlist[i][0])and(zeropoints[1]==surfacepointlist[i][2])))
        {
          // different sign of Phi at the remaining vertices
          if (gfuncvalues[surfacepointlist[i][1]]*gfuncvalues[surfacepointlist[i][3]]<0)
          {
            //store in segmentlist
            segmentlist.insert(pair<int,std::vector<int> >(i,zeropoints));
          }
        }
        else if (((zeropoints[0]==surfacepointlist[i][1])and(zeropoints[1]==surfacepointlist[i][3])))
        {
          if (gfuncvalues[surfacepointlist[i][0]]*gfuncvalues[surfacepointlist[i][2]]<0)
          {
            segmentlist.insert(pair<int,std::vector<int> >(i,zeropoints));
          }
        }
        else
        {
          // no intersection
          // possibly segment limits the interface
          // store in segment list (key=-1), if it is not already contained in the segmentlist
          if (zeropoints[0]>zeropoints[1])
          {
            int temp =zeropoints[1];
            zeropoints[1] = zeropoints[0];
            zeropoints[0] = temp;
          }
          //already in segmentlist?
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
      case 3:
      {
        // surface is not intersected, but two following edges are aligned with the interface
        // store in segment list (key=-1), if it is not already contained in the segmentlist

        // zeropoints are assigned to segments
        // checke which vertex is missing
        if(zeropoints[0]==surfacepointlist[i][0])
        {
          if(zeropoints[1]==surfacepointlist[i][1])
          {
            if(zeropoints[2]==surfacepointlist[i][2])
            {
              //vertex 4 is missing -> no problem
            }
            else
            {
              //vertex 3 is missing -> rearrange zeropoints
              int temp = zeropoints[0];
              zeropoints[0] = zeropoints[2];
              zeropoints[2] = zeropoints[1];
              zeropoints[1] = temp;
            }
          }
          else
          {
            //vertex 2 is missing -> rearrange zeropoints
            int temp = zeropoints[0];
            zeropoints[0] = zeropoints[1];
            zeropoints[1] = zeropoints[2];
            zeropoints[2] = temp;
          }
        }
        else
        {
          //vertex 1 is missing -> no problem
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
      case 4:
        // suface is aligned with the interface -> cell is not intersected
        // surface has to be added to BoundaryIntCells
      {

        // as all segments, limiting the interface, are already found in case 2 and 3,
        // nothing remains to do here

        break;
      }
      default:
        dserror("impossible number of zero values");
      }
      break;
    }
    case 1:
    {
      //std::cout << "one intersection point" << std::endl;
      //find zeropoint to build segment
      std::vector<int> zeropoints;
      std::vector<std::vector<int> > surfacepointlist = DRT::UTILS::getEleNodeNumberingSurfaces(DRT::Element::hex8);
      for(int k=0; k<4; k++)//loop over nodes
      {
        if (gfuncvalues[surfacepointlist[i][k]]==0)
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
      //std::cout << "4 intersection points " << std::endl;
      // 2 segments per surface are possible for hex8 -> 4 intersection points
      // idea to assign the intersection points to the segments:
      //       look for maximal distance
      //       this combination is not possible
      //       the segments are now obtained
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

      //build segment
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


/*------------------------------------------------------------------------------------------------*
 | build piecewise linear complex (PLC) in Tetgen format                              henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::buildPLC(
        //const Teuchos::RCP<const COMBUST::RefinementCell> cell,
        const COMBUST::RefinementCell* cell,
        const std::vector<double>& gfuncvaluesrootcell,
        GEO::DomainIntCells& domainintcelllist,
        GEO::BoundaryIntCells& boundaryintcelllist)
{
 /*---------------------------------------------------------------------------------
  * a PLC contains
  * - points -> vertices, intersection points, center of interface patch
  * - XFEMSurfaces -> surfaces of refinement cell
  * - trianglelist -> triangles approximating the interface: TriangulateFlameFront()
  * - segmentlist -> intersection curve between interface and surfaces of cell: buildFlameFrontSegments()
  * input of CreateIntegrationCells()
  *---------------------------------------------------------------------------------*/

  std::vector<std::vector<double> >    pointlist;
  std::multimap<int,std::vector<int> > segmentlist;
  std::vector<std::vector<int> >       trianglelist;
  // TODO: sollte hier eigentlich nicht mehr nötig sein da ich nur die Nummerierung der Ecken möchte
  //       ich lasse es trotzdem mal drin

  const DRT::Element::DiscretizationType distype = cell->Ele()->Shape();
  // global coordinates of element this cell belongs to
  const int numnode = cell->Ele()->NumNode();
  LINALG::SerialDenseMatrix  xyze(3,numnode);
  for(int inode=0;inode<numnode;inode++)
  {
    xyze(0,inode) = cell->Ele()->Nodes()[inode]->X()[0];
    xyze(1,inode) = cell->Ele()->Nodes()[inode]->X()[1];
    xyze(2,inode) = cell->Ele()->Nodes()[inode]->X()[2];
  }

  // node ids of each surface
  // here used as point ids
  const std::vector<std::vector<int> >& xfemsurfacelist = DRT::UTILS::getEleNodeNumberingSurfaces(distype);
  // corresponds to intersection points, contains point ids instead of coordinates
  std::map<int,int> intersectionpointsids;

  // get vertex coordinates (local fluid element coordinates) from refinement cell
  const std::vector<std::vector<double> >& vertexcoord = cell->GetVertexCoord();
  // get intersection points from refinement cell
  const std::map<int,std::vector<double> >& intersectionpoints = cell->intersectionpoints_;

  // get G-function values from refinement cell
  const std::vector<double>& gfuncvalues = cell->GetGfuncValues();

//TEST
//  std::cout << cell->Ele()->Id() << std::endl;
//  std::cout<< "G-Werte" << std::endl;
//  for (std::size_t i=0; i<gfuncvalues.size(); i++)
//  {
//    std::cout<< gfuncvalues[i] << std::endl;
//  }

  //---------------------------------------------
  // fill list of points and intersection points
  //---------------------------------------------
  // das ist allgemein gültig, egal welcher Elementtyp
  // vertices of the cell
  int numofpoints = 0;
  int numvertex = DRT::UTILS::getNumberOfElementCornerNodes(distype);
  for (int ivertex=0; ivertex<numvertex; ivertex++)
  {
    // remark: Hex hat 8 Ecken, weitere Koord in vertexcoord bei hex20-27 sind innere Knoten
    //         und keine Ecken!
    pointlist.push_back(vertexcoord[ivertex]);
    numofpoints++;
  }
  // intersection points
  for (std::map<int,std::vector<double> >::const_iterator iter = intersectionpoints.begin(); iter != intersectionpoints.end(); ++iter)
  {
    pointlist.push_back(iter->second);
    intersectionpointsids[iter->first] = numofpoints;
    numofpoints++;
  }
  //TEST
//    for (std::size_t iter=0; iter<pointlist.size(); ++iter)
//    {
//      std::cout<< iter << std::endl;
//      std::vector<double> coord = pointlist[iter];
//      for (std::size_t isd=0; isd<3; isd++)
//      {
//        std::cout<< coord[isd] << std::endl;
//      }
//    }

  //------------------------------------------------------------------
  // build segments enclosing interface patches within a cell/element
  //------------------------------------------------------------------
  /*
   * Segments form the boundary of the interfacepatch in the element.
   * Here, some special cases have to be distinguished
   */
  if(distype == DRT::Element::hex8) {
    buildFlameFrontSegments(intersectionpointsids, segmentlist, gfuncvalues, pointlist);
  }
  else {
    dserror("buildFlameFrontSegments() relies on hex8 elements!");
  }

  //TEST
//  std::cout<<"Segments"<< std::endl;
//  for (std::map<int,std::vector<int> >::const_iterator iter = segmentlist.begin(); iter != segmentlist.end(); ++iter)
//  {
//    std::cout<<"Segment "<< iter->first << std::endl;
//    std::vector<int> point = iter->second;
//    for (std::size_t isd=0; isd<2; isd++)
//    {
//      std::cout<< point[isd] << std::endl;
//    }
//  }

  //-----------------------------------------------
  // triangulate flame front within a cell/element
  //-----------------------------------------------
  if(distype == DRT::Element::hex8) {
    TriangulateFlameFront(trianglelist, segmentlist, pointlist, intersectionpointsids, gfuncvalues);
  }
  else {
    dserror("TriangulateFlameFront() relies on hex8 elements!");
  }

  //TEST
//  for (std::size_t iter=0; iter<pointlist.size(); ++iter)
//  {
//    std::cout<< iter << std::endl;
//    std::vector<double> coord = pointlist[iter];
//    for (std::size_t isd=0; isd<3; isd++)
//    {
//      std::cout<< coord[isd] << std::endl;
//    }
//  }

  //--------------------------
  // cell does have triangles
  //--------------------------
  // remark: Normally, a refinement cell will have triangles if it is cut. However, if a single
  //         point lies on the cell surface, or the flame front is aligned with the cell surface, it
  //         is "cut" (cell->Intersected() == true), but has no triangles in its list.
  if (trianglelist.size() > 0)
  {
    // if cell is intersected, there have to be triangles in the list
    dsassert(trianglelist.size() != 0,"The triangle list for this cell is empty!\n");
#ifdef QHULL
    //--------------------------------------------------------------
    // create and store domain integration cells for a cell/element
    //--------------------------------------------------------------
    // remark: if cell is cut, tetrahedra have to be created serving as integration cells
    CreateIntegrationCells(pointlist, segmentlist, xfemsurfacelist, trianglelist,
                           domainintcelllist, xyze, gfuncvaluesrootcell);
#else
    dserror("Set QHULL flag to use Tetgen!");
#endif
    //----------------------------------
    // store boundary integration cells
    //----------------------------------
    // remark: Although Tetgen also delivers triangles as interface patches, these are not used,
    //         because their direction of rotation (Umlaufsinn) is not known. Instead, the triangles
    //         from the trianglelist are used directly.
    for (std::size_t itriangle=0; itriangle<trianglelist.size(); itriangle++)
    {
      LINALG::SerialDenseMatrix trianglecoord(3,3); //3 directions, 3 nodes
      LINALG::SerialDenseMatrix phystrianglecoord(3,3);
      for (int inode=0; inode<3; inode++)
      {
        static LINALG::Matrix<3,1> tcoord;
        for (int dim=0; dim<3; dim++)
        {
          trianglecoord(dim,inode) = pointlist[trianglelist[itriangle][inode]][dim];
          tcoord(dim) = trianglecoord(dim,inode);
        }
        GEO::elementToCurrentCoordinatesInPlace(distype, xyze, tcoord);
        for(int  dim=0; dim<3; dim++)
          phystrianglecoord(dim,inode) = tcoord(dim);
      }
      //store boundary integration cells in boundaryintcelllist
      boundaryintcelllist.push_back(GEO::BoundaryIntCell(DRT::Element::tri3, -1, trianglecoord,
                                    Teuchos::null, phystrianglecoord));
    }
  }
  //------------------------------
  // cell does not have triangles
  //------------------------------
  // remark: if the refinement cell is not cut (does not have triangles), it will be an integration
  //         cell of either side of the domain (refinement cell == integration cell) refinement cell
  //         can be identical with element, if refinement strategy is turned off
  else
  {
//    std::cout<< "add DomainIntCell for uncut refinement cell" <<std::endl;
    //--------------------------------------------------------------
    // create and store domain integration cells for a cell/element
    //--------------------------------------------------------------
    //corresponds to vertexcoord
    LINALG::SerialDenseMatrix cellcoord(3, numnode);
    //physical cell coordinates
    LINALG::SerialDenseMatrix physcellcoord(3, numnode);
    for(int inode=0; inode<numnode; inode++)
    {
      static LINALG::Matrix<3,1> ccoord;
      for(int dim=0; dim<3; dim++)
      {
        cellcoord(dim,inode) = vertexcoord[inode][dim];
        ccoord(dim) = cellcoord(dim,inode);
      }
      GEO::elementToCurrentCoordinatesInPlace(distype, xyze, ccoord);
      for(int  dim=0; dim<3; dim++)
        physcellcoord(dim,inode) = ccoord(dim);
    }
    // determine which domain of flame
    bool inGplus = GetIntCellDomain(cellcoord, gfuncvaluesrootcell, distype, DRT::Element::hex8);
    //TEST
//    std::cout << "physcellcoord " << physcellcoord(0,3) << physcellcoord(1,3) << std::endl;
//    if(inGplus) {
//      std::cout << "In G plus" << std::endl;
//    }
//    else {
//      std::cout << "In G minus" << std::endl;
//    }
    domainintcelllist.push_back(GEO::DomainIntCell(distype, cellcoord, physcellcoord, inGplus));

    //----------------------------------
    // store boundary integration cells
    //----------------------------------
    // remark: Boundary integartion cells are available, if one or more surfaces of the cell are aligned
    //         with the flame front. However, to avoid to have boundary inegration cells twice, they are
    //         only stored, if their normal vector points from + to -.
    if(segmentlist.size()>=4) //hex8 only
    {
      if (cell->Ele()->Shape()!=DRT::Element::hex8)
        dserror("buildPLC() not supported for this element shape!");
      if(inGplus)
      {
        int numsurf = DRT::UTILS::getNumberOfElementSurfaces(DRT::Element::hex8);
        for(int isurf=0; isurf<numsurf; isurf++)
        {
          int numzeronodes = 0;
          for(int k=0; k<4; k++)//loop over nodes
          {
            if (gfuncvalues[xfemsurfacelist[isurf][k]]==0)
              numzeronodes++;
            //if (gfuncvalues[xfemsurfacelist[isurf][k]]<0)
            //  influidplus = false;
          }
          if (numzeronodes==4)
          {
            //store boundary integration cells in boundaryintcelllist
//            std::cout<< "store  BoundIntCell " <<std::endl;
            //std::cout<< isurf <<std::endl;
            LINALG::SerialDenseMatrix quadcoord(3,4); //3 directions, 4 nodes
            LINALG::SerialDenseMatrix physquadcoord(3,4);
            for (int inode=0; inode<4; inode++)
            {
              static LINALG::Matrix<3,1> qcoord;
              for (int dim=0; dim<3; dim++)
              {
                quadcoord(dim,inode) = pointlist[xfemsurfacelist[isurf][inode]][dim];
                qcoord(dim) = quadcoord(dim,inode);
              }
              GEO::elementToCurrentCoordinatesInPlace(DRT::Element::hex8, xyze, qcoord);
              for(int  dim=0; dim<3; dim++)
                physquadcoord(dim,inode) = qcoord(dim);
            }
            boundaryintcelllist.push_back(GEO::BoundaryIntCell(DRT::Element::quad4, -1, quadcoord,
                                          Teuchos::null, physquadcoord));
          }
        }
      }
    }
  }
  return;
}


#ifdef QHULL
/*------------------------------------------------------------------------------------------------*
 | calls the CDT to create burnt and unburnt integration cells                        henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::CreateIntegrationCells(
       const std::vector<std::vector<double> >&  pointlist,
       std::multimap<int,std::vector<int> >&     segmentlist,
       const std::vector<std::vector<int> >&     xfemsurfacelist,
       const std::vector<std::vector<int> >&     trianglelist,
       GEO::DomainIntCells&                      domainintcelllist,
       //GEO::BoundaryIntCells&                    boundaryintcelllist,
       const LINALG::SerialDenseMatrix           xyze,
       const std::vector<double>&                gfuncvalues
        ) // output: integration cells in Tetgen format
{
  /*---------------------------------------------------------------------------
   * tetgenio in is filled with the PLC
   * tetgenio out is filled by TetGen with tetrahedras
   *---------------------------------------------------------------------------*/

  const int dim = 3;
  tetgenio in;
  tetgenio out;
  char switches[] = "pQ";    //- p     tetrahedralizes a PLC
                               //-Q      no terminal output except errors
                               // YY     do not generate additional points on surfaces -> fewer cells
  tetgenio::facet *f;
  tetgenio::polygon *p;

  //allocate point list
  in.numberofpoints = pointlist.size();
  in.pointlist = new REAL[in.numberofpoints * dim];

  // TODO: should scale factor be an input parameter? Is it neccessary?
  const double scalefactor = 1.0E7;

  // fill point list
  int fill = 0;
  for(int i = 0; i <  in.numberofpoints; i++)
  {
    for(int j = 0; j < dim; j++)
    {
      double coord = pointlist[i][j] * scalefactor;
      in.pointlist[fill] = (REAL) coord;
      fill++;
    }
  }

  in.numberoffacets = xfemsurfacelist.size() + trianglelist.size();
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];

  // loop over all element surfaces
  for(std::size_t i=0; i<xfemsurfacelist.size(); i++)
  {
    f = &in.facetlist[i];
    //map
    // int numsegments = 0;
    // std::map<int,std::vector<int> >::const_iterator it_segmentlist = segmentlist.find(i);
    // if (it_segmentlist!=segmentlist.end()) //check: is surface intersected?
    // {
    //   numsegments = 1;
    // }
    // falls zweimal geschnitten : numpolygons = 3
    //multimap
    int numsegments = (int) segmentlist.count(i);
    f->numberofpolygons = 1 + numsegments;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    // face of hex
    p = &f->polygonlist[0];
    p->numberofvertices = 4; //hexahedra has 4 verices
    p->vertexlist = new int[p->numberofvertices];
    for(int ivertex = 0; ivertex < p->numberofvertices; ivertex++)
    {
      p->vertexlist[ivertex] = xfemsurfacelist[i][ivertex];
    }
    //store segments if present
    if (numsegments)
    {
      //map
      // for(int j=1; j<(numsegments+1); j++) //vorsicht bei Obergrenze
      // {
      //   p = &f->polygonlist[j];
      //   p->numberofvertices = 2;
      //   p->vertexlist = new int[p->numberofvertices];
      //   for(int k=0; k<2; k++)
      //   {
      //     p->vertexlist[k] = segmentlist[i][k]; //surface i, point k
      //TEST
      //     std::cout << "Surface " << i << std::endl;
      //     std::cout << p->vertexlist[k] << std::endl;
      //   }
      // }
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
  std::size_t k = xfemsurfacelist.size();
  for(std::size_t i=0; i<trianglelist.size(); i++)
  {
    f = &in.facetlist[k];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 3; //triangle has 3 vertices
    p->vertexlist = new int[p->numberofvertices];
    for(int j=0; j<3; j++)
    {
      p->vertexlist[j] = trianglelist[i][j];
      //TEST
// std::cout << "Triangle " << i << std::endl;
// std::cout << k << std::endl;
// std::cout << p->vertexlist[j] << std::endl;
    }
    //in.facetmarkerlist[k] = 1; nur für BoundaryIntCells
    k++;
  }

//  std::cout << "----- Start Tetgen -----" << std::endl;
  //in.save_nodes("tetin");
  //in.save_poly("tetin");
  tetrahedralize(switches, &in,  &out);
  //out.save_nodes("tetout");
  //out.save_elements("tetout");
  //out.save_faces("tetout");
//  std::cout << "----- End Tetgen -----" << std::endl;

  // store interface triangles (+ recovery of higher order meshes)
  fill = 0;
  for(int i = 0; i <  out.numberofpoints; i++)
    for(int j = 0; j < dim; j++)
    {
      out.pointlist[fill] = (REAL) (out.pointlist[fill] * (1.0/scalefactor));
      fill++;
    }

  TransformIntegrationCells(out, domainintcelllist, xyze, gfuncvalues);

  return;
}
#endif


#ifdef QHULL
/*------------------------------------------------------------------------------------------------*
 | this could easlily be integrated into CreateIntegrationCells()                     henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::TransformIntegrationCells(
       tetgenio&                       out,
       GEO::DomainIntCells&            domainintcelllist,
       //GEO::BoundaryIntCells&          boundaryintcelllist,
       const LINALG::SerialDenseMatrix xyze,
       const std::vector<double>&      gfuncvalues
       ) // output: integration cells in baci format
{
  DRT::Element::DiscretizationType distype = DRT::Element::tet4;
  const int numTetNodes = DRT::UTILS::getNumberOfElementNodes(distype);

  //std::cout << "number of Tertrahedra " << out.numberoftetrahedra << std::endl;
  for(int i=0; i<out.numberoftetrahedra; i++ )
  {
    LINALG::SerialDenseMatrix tetrahedroncoord(3, numTetNodes);
    LINALG::SerialDenseMatrix phystetrahedroncoord(3, numTetNodes);
    for(int j = 0; j < numTetNodes; j++)
    {
      static LINALG::Matrix<3,1> tetcoord;
      for(int k = 0; k < 3; k++)
      {
        tetrahedroncoord(k,j) = out.pointlist[out.tetrahedronlist[i*out.numberofcorners+j]*3+k];
        // k: drei aufeinanderfolgende Einträge in der pointlist, entspricht den 3 Richtungen
        // *3: Richtungen je Knoten
        // tetraherdronlist[i*out.numberofcorners+j]: node ID for ith element, jth node
        tetcoord(k) = tetrahedroncoord(k,j);
      }
      // compute physical coordinates
      GEO::elementToCurrentCoordinatesInPlace(DRT::Element::hex8, xyze, tetcoord);
      for(int k = 0; k < 3; k++)
        phystetrahedroncoord(k,j) = tetcoord(k);
    }

    // if degenerated don't store
    if(!GEO::checkDegenerateTet(numTetNodes, tetrahedroncoord, phystetrahedroncoord))
    {
       bool inGplus = GetIntCellDomain(tetrahedroncoord, gfuncvalues, DRT::Element::hex8, distype);
       if(inGplus)
       {
//         std::cout << "In G plus" << std::endl;
       }
       else
       {
//         std::cout << "In G minus" << std::endl;
       }
       domainintcelllist.push_back(GEO::DomainIntCell(distype, tetrahedroncoord, phystetrahedroncoord, inGplus));
    }
  }

  return;
}
#endif


/*------------------------------------------------------------------------------------------------*
 | compute average GfuncValue for integration cell                                                |
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::FlameFront::GetIntCellDomain(
       const LINALG::SerialDenseMatrix        IntCellCoord,
       const std::vector<double>&             gfuncvalues,
       const DRT::Element::DiscretizationType xfem_distype,
       const DRT::Element::DiscretizationType cell_distype
       )
{
  /*------------------------------------------------------------------------------
   * - compute phi at each node of domain integration cell
   * - compute average G-value for each domain integration cell
   * -> >=0 in plus == true
   * _> <0  in minus == false
   * ----------------------------------------------------------------------------*/

  bool inGplus = false;

  int numcellnodes = 0;
  switch(cell_distype)
  {
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
  if (xfem_distype==DRT::Element::hex8) {
    numelenodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
  }
  else {
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
      gvalcellnodes[icellnode] += funct(ielenode) * gfuncvalues[ielenode];
    }
  }

  //calculate average Gfunc value
  double averageGvalue = 0.0;

  for (int icellnode=0; icellnode<numcellnodes; icellnode++)
    averageGvalue += gvalcellnodes[icellnode];

  if (numcellnodes == 0.0)
    dserror("division by zero: number of vertices per cell is 0!");
  averageGvalue /= numcellnodes;

  //determine DomainIntCell position
  //URSULA 030809 ----- Modfication
  // >= ; because of the approximation of the interface,
  // also averageGvalue=0 is possible; as always approximation errors are made,
  // this is just one more
  if(averageGvalue>=0.0)
    inGplus = true;
//  if(averageGvalue==0)
//    dserror("can't determine DomainIntCell position");

  return inGplus;
}


#ifdef PARALLEL
/*------------------------------------------------------------------------------------------------*
 | export flame front to all processors                                               henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ExportFlameFront(std::map<int, GEO::BoundaryIntCells>& myflamefront)
{
  //-------------------------------
  // prepare parallel communication
  //-------------------------------
  // get communicator, e.g. from fluid discretization
  const Epetra_Comm& comm = fluiddis_->Comm();
  const int myrank = comm.MyPID();
  const int numproc = comm.NumProc();

  int size_one = 1;

  DRT::Exporter exporter(comm);

  // destination proc (the "next" one)
  int dest = myrank+1;
  if(myrank == (numproc-1))
    dest = 0;

  // source proc (the "last" one)
  int source = myrank-1;
  if(myrank == 0)
    source = numproc-1;

//#ifdef DEBUG
//  cout << "number of cell groups (cut column elements) on proc " << myrank << " " << myflamefront.size() << endl;
//  for(std::map<int, GEO::BoundaryIntCells>::const_iterator cellgroup=myflamefront.begin(); cellgroup != myflamefront.end(); ++cellgroup)
//    {
//      // put ID of cut element here
//      if (cellgroup->first == 1274)
//      {
//        cout << "output for element 1274 packing" << endl;
//        const size_t numcells = (cellgroup->second).size();
//        cout << "proc " << myrank << " number of integration cells: " << numcells << endl;
//        for (size_t icell=0; icell<numcells; ++icell)
//        {
//          GEO::BoundaryIntCell cell = cellgroup->second[icell];
//          cout << "proc " << myrank << " cell " << icell << " vertexcoord " << cell.CellNodalPosXYZ();
//        }
//      }
//    }
//#endif

#ifdef DEBUG
  std::cout << "proc " << myrank << " number of flame front pieces available before export " << myflamefront.size() << std::endl;
#endif

  //-----------------------------------------------------------------
  // pack data (my boundary integration cell groups) for initial send
  //-----------------------------------------------------------------
  vector<char> dataSend;
  COMBUST::FlameFront::packBoundaryIntCells(myflamefront, dataSend);

  //-----------------------------------------------
  // send data around in a circle to all processors
  //-----------------------------------------------
  // loop over processors
  for(int num = 0; num < numproc-1; num++)
  {
    vector<int> lengthSend(1,0);
    lengthSend[0] = dataSend.size();

#ifdef DEBUG
    cout << "--- sending "<< lengthSend[0] << " bytes: from proc " << myrank << " to proc " << dest << endl;
#endif

    // send length of the data to be received ...
    MPI_Request req_length_data;
    int length_tag = 0;
    exporter.ISend(myrank, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);
    // ... and receive length
    vector<int> lengthRecv(1,0);
    exporter.Receive(source, length_tag, lengthRecv, size_one);
    exporter.Wait(req_length_data);

    // send actual data ...
    int data_tag = 4;
    MPI_Request req_data;
    exporter.ISend(myrank, dest, &(dataSend[0]), lengthSend[0], data_tag, req_data);
    // ... and receive data
    vector<char> dataRecv(lengthRecv[0]);
    exporter.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
    exporter.Wait(req_data);

#ifdef DEBUG
    cout << "--- receiving "<< lengthRecv[0] << " bytes: to proc " << myrank << " from proc " << source << endl;
#endif

    //-----------------------------------------------
    // unpack data (boundary integration cell groups)
    //-----------------------------------------------
    std::map<int, GEO::BoundaryIntCells> flamefront_recv;

    COMBUST::FlameFront::unpackBoundaryIntCells(dataRecv, flamefront_recv);

//#ifdef DEBUG
//    cout << "proc " << myrank << " receiving "<< lengthRecv[0] << " bytes from proc " << source << endl;
//    cout << "proc " << myrank << " receiving "<< flamefront_recv.size() << " flame front pieces from proc " << source << endl;
//
//    for(std::map<int, GEO::BoundaryIntCells>::const_iterator cellgroup=flamefront_recv.begin(); cellgroup != flamefront_recv.end(); ++cellgroup)
//    {
//      // put ID of cut element here
//      if (cellgroup->first == 1274)
//      {
//        cout << "output for element 1274 unpacking" << endl;
//        const size_t numcells = (cellgroup->second).size();
//        cout << "proc " << myrank << " number of integration cells: " << numcells << endl;
//        for (size_t icell=0; icell<numcells; ++icell)
//        {
//          GEO::BoundaryIntCell cell = cellgroup->second[icell];
//          cout << "proc " << myrank << " cell " << icell << " vertexcoord " << cell.CellNodalPosXYZ();
//        }
//      }
//    }
//#endif

    // add group of cells to my flame front map
    // remark: all groups of boundary integration cells (flame front pieces within an element) are collected here
    for (std::map<int, GEO::BoundaryIntCells>::const_iterator cellgroup = flamefront_recv.begin(); cellgroup != flamefront_recv.end(); ++cellgroup)
    {
      myflamefront.insert(*cellgroup);
    }

    // make received data the new 'to be sent' data
    dataSend = dataRecv;

    // processors wait for each other
    comm.Barrier();
  }
#ifdef DEBUG
  std::cout << "proc " << myrank << " number of flame front pieces now available " << myflamefront.size() << std::endl;
#endif
}
#endif


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void COMBUST::FlameFront::packBoundaryIntCells(
    const std::map<int, GEO::BoundaryIntCells>& intcellmap,
    vector<char>&                               dataSend)
{
  // pack data on all processors
  // loop entries of map (groups of boundary integration cells)
  for(std::map<int, GEO::BoundaryIntCells>::const_iterator cellgroup=intcellmap.begin(); cellgroup != intcellmap.end(); ++cellgroup)
  {
    vector<char> data;
    data.resize(0);

    // pack data of all boundary integrations cells belonging to an element
    const int elegid = cellgroup->first;
    DRT::ParObject::AddtoPack(data,elegid);

    const size_t numcells = (cellgroup->second).size();
    DRT::ParObject::AddtoPack(data,numcells);

    for (size_t icell=0; icell<numcells; ++icell)
    {
      GEO::BoundaryIntCell cell = cellgroup->second[icell];
      // get all member variables from a single boundary integration cell
      // distype of cell
      const DRT::Element::DiscretizationType distype = cell.Shape();
      DRT::ParObject::AddtoPack(data,distype);

      // coordinates of cell vertices in (fluid) element parameter space
//      const Epetra_SerialDenseMatrix& vertices_xi = cell.CellNodalPosXiDomain();
      const LINALG::SerialDenseMatrix& vertices_xi = cell.CellNodalPosXiDomain();
      DRT::ParObject::AddtoPack(data,vertices_xi);

      // coordinates of cell vertices in physical space
//      const Epetra_SerialDenseMatrix& vertices_xyz = cell.CellNodalPosXYZ();
      const LINALG::SerialDenseMatrix& vertices_xyz = cell.CellNodalPosXYZ();
      DRT::ParObject::AddtoPack(data,vertices_xyz);
    }

    // write packed boundary integration cell to dataSend
    DRT::ParObject::AddtoPack(dataSend,data);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void COMBUST::FlameFront::unpackBoundaryIntCells(
    const vector<char>&                     dataRecv,
    std::map<int, GEO::BoundaryIntCells>&   intcellmap)
{
  // pointer to current position of group of cells in global string (counts bytes)
  int posofgroup = 0;

  while (posofgroup < (int) dataRecv.size())
  {
    // pointer to current position in a group of cells in local string (counts bytes)
    int posingroup = 0;
    vector<char> data;
    DRT::ParObject::ExtractfromPack(posofgroup, dataRecv, data);

    // extract fluid element gid
    int elegid = -1;
    DRT::ParObject::ExtractfromPack(posingroup,data,elegid);
    if (elegid < 0) dserror("extraction of element gid failed");

    //extract number of boundary integration cells for this element
    int numvecs = -1;
    DRT::ParObject::ExtractfromPack(posingroup,data,numvecs);

    // vector holding group of boundary integration cells belonging to this element
    GEO::BoundaryIntCells intcellvector;

    for (int icell=0; icell<numvecs; ++icell)
    {
      //--------------------------------------------------------------------
      // extract all member variables for a single boundary integration cell
      //--------------------------------------------------------------------
      // distype of cell
      DRT::Element::DiscretizationType distype;
      int distypeint = -1;
      DRT::ParObject::ExtractfromPack(posingroup,data,distypeint);
      if (distypeint == 4)
        distype = DRT::Element::tri3;
      else if (distypeint == 1)
        distype = DRT::Element::quad4;
      else
        dserror("unexpected distype");

      LINALG::SerialDenseMatrix vertices_xi;
      DRT::ParObject::ExtractfromPack(posingroup,data,vertices_xi);

      // coordinates of cell vertices in physical space
      LINALG::SerialDenseMatrix vertices_xyz;
      DRT::ParObject::ExtractfromPack(posingroup,data,vertices_xyz);

      //store boundary integration cells in boundaryintcelllist
      intcellvector.push_back(GEO::BoundaryIntCell(distype, -1, vertices_xi,
                              Teuchos::null, vertices_xyz));
    }

    // add group of cells for this element to the map
    intcellmap.insert(make_pair(elegid,intcellvector));

    // check correct reading
    if (posingroup != (int)data.size())
      dserror("mismatch in size of data %d <-> %d",(int)data.size(),posingroup);
  }
  // check correct reading
  if (posofgroup != (int)dataRecv.size())
    dserror("mismatch in size of data %d <-> %d",(int)dataRecv.size(),posofgroup);
}

#endif // #ifdef CCADISCRET
