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

// extern struct _FILES  allfiles;

/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::FlameFront::FlameFront(
    const Teuchos::RCP<const DRT::Discretization> fluiddis,
    const Teuchos::RCP<DRT::Discretization> gfuncdis
    ) :
    fluiddis_(fluiddis),
    gfuncdis_(gfuncdis)
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

  // loop over fluid (combustion) row elements
  for (int iele=0; iele<fluiddis_->NumMyRowElements(); ++iele)
  {
    // get element from discretization
    const DRT::Element *ele = fluiddis_->lRowElement(iele);

//#ifdef DEBUG
    if(ele->Type() != DRT::Element::element_combust3)
      // this is not compulsory, but combust3 elements are expected here!
      dserror("unexpected element type: this should be of combust3 type!");
//#endif

    // create refinement cell from a fluid element -> cell will have same geometry as element!
    const Teuchos::RCP<COMBUST::RefinementCell> rootcell = rcp(new COMBUST::RefinementCell(ele));

    // refinement strategy is turned on
    if (Teuchos::getIntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REFINEMENT") == true)
      RefineFlameFront(rootcell);

    else // refinement strategy is turned off
      FindFlameFront(phinpcol,rootcell);

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

  //#ifdef DEBUG
        if (lm.size() != mygfuncvalues.size() )
          dserror("unexpected number of nodes: there should be 8 nodes per element!");
  //#endif

      //--------------------------------------------------------------------------------------------
      // store G-function values in refinement cell and decide whether this cell is intersected
      //--------------------------------------------------------------------------------------------
      // node numbering in element has to match data structure in cell
      cell->SetGfuncValues(mygfuncvalues);

//      std::cout << "G-Funktionswerte in Zelle verpackt für Element: " << cell->Ele()->Id() << " - geschnitten?: " << cell->Intersected() << endl;
//      for (unsigned i=0; i < mygfuncvalues.size(); i++)
//        std::cout << "Wert am Knoten " << i << ": " << mygfuncvalues[i] << endl;
    }
    //----------------------------------------------------------------------------------------------
    // higher level of refinement
    //----------------------------------------------------------------------------------------------
    else
      dserror("refinement not working yet!");
      /* hole GID dieses Fluid Elementes
           // const int gid = ele->Id();
         für alle Eckpunkte der Zelle
           werte G-Funktion an lokalen Koordinaten des Elementes GID aus
             // GID der Elemente beider Dis identisch!
             // lokale Fluid und G-Funktion Koordinaten identisch!
           speichere den Wert am Eckpunkt */
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
    // fertig (-> nächstes Fluidelement)
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
  std::vector<double> gfuncvalues = cell->GetGfuncValues();

  // fill vector with coordinates of intersection point
  std::vector<double> coordinates(3);
  for (int dim = 0; dim < 3; dim++)
  {
    coordinates[dim] = 0.0;
  }

  // fill map with 12 (stupid default: hex cell has 12 edges!) vectors holding coordinates
  std::map<int,std::vector<double> > intersectionpoints;
  for (int ipoint = 0; ipoint < 12; ipoint++)
  {
    intersectionpoints[ipoint] = coordinates;
  }

  // store intersection points in refinement cell
  cell->intersectionpoints_ = intersectionpoints;
  return;
}


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
     CreateIntegrationCells()
     TransformIntegrationCells()
  */

  // hier muss jetzt die map "flamefrontpatches_" mit Oberflächenstückchen gefüllt werden

  // eventuell macht es keinen Sinn die übergebene Zelle "cell" als doppelt "const" zu deklarieren

  flamefrontpatches_.insert(make_pair(cell->Ele()->Id(),cell));

  return;
}

/*------------------------------------------------------------------------------------------------*
 | triangulate the intersection surface of a refinement cell                          henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::TriangulateFlameFront() //input intersection points for a refinement cell
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 | build piecewise linear complex (PLC) in Tetgen format,based on Ursulas preparePLC() henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::buildPLC()
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 | calls the CDT to create burnt and unburnt integration cells                        henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::CreateIntegrationCells() // output: integration cells in Tetgen format
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 | this could easlily be integrated into CreateIntegrationCells()                     henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::TransformIntegrationCells() // output: integration cells in baci format
{
  return;
}

#endif // #ifdef CCADISCRET
