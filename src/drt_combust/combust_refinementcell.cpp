/*!-----------------------------------------------------------------------------------------------*
 \file combust_refinementcell.cpp

 \brief

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#include "combust_refinementcell.H"
#include "../drt_io/io_gmsh.H"


/*------------------------------------------------------------------------------------------------*
 | constructor: create refinement cell from element                                   henke 12/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::RefinementCell::RefinementCell(const DRT::Element* ele)
:
ele_(ele),
distype_(ele_->Shape()),
refinementlevel_(0),
vertexcoord_(DRT::UTILS::getNumberOfElementNodes(distype_),std::vector<double>(3,0.0)),
gfuncvalues_(DRT::UTILS::getNumberOfElementNodes(distype_),1.0),
bisected_(false),
touched_(false),
parent_(NULL)
{
  const int numnode = DRT::UTILS::getNumberOfElementNodes(distype_);

  // fill matrix holding coordinates of vertices of cell with local node coordinates of element
  switch (distype_)
  {
    case DRT::Element::hex8:
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    {
      for(int inode = 0; inode < numnode; ++inode)
      {
        for(int k = 0; k < 3; ++k) // problem dimension = 3!
        {
          vertexcoord_[inode][k] = DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[inode][k];
        }
      }
      break;
    }
    // more distypes could be added to switch; in the long term however, this class should be
    // templated with respect to the distype
    default:
      dserror("not supported for this element shape!");
  }

  children_.clear();

  return;
}

/*------------------------------------------------------------------------------------------------*
 | constructor: create refinement cell from refinement cell                           henke 12/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::RefinementCell::RefinementCell(const COMBUST::RefinementCell* cell, const std::vector<std::vector<double> > vertexcoord)
:
ele_(cell->ele_),
distype_(ele_->Shape()),
refinementlevel_(cell->RefinementLevel()+1),
vertexcoord_(vertexcoord),
gfuncvalues_(DRT::UTILS::getNumberOfElementNodes(distype_),1.0),
bisected_(false),
touched_(false),
parent_(cell)
{
//  std::cout << "Constructor Refinementcell" << std::endl;
  children_.clear();
  return;
}

/*------------------------------------------------------------------------------------------------*
 | constructor: create refinement cell from coordinates                               henke 11/11 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::RefinementCell::RefinementCell(
    const DRT::Element* ele,
    const DRT::Element::DiscretizationType distype,
    const std::vector<std::vector<double> > vertexcoord)
:
    ele_(ele),
    distype_(distype),
    refinementlevel_(0),
    vertexcoord_(vertexcoord),
    gfuncvalues_(DRT::UTILS::getNumberOfElementNodes(distype_),1.0),
    bisected_(false),
    touched_(false),
    parent_(NULL)
{
  //  std::cout << "Constructor Refinementcell" << std::endl;
  children_.clear();
  return;
}

/*------------------------------------------------------------------------------------------------*
 | set G-function values for refinement cell from outside the cell                    henke 03/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::RefinementCell::SetGfuncValues(std::vector<double> gfuncvalues)
{
  if (gfuncvalues.size() != gfuncvalues_.size())
    dserror("G-function value vectors don't have same length");

  gfuncvalues_ = gfuncvalues;
//#ifndef COMBUST_CUT
  // needed for refinement
  IdentifyIntersectionStatus();
//#endif
}

/*------------------------------------------------------------------------------------------------*
 | find out, whether refinement cell is intersected by the interface                  henke 03/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::RefinementCell::IdentifyIntersectionStatus()
{

//URSULA 270709
  /* When the scatra-field is initialized, nodes may be not exactly equal zero.
   * This can lead to incorrect intersection points and , hence, to a failure of the intersection
   * algorithm.
   * --------------------------------------------------------------------------------------------
   * Kann man das vielleicht schon bei der Initialisierung abfangen?
   * --------------------------------------------------------------------------------------------
   * Moreover, very small G-Function values can yield intersection points which are approximately
   * equal the vertices. This can lead small element parts as well as problems with TetGen. It remains
   * to find a definition of very small G-Function values which can then be reset to 0.
   * Preliminary,to avoid these problems, G-funcion values smaller than 10E-8 (?) are reset to zero.
   * ----------------------------------------------------------------------------------------------
   * G-Function < 10E-8 ist möglicherweise noch nicht optimal, erfüllt aber im Moment seinen Zweck
   * ----------------------------------------------------------------------------------------------
   */
  //TEST
  //std::cout << "G-Func" << std::endl;
  //for (std::size_t i=0; i<gfuncvalues_.size(); i++ )
  //  std::cout << gfuncvalues_[i] << std::endl;

  // schott Aug 6, 2010

  // REMARK:
  // boolian variables:
  // bisected_ : if we have at least one node with (+) and one node with (-)
  // touched_: if at least one node is 0.0
  //
  // here in refinementcell:
  // cell is touched, also if only touched at one node or edge
  // => boundary integrations cells are build in COMBUST::FlameFront::buildPLC only
  //    for really touched (touched at a whole face) cells

  // reset booleans, just in case they have been modified by mistake
  bisected_ = false;
  touched_ = false;

  // idea: Since the interface is defined by the zero iso-surface of the G-function, we look for
  // sign changes among the G-function values at the vertices of the refinement cell.
  unsigned counter = 0;
  bool zeros = false;

  // advance to first non-zero G-function value; stop one node before the last node
  while ((gfuncvalues_[counter] == 0.0) and (counter < (gfuncvalues_.size()-1)))
  {
    counter++;
    zeros = true; // if there are one or more zero values
  }
  // if there are values which are maybe not zero and have maybe the other sign
  if (counter < (gfuncvalues_.size()-1))
  {
  // first non-zero G-function value is negative
  if (gfuncvalues_[counter] < 0.0)
  {
    while (counter < (gfuncvalues_.size()-1) and
          (bisected_ == false)) // if we are bisected we can stop
    {
      counter++;
      // if next G-function value is positive
      if(gfuncvalues_[counter] >= 0.0)
      {
        if (gfuncvalues_[counter] > 0.0)
        {
          bisected_ = true;
        }
        else
          zeros = true;
      }
      // else -> G-function value is also negative
    }
  }
  // first non-zero G-function value is positive
  else if (gfuncvalues_[counter] > 0.0)
  {
    while (counter < (gfuncvalues_.size()-1) and
          (bisected_ == false)) // if we are bisected we can stop
   {
      counter++;
      // if next G-function value is negative
      if(gfuncvalues_[counter] <= 0.0)
      {
        if (gfuncvalues_[counter] < 0.0)
        {
          bisected_ = true;
        }
        else
          zeros = true;
      }
      // else -> G-function value is also positive
    }
  }
  else
  {
    std::cout << gfuncvalues_[counter] << std::endl;
    dserror("impossible!");
  }
  }
  // element is not bisected but interface touches the element
  if((zeros == true) and (bisected_ == false))
  {
    touched_ = true;
    //std::cout << "Interface berührt Element: " << ele_->Id() << std::endl;
  }
  if (bisected_ and touched_) dserror("impossible!");

  //std::cout << "Element " <<  ele_->Id() << " -- bisected status: " << bisected_ <<
  //                                  "\t" << " -- touched_ status: " << touched_  << std::endl;

  if (distype_ == DRT::Element::hex20)
  {
    // TODO complete the algortihm
    //dserror("read comment first!");
    //-------------------------------------------------------------------------------------------
    // Im urspruenglichen, ersten Teil werden nur die die Vorzeichen der Konten des Elements
    // auf einen Vorzeichenwechsel hin ueberprueft. Das ist fuer Hex8 absolut ok, reicht aber fuer
    // Hex20 nicht aus. Diese Kiterium reicht dann nicht, wenn alle Koten das gleiche Vorzeichen
    // haben. Fuer diesen Fall ist es moegchlich, dass zwei Nullstellen zwischen einem der
    // aeusseren Knoten und dem Mittelknoten liegen. Nur eine Nullstelle ist auch moeglich und
    // entspricht dann einem Beruehren.
    //
    // Dazu wurde von einem Studenten ein zusaetzliches Kriterium implementiert. Das, was der
    // Student hier gemacht hat, muss noch uberprueft werden.
    //
    // In den meisten Schnittfaellen (ich glaube in allen mit einem Interface je Element) sollte
    // die obige Ueberpruefung des IntersectionStatus vollig ausreichend sein. Daher habe ich diesen
    // Teil mal raus genommen.
    //
    // TODO: @Kilian: Ich glaube, du kannst das erstmal weglassen.
    //-------------------------------------------------------------------------------------------
# if 0
    //changed_grundl Funktion teilweise geändert (weiterhin keine Garantie ...) (beim Kopieren darauf achten den alten Code stehen zu lassen)
    if ((distype_ == DRT::Element::hex20) and (bisected_ == false))
    {
      //  dserror("This part has been done by a student and should be check before use!");
      // FARAH
      // Grundl 11/2010

      std::vector<std::vector<int> > lines = DRT::UTILS::getEleNodeNumberingLines(DRT::Element::hex20);

      int ecounter = 2;  //Grad des Polynoms

      for (std::size_t iline=0; iline < lines.size(); iline++)
      {

        //the function is f(x)=a1*xi^2+a2*xi+a3;
        //lines[iline][0] is the first value (xi = -1)
        //lines[iline][1] is the last value (xi = +1)
        //lines[iline][2] is the value in the middle (xi = 0)
        double a1 = (gfuncvalues_[lines [iline][1]]+gfuncvalues_[lines [iline][0]])/2 - gfuncvalues_[lines [iline][2]];
        double a2 = (gfuncvalues_[lines [iline][1]]-gfuncvalues_[lines [iline][0]])/2;
        double a3 = gfuncvalues_[lines [iline][2]];

        //        std::cout << "G-Functionvalues for identification of intersection-status:" <<  std::endl;
        //        std::cout << "First-Point:  " << gfuncvalues_[lines [iline][0]]  << std::endl;
        //        std::cout << "Middle-Point:  " << gfuncvalues_[lines [iline][2]]  << std::endl;
        //        std::cout << "End-Point:  " << gfuncvalues_[lines [iline][1]]  << std::endl;

        //Diskriminante
        double D = a2*a2-4*a1*a3;

        int p = 0;  //counter for positive values
        int n = 0;  //counter for negative values
        for (int b = 0; b <= ecounter; b++)
        {
          if ( gfuncvalues_[lines[iline][b]] >0)
          {
            p++;
          }
          if ( gfuncvalues_[lines[iline][b]] < 0)
          {
            n++;
          }
        }

        if ( p == 3 or n == 3 ) //three values with the same sign on one line
        {
          if (D >= 0)
          {
            double x1= (-a2+ sqrt(D))/(2*a1);
            double x2= (-a2-sqrt(D))/(2*a1);

            if ( abs(x1) < 1)
              if (x1 != x2 )
              {
                bisected_ = true;
              }
              else //double zero
              {
                touched_ = true;
              }
            else if (abs(x2) < 1)
            {
              bisected_ = true;
            }
          }
        }

      }

      // 1 Wert != 0 und alle anderen 0
      //Kann das überhaupt vorkommen?

      if (bisected_ == false)
      {

        int counter2= 0;

        for(int counter1 = 0; counter1 < gfuncvalues_.size()-1; counter1++)
        {
          if(gfuncvalues_[counter1] != 0.0)
          {
            counter2++;
          }
        }

        if(counter2==1)
        {
          bisected_ = true;
          std::cout << "please check code here in:" << std::endl;
          std::cout << "COMBUST::RefinementCell::IdentifyIntersectionStatus()" << std::endl;
        }
      }
    }
#endif

#if 0
    double D = 0;
    int ecounter = 2;

    //  switch (cell->Ele()->Shape())     //???
    //  {
    //  case DRT::Element::hex8:
    //        {ecounter = 1;}
    //
    //  case DRT::Element::hex20:
    //        {ecounter = 2;}
    //
    //  default:
    //      dserror("IdentifyIntersectionStatus() does not support this element shape!");

    // }

    // 3 positive Werte auf kante
    if ((intersected_ == false) and (distype_==DRT::Element::hex20))
    {


      for (std::size_t iline=0; iline < lines.size(); iline++)
      {  int c = 0;


      for (int b = 0; b <= ecounter; b++)
      {


        if ( gfuncvalues_[lines[iline][b]] >0)
        {
          c++;

        }
        //else
        // break;
      }

      if ( c == 3 )
        // Diskriminante
      {//D=(((5*gfuncvalues_[lines [iline][1]]-3*gfuncvalues_[lines [iline][0]]-2*gfuncvalues_[lines [iline][2]])(5*gfuncvalues_[lines [iline][1]]-3*gfuncvalues_[lines [iline][0]]-2*gfuncvalues_[lines [iline][2]]))-4*gfuncvalues_[lines [iline][0]]*(0.5*gfuncvalues_[lines [iline][2]]-gfuncvalues_[lines [iline][1]]+0.5*gfuncvalues_[lines [iline][0]]));

        double a1= 0.5*gfuncvalues_[lines [iline][2]]-gfuncvalues_[lines [iline][1]]+0.5*gfuncvalues_[lines [iline][0]];
        double a2= 2*gfuncvalues_[lines [iline][1]]-1.5*gfuncvalues_[lines [iline][0]]-0.5*gfuncvalues_[lines [iline][2]];
        double a3= gfuncvalues_[lines [iline][0]];
        D= a2*a2-4*a1*a3;
        if (D!=0)
        {
          double x1= (-a2+ sqrt(a2*a2-4*a1*a3))/(2*a1);
          double x2= (-a2-sqrt(a2*a2-4*a1*a3))/(2*a1);


          if ((D==0)and (0<=x1)and (2>=x1))
          {touched_ = true;
          //break;
          }
          if ((0 < D)and (((0<=x1)and (2>=x1))or((0<=x2)and (2>=x2))))
          {intersected_ = true;
          //break;
          }
        }
      }

      }
    }
    // 3 negative Werte auf Kante
    if ((intersected_ == false) and (distype_==DRT::Element::hex20))
    {

      for (std::size_t iline=0; iline < lines.size(); iline++)
      {
        int c = 0;

        for (int b = 0; b <= ecounter; b++)
        {
          if (0 > gfuncvalues_[lines [iline][b]])
          {
            c++;
          }
        }
        //else
        //break;
        //};

        if ( c == ecounter )
          //Diskriminante
        {//D=(((5*gfuncvalues_[lines [iline][1]]-3*gfuncvalues_[lines [iline][0]]-2*gfuncvalues_[lines [iline][2]])(5*gfuncvalues_[lines [iline][1]]-3*gfuncvalues_[lines [iline][0]]-2*gfuncvalues_[lines [iline][2]]))-4*gfuncvalues_[lines [iline][0]]*(0.5*gfuncvalues_[lines [iline][2]]-gfuncvalues_[lines [iline][1]]+0.5*gfuncvalues_[lines [iline][0]]));
          double a1= 0.5*gfuncvalues_[lines [iline][2]]-gfuncvalues_[lines [iline][1]]+0.5*gfuncvalues_[lines [iline][0]];
          double a2= 2*gfuncvalues_[lines [iline][1]]-1.5*gfuncvalues_[lines [iline][0]]-0.5*gfuncvalues_[lines [iline][2]];
          double a3= gfuncvalues_[lines [iline][0]];
          D= a2*a2-4*a1*a3;
          if(D!=0)
          {
            double x1= (-a2+ sqrt(a2*a2-4*a1*a3))/(2*a1);
            double x2= (-a2-sqrt(a2*a2-4*a1*a3))/(2*a1);


            if ((D==0)and (0<=x1)and (2>=x1))
            {
              touched_ = true;
            }
            if ((0 < D)and (((0<=x1)and (2>=x1))or((0<=x2)and (2>=x2))))
            {
              intersected_ = true;
            }
          }
        }
      }
    }

    // 1 Wert != 0 und alle anderen 0

    if (intersected_ == false)
    {

      int counter2= 0;

      for(int counter1 = 0; counter1 < gfuncvalues_.size()-1; counter1++)
      {
        if(gfuncvalues_[counter1] != 0.0)
        {
          counter2++;
        }
      }

      if(counter2==1)
      {
        intersected_ = true;
      }
    }

    // Prüfe ob Nullstelle(n) auf kante liegen

    //double a1= 0.5*gfuncvalues_[lines [iline][2]]-gfuncvalues_[lines [iline][1]]+0.5*gfuncvalues_[lines [iline][0]];
    //double a2= 5*gfuncvalues_[lines [iline][1]]-3*gfuncvalues_[lines [iline][0]]-2*gfuncvalues_[lines [iline][2]];
    //double a3= gfuncvalues_[lines [iline][0]];
    //
    //double x1= (-a2+(a2*a2-4*a1*a3)^0,5)/(2*a1);
    //double x2= (-a2-(a2*a2-4*a1*a3)^0,5)/(2*a1);


    //FARAH
#endif
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
 | refine a cell embedded in the octree at the level above                        rasthofer 08/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::RefinementCell::RefineCell() // input: a refinement cell; output: geteilte Zellen
{
/*
  Die Funktionalität dieser Funktion (teile eine Zelle und gib die Unterzellen zurück) entspricht
  der Klasse "treenode" wie sie in Ursulas Searchtree existiert. Die diese Funktion sollte also
  irgendwann als member function der Klasse RefinementCell() implementiert werden.       henke 12/08
*/
//  std::cout<<"RefineCell"<<std::endl;
  std::vector<std::vector<double> > refinevertexcoord;
  const int numnode = ele_->NumNode();

/* Var1:
  std::vector<double> actrefinevertexcoord(3);// problem dimension = 3
  const int numvertex = 28;
  for (int inode = 0; inode<numnode; inode++)
    refinevertexcoord.push_back(vertexcoord_[inode]);
  for(int ivertex = numnode; ivertex < numvertex; ++ivertex)
  {
 // get vertex coordinate of refinement cell in the coordinate system of the refinement cell
 for (int k = 0; k < 3; ++k)
    {
      actrefinevertexcoord[k] = DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[ivertex][k];
    }
 // evaluate shape function at the coordinates of the current refinement cell vertex
 Epetra_SerialDenseVector funct(numnode);
 DRT::UTILS::shape_function_3D(funct,actrefinevertexcoord[0],actrefinevertexcoord[1],actrefinevertexcoord[2],DRT::Element::hex8);
    // compute the coordinates of the actual refinement cell vertex in element coordinates
 std::vector<double> actvertexcoord(3);
    for (int k = 0; k < 3; ++k)
    {
     actvertexcoord[k] = 0.0;
        for (int inode = 0; inode < numnode; inode++)
        {
          actvertexcoord[k] += vertexcoord_[inode][k] * funct(inode);
        }
    }
    refinevertexcoord.push_back(actvertexcoord);
  }
*/

  // first all vertices of the refinement cells are computed from the vertices of
  // the cell which is refined
  // these vertices are stored in refinevertexcoord
  // therefore, the midpoint of each line, surface and the cell itself is computed

  // assignment of vertices to lines -> compute midpoint of line
  std::vector<std::vector<int> > lines = DRT::UTILS::getEleNodeNumberingLines(distype_);
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
  // assignment of vertices to surface -> compute midpoint of surface
  int surf[6][2] = { { 8, 10},
                     { 8, 16},
                     { 9, 17},
                     {10, 18},
                     {11, 19},
                     {19, 17} };
  // assignment of vertices to cell -> compute midpoint of cell
  int center[2]= {20, 25};

  std::vector<double> actvertexcoord(3);

  // store the vertices of cell which is refined
  for (int inode = 0; inode < numnode; inode++)
    refinevertexcoord.push_back(vertexcoord_[inode]);

  //compute and store midpoint of lines
  for (std::size_t iline = 0; iline < lines.size(); iline++)
  {
    for (int k = 0; k < 3; ++k)
      actvertexcoord[k] = 0.5 * (refinevertexcoord[lines[iline][0]][k] + refinevertexcoord[lines[iline][1]][k]);
    refinevertexcoord.push_back(actvertexcoord);
  }

  //compute and store midpoint of surfaces
  for (int isurf = 0; isurf < 6; isurf++)
  {
    for (int k = 0; k < 3; ++k)
      actvertexcoord[k] = 0.5 * (refinevertexcoord[surf[isurf][0]][k] + refinevertexcoord[surf[isurf][1]][k]);
    refinevertexcoord.push_back(actvertexcoord);
  }

  //compute and store midpoint of cell which is refined
  for (int k = 0; k < 3; ++k)
    actvertexcoord[k] = 0.5 * (refinevertexcoord[center[0]][k] + refinevertexcoord[center[1]][k]);
  refinevertexcoord.push_back(actvertexcoord);

  // assigment of the vertices in refinevertexcoord to the new refinement cells
  int verticesorder[8][8] = { { 0,  8, 20, 11, 12, 21, 26, 24},
                              { 8,  1,  9, 20, 21, 13, 22, 26},
                              {20,  9,  2, 10, 26, 22, 14, 23},
                              {11, 20, 10,  3, 24, 26, 23, 15},
                              {12, 21, 26, 24,  4, 16, 25, 19},
                              {21, 13, 22, 26, 16,  5, 17, 25},
                              {26, 22, 14, 23, 25, 17,  6, 18},
                              {24, 26, 23, 15, 19, 25, 18,  7} };

  // loop over all new cells
  for (int icell = 0; icell < 8; icell++)
  {
    std::vector<std::vector<double> > refinecellvertexcoord;
    for (int icellvertex = 0; icellvertex < 8; icellvertex++)
    {
      // select corresponding vertices from verticesorder
          refinecellvertexcoord.push_back(refinevertexcoord[verticesorder[icell][icellvertex]]);
//          //TEST
//          std::cout << "Cell " << icell << std::endl;
//          for (int k=0; k<3; k++)
//          {
//            std::cout << refinecellvertexcoord[icellvertex][k] << std::endl;
//          }
    }

    // create new refinement cell
    children_.push_back(Teuchos::rcp(new RefinementCell(this, refinecellvertexcoord)));
//    std::cout << "Cell  build" << std::endl;
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | returns the refinement cells which are not further intersected,                                |
 | i.e. the leafs of the search tree                                              rasthofer 08/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::RefinementCell::SearchRefinementCells (std::vector<const COMBUST::RefinementCell* >& Cells)const
{
  if (children_.size()>0)// cell is refined
  {
    // got to the next level
    for (std::size_t ichild=0; ichild < children_.size(); ichild++)
      children_[ichild]->SearchRefinementCells (Cells);
  }
  else // cell is not refined -> leaf of tree -> store
    Cells.push_back(this);
  return;
}

/*------------------------------------------------------------------------------------------------*
 | returns the root cell, i.e. the root of the search tree                        rasthofer 08/09 |
 *------------------------------------------------------------------------------------------------*/
const COMBUST::RefinementCell* COMBUST::RefinementCell::ReturnRootCell() const
{
  const COMBUST::RefinementCell* tmp = NULL;
  if (parent_ == NULL) // cell is root
    tmp = this;
  else // cell is not root -> go back one level
    tmp = parent_->ReturnRootCell();
  return tmp;
}

/*------------------------------------------------------------------------------------------------*
 | deletes all children, i.e. all refinement cells of the root cell               rasthofer 05/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::RefinementCell::Clear()
{
  if (children_.size()>0)
  {
    for (size_t i=0; i<children_.size(); i++)
    {
      if (children_[i]->NumOfChildren()>0)
      {
        children_[i]->Clear();
      }
    }
//    std::cout << children_.size() << " Cells deleted" << std::endl;
    children_.clear();
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
 | returns the refinement cell of refined cell                                    rasthofer 08/09 |
 *---------------------------------------------------------------------- -------------------------*/
const Teuchos::RCP<COMBUST::RefinementCell> COMBUST::RefinementCell::GetRefinementCell (const int index) const
{
  return children_[index];
}

void COMBUST::RefinementCell::RootCellToGmsh() const
{
  if (parent_ != NULL)
     dserror("Gmsh output only for rootcells!");

  const bool screen_out = false;

  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("Rootcell_and_Intersectionpoints", 999, 5, screen_out, 0);
  std::ofstream gmshfilecontent(filename.c_str());
  {
    gmshfilecontent << "View \" " << "Rootcell \" {\n";
    int Id = ele_->Id();
    int numnode = ele_->NumNode();
    LINALG::SerialDenseMatrix Pos(3,numnode);
    for (int i=0; i<numnode; i++)
    {
      for (size_t k=0; k<3; k++)
      {
        Pos(k,i) = vertexcoord_[i][k];
      }
    }
    IO::GMSH::cellWithScalarToStream(distype_, Id, Pos, gmshfilecontent);
    gmshfilecontent << "};\n";
  }
  {
    gmshfilecontent << "View \" " << "Intersectionpoints \" {\n";
    for(std::map<int,std::vector<double> >::const_iterator iter = intersectionpoints_.begin(); iter != intersectionpoints_.end(); ++iter)
    {
      int id = iter->first;
      const std::vector<double> point = iter->second;
      LINALG::Matrix<3,1> pos;
      for (size_t k=0; k<point.size(); k++)
         pos(k) = point[k];

      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, id, pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }

  return;
}

