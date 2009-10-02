/*!-----------------------------------------------------------------------------------------------*
 \file combust_refinementcell.cpp

 \brief

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>

 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "combust_refinementcell.H"

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
intersected_(false),
parent_(NULL)
{
  const int numnode = DRT::UTILS::getNumberOfElementNodes(distype_);

  // fill matrix holding coordinates of vertices of cell with local node coordinates of element
  switch (distype_)
  {
    case DRT::Element::hex8: // hex20, hex27 genauso!
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
intersected_(false),
parent_(cell)
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

  IdentifyIntersectionStatus();
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
//  std::cout << "G-Func" << std::endl;
//  for (std::size_t i=0; i<gfuncvalues_.size(); i++ )
//	  std::cout << gfuncvalues_[i] << std::endl;

  for (std::size_t i=0; i<gfuncvalues_.size(); i++ )
  {
    if (fabs(gfuncvalues_[i])<1.0E-7)
    {
      gfuncvalues_[i] = 0.0;
//      std::cout << " G-Function value  reset to 0 " << std::endl;
    }
  }

  // idea: Since the interface is defined by the zero iso-surface of the G-function, we look for
  // sign changes among the G-function values at the vertices of the refinement cell.
  unsigned counter = 0;

  // advance to first non-zero G-function value
  while ((gfuncvalues_[counter] == 0.0) and (counter < (gfuncvalues_.size()-1)))
  {
    counter++;
    intersected_ = true;
  }

  // first non-zero G-function value is negative
  if (gfuncvalues_[counter] < 0.0)
  {
    while (counter < (gfuncvalues_.size()-1) and (intersected_ == false))
    {
      counter++;
      // if next G-function value is positive
      if(gfuncvalues_[counter] >= 0.0)
        intersected_ = true;
    }
  }
  // first non-zero G-function value is positive
  else if (gfuncvalues_[counter] > 0.0)
  {
    while (counter < (gfuncvalues_.size()-1) and (intersected_ == false))
   {
      counter++;
      // if next G-function value is negative
      if(gfuncvalues_[counter] <= 0.0)
        intersected_ = true;
    }
  }
  else
    dserror("impossible!");
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
    children_.push_back(rcp(new RefinementCell(this, refinecellvertexcoord)));
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
 | returns the refinement cell of refined cell                                    rasthofer 08/09 |
 *---------------------------------------------------------------------- -------------------------*/
const Teuchos::RCP<COMBUST::RefinementCell> COMBUST::RefinementCell::GetRefinementCell (const int index) const
{
  return children_[index];
}


#endif // #ifdef CCADISCRET
