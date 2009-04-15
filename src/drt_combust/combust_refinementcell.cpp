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
intersected_(false)
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

  return;
}

/*------------------------------------------------------------------------------------------------*
 | constructor: create refinement cell from refinement cell                           henke 12/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::RefinementCell::RefinementCell(const COMBUST::RefinementCell* cell)
:
ele_(cell->ele_),
distype_(ele_->Shape()),
refinementlevel_(cell->RefinementLevel()+1),
//vertexcoord_(?),
//gfuncvalues_(?)
intersected_(false)
{
  dserror("Thou shalt not call this constructor, cause it ain't workin!");
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
 | find out, whether refinement cell is intersected by the zero iso-surface of the G-function by  |
 | looking for sign changes among the G-function values at the vertices               henke 03/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::RefinementCell::IdentifyIntersectionStatus()
{
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
 | refine a cell embedded in the octree at the level above                            henke 12/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::RefinementCell::RefineCell() // input: a refinement cell; output: geteilte Zellen
{
/*
  Die Funktionalität dieser Funktion (teile eine Zelle und gib die Unterzellen zurück) entspricht
  der Klasse "treenode" wie sie in Ursulas Searchtree existiert. Die diese Funktion sollte also
  irgendwann als member function der Klasse RefinementCell() implementiert werden.       henke 12/08
*/
  return;
}
