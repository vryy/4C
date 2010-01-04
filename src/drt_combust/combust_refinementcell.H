/*!-----------------------------------------------------------------------------------------------*
 \file combust_refinementcell.H

 \brief

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef COMBUST_REFINEMENTCELL_
#define COMBUST_REFINEMENTCELL_

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
//#include "../drt_lib/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

namespace COMBUST
{
  /*!
    class for a refinement cell

   Bemerkung:
   Die Klasse RefinementCell() sollte irgendwann per Vererbung von Treenode aus dem Searchtree
   abgeleitet werden.

    \author henke
    \date 12/08
   */
  class RefinementCell
  {

  public:

    //! constructor
    explicit RefinementCell(const DRT::Element *ele);

    //! constructor
    explicit RefinementCell(const COMBUST::RefinementCell *cell, const std::vector<std::vector<double> > vertexcoord);

    //! refine a cell (split it into subcells)
    void RefineCell();

    //! remote access method to root fluid element
    const DRT::Element* Ele() const { return ele_; }

    //! remote access method to refinement level
    const int RefinementLevel() const { return refinementlevel_; }

    //! remote access method to intersection state
    bool Intersected() const { return intersected_; }

    //! remote access method to G-function values; could be spared by declaring "gfuncvalues" protected/public
    void SetGfuncValues(std::vector<double> gfuncvalues);

    //! remote access method to G-function values; could be spared by declaring "gfuncvalues" protected/public
    std::vector<double>  GetGfuncValues() const { return gfuncvalues_; }

    //! find out, whether this refinement cell is intersected by flame G-function or not
    void IdentifyIntersectionStatus();

    //! remote access method to coordinates of a refinement cell
    std::vector<std::vector<double> > GetVertexCoord() const { return vertexcoord_; }

    //! remote access method to number of children, depends on the tree type, here 8
    const int NumOfChildren() const { return children_.size(); }

    //! remote access method to children
    const Teuchos::RCP<COMBUST::RefinementCell> GetRefinementCell (const int index) const;

    //! fills a vector with pointers to all refinement cells which are not further intersected
    void SearchRefinementCells (std::vector<const COMBUST::RefinementCell* >& Cells) const;

    //! remote access method to rootcell
    const COMBUST::RefinementCell* ReturnRootCell() const;

    //! public: map of cell edge ID and intersection point coodinates (zero G-function)
    std::map<int,std::vector<double> > intersectionpoints_;

  private:

    //! pointer to the fluid (combustion) element this cell belongs to
    const DRT::Element* ele_;

    //! geometrical shape of the refinement cell corresponds to shape of element
    const DRT::Element::DiscretizationType distype_;

    //! counter for level of refinement - corresponds to number of refinement steps done to get this cell
    const int refinementlevel_;

    //! matrix of vertex coordinates defining this refinementcell (local fluid element coordinates!)
    std::vector<std::vector<double> > vertexcoord_;

    //! vector of G-function values at vertices of this refinementcell
    std::vector<double> gfuncvalues_;

    //! boolean indicating whether a refinement cell is intersected by the interface
    bool intersected_;

    //! pointer to parent node
    const COMBUST::RefinementCell* parent_;

    //! vector of pointers to the children nodes
    std::vector< Teuchos::RCP<COMBUST::RefinementCell> > children_;

  }; // class RefinementCell
} // namespace COMBUST

#endif // COMBUST_REFINEMENTCELL_H
#endif // CCADISCRET
