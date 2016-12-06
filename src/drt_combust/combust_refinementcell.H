/*!-----------------------------------------------------------------------------------------------*
 \file combust_refinementcell.H

 \brief

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef COMBUST_REFINEMENTCELL_
#define COMBUST_REFINEMENTCELL_


#include "../drt_lib/drt_element.H"


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

    //! constructor
    explicit RefinementCell(
        const DRT::Element* ele,
        const DRT::Element::DiscretizationType distype,
        const std::vector<std::vector<double> > vertexcoord
    );

    //! refine a cell (split it into subcells)
    void RefineCell();

    //! remote access method to root fluid element
    const DRT::Element* Ele() const { return ele_; }

    //! remote access method to distype
    DRT::Element::DiscretizationType Shape() const { return distype_; }

    //! remote access method to refinement level
    int RefinementLevel() const { return refinementlevel_; }

    //! remote access method to bisection state
    bool Bisected() const { return bisected_; }

    //! remote access method to touched state
    bool Touched() const { return touched_; }

    //! remote access method to G-function values; could be spared by declaring "gfuncvalues" protected/public
    void SetGfuncValues(std::vector<double> gfuncvalues);

    //! remote access method to G-function values; could be spared by declaring "gfuncvalues" protected/public
    const std::vector<double> & GetGfuncValues() const { return gfuncvalues_; }

    //! find out, whether this refinement cell is intersected by flame G-function or not
    void IdentifyIntersectionStatus();

    //! remote access method to coordinates of a refinement cell
    std::vector<std::vector<double> > GetVertexCoord() const { return vertexcoord_; }

    //! remote access method to number of children, depends on the tree type, here 8
    int NumOfChildren() const { return children_.size(); }

    //! remote access method to children
    const Teuchos::RCP<COMBUST::RefinementCell> GetRefinementCell (const int index) const;

    //! fills a vector with pointers to all refinement cells which are not further intersected
    void SearchRefinementCells (std::vector<const COMBUST::RefinementCell* >& Cells) const;

    //! remote access method to rootcell
    const COMBUST::RefinementCell* ReturnRootCell() const;

    //! deletes all children, i.e. all refinement cells of the root cell
    void Clear();

    //! create Gmsh output of rootcell
    void RootCellToGmsh() const;

    //! public: map of cell edge ID and intersection point coodinates (zero G-function)
    // remark: for hex8 cells a map would be sufficient; for hex20 cells a multimap is required
    std::multimap<int,std::vector<double> > intersectionpoints_;

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

    //! boolean indicating whether a refinement cell is bisected by the interface
    bool bisected_;

    //! boolean indicating whether a refinement cell is touched by the interface,
    //! at least one G-function value = 0
    bool touched_;

    //! pointer to parent node
    const COMBUST::RefinementCell* parent_;

    //! vector of pointers to the children nodes
    std::vector< Teuchos::RCP<COMBUST::RefinementCell> > children_;

  }; // class RefinementCell

} // namespace COMBUST

#endif // COMBUST_REFINEMENTCELL_H
