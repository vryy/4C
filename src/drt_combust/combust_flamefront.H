/*!-----------------------------------------------------------------------------------------------*
 \file combust_flamefront.H

 \brief

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef COMBUST_FLAMEFRONT_H
#define COMBUST_FLAMEFRONT_H

#include "../drt_lib/drt_globalproblem.H"
#include "combust_refinementcell.H"

namespace COMBUST
{

  class DRT::Discretization;

  /// COMBUST FlameFront class
  /*!

    The purpose of this class is to perform the intersection of the fluid domain by the flame front
    and to provide integration cells (in the burnt and unburnt domain) for the fluid time integration
    scheme. A refinement strategy is used to increase accuracy and simplify intersection patterns.

    The following tasks are peformed:

    refinement
    - initialize the searchtree
    - start a recursive loop to refine the region around the flame front until the maximal number
      of refinements is reached
      {
        - identify the intersection points (zeros of the level set function (G-equation) on the
          refinement cell edges)
        - if a refinement cell is intersected, it will be refined once more
      }

    capture flame front
    - triangulate flame front surface (spanned by the intersection points) in every refinement cell
    - create a piecewise linear complex (PLC) in Tetgen format from every refinement cell
    - call the Constraint Delaunay Tetrahedrization (CDT) to produce burnt and unburnt integration cells
    - transform integration cells from Tetgen format to baci format

    \author henke
    \date 10/08
   */
  class FlameFront
  {
  public:
    /// constructor
    explicit FlameFront(
             const Teuchos::RCP<const DRT::Discretization>  fluiddis, ///< fluid discretization (combust3 elements) which will be intersected
             Teuchos::RCP<DRT::Discretization>  gfuncdis  ///< G-function discretization (condif3 elements) describing the level set function
             );
    /// destructor
    ~FlameFront();

    /// remote access method to flame front patches
    const std::map< int,const Teuchos::RCP<const COMBUST::RefinementCell> > FlameFrontPatches() const
        { return flamefrontpatches_; }

    /// flame front control routine
    void ProcessFlameFront(
           const Teuchos::ParameterList& combustdyn,
           const Teuchos::RCP<Epetra_Vector> phinp);

  private:

    /// don't want default constructor
    explicit FlameFront();

    /// don't want copy constructor
    explicit FlameFront(const FlameFront&);

    /// refine the region around flame front
    void RefineFlameFront(const Teuchos::RCP<const COMBUST::RefinementCell> cell);

    /// find the flame front within a refinement cell according to G-function field
    void FindFlameFront(
           const Teuchos::RCP<Epetra_Vector>           phi,
           const Teuchos::RCP<COMBUST::RefinementCell> cell);

    /// find intersection points of G-function (level set zero iso-surface) with refinement cell edges
    void FindIntersectionPoints(const Teuchos::RCP<COMBUST::RefinementCell> cell);

    /// capture flame front within one refinement cell to provide integration cells
    void CaptureFlameFront(const Teuchos::RCP<const COMBUST::RefinementCell> cell);

    /// not documented yet
    void TriangulateFlameFront();

    // these 3 should be moved to a class accepting the 2D flame front as an input and producing 3D
    // integration cells as an output
    /// not documented yet
    void buildPLC();
    /// not documented yet
    void CreateIntegrationCells();
    /// not documented yet
    void TransformIntegrationCells();

    /// pointer to fluid discretization
    const Teuchos::RCP<const DRT::Discretization> fluiddis_;

    /// pointer to level set (G-function) discretization that intersects the fluiddis_
    const Teuchos::RCP<DRT::Discretization> gfuncdis_;

    /// vector on fluid NodeColMap holding G-function values originally stored in phinp
    const Teuchos::RCP<Epetra_Vector> phinp_;

    /// map containing for each fluid element flame front geometry to be exported to the interfacehandle
    // this is preliminary, just to have a "flame front" for now: it should hold the interface information
    std::map< int,const Teuchos::RCP<const COMBUST::RefinementCell> > flamefrontpatches_;
//    std::map< int, PLC_per_element > flamefrontpatches_;

  }; // class FlameFront
} // namespace COMBUST

#endif // COMBUST_FLAMEFRONT_H
#endif // CCADISCRET
