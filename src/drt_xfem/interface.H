/*!
 \file interface.H

 \brief interface handle that transports the intersection related things around

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */

#ifdef CCADISCRET
#ifndef INTERFACE_H
#define INTERFACE_H

#include "../drt_geometry/integrationcell.H"
#include "../drt_geometry/searchtree.H"

namespace XFEM
{

  /*!
 * It knows everything about the XFEM interface geometry
 *
 * It does not deal with DOFs in any way. DOF management is done by
 * the DofManager and ElementDofManager classes.
 *
 */
  class InterfaceHandle
  {
  public:

    //! constructor
    explicit InterfaceHandle(
            const Teuchos::RCP<DRT::Discretization>&  xfemdis   ///< discretization that will be intersected
            );

    //! destructor
    virtual ~InterfaceHandle();

    //! return string representation of the interface handle
    virtual std::string toString() const;

    //! output interface handle to Gmsh postprocesing
    virtual void toGmsh(const int step) const;

    //! return pointer to the intersected discretization
    Teuchos::RCP<DRT::Discretization> xfemdis() const {return xfemdis_;}

    //! get pointer to oct tree
    Teuchos::RCP<GEO::SearchTree> getSearchTree() const {return octTreenp_;}

    /// Print Information about the contained search tree
    void PrintTreeInformation(const int step) const
    {
        octTreenp_->evaluateTreeMetrics(step);
        return;
    }

//    //! get pointer to elementsByLabel_
//    const std::map<int,std::set<int> >& BoundaryElementsByLabel() const {return boundaryElementsByLabel_;}

    //! get set of all available labels
    std::set<int> GetAvailableBoundaryLabels() const;

    //! tell me, whether the element is intersected or not
    bool ElementIntersected(const int element_gid) const;

    //! check, whether given element id has been intersected by the given interface (label)
    bool ElementHasLabel(
        const int element_gid,
        const int label
    ) const;

    //! get set of boundary element GIDs that intersect an element
    std::set<int> GetIntersectingBoundaryElementsGID(
        const int element_gid  ///< element GID of intersected element
        ) const;

    //! get labels for each element
    std::set<int> LabelsPerElement( const int element_gid) const;

    //! return list of domain integrationcells for a given global element id, returned list can be empty -> always terminates
    GEO::DomainIntCells GetDomainIntCells( const DRT::Element*  xfemElement) const;

    //! return number of domain integration cells for a given element
    std::size_t GetNumDomainIntCells( const DRT::Element*  xfemElement) const;

    //! return map of elemental boundary integration cells
    std::map<int, GEO::BoundaryIntCells> GetElementalBoundaryIntCells () const { return elementalBoundaryIntCells_; };

    //! return list of boundary integrationcells for a given global element id, returned list can be empty -> always terminates
    GEO::BoundaryIntCells GetBoundaryIntCells(const int   gid) const;

    //! get domain information at x_in at timestep n+1
    virtual int PositionWithinConditionNP(const LINALG::Matrix<3,1>& x_in) const;

    //! get domain information at x_in at timestep n
    virtual int PositionWithinConditionN(const LINALG::Matrix<3,1>& x_in) const;

    //! get domain information at x_in at timestep n+1
    virtual int PositionWithinConditionNP(const LINALG::Matrix<3,1>&     x_in,
                                          GEO::NearestObject&  nearestobject) const;

    //! get domain information at x_in at timestep n
    virtual int PositionWithinConditionN(const LINALG::Matrix<3,1>&     x_in,
                                         GEO::NearestObject&  nearestobject) const;

    //! get the label to which the boundary element belongs
    virtual int GetLabelPerBoundaryElementId(const int boundaryEleId) const
    {
      return labelPerBoundaryElementId_.find(boundaryEleId)->second;
    }

  protected:

    //! fill label per element list
    void InvertElementsPerLabel();

    //! intersected domain
    const Teuchos::RCP<DRT::Discretization>   xfemdis_;

    //! integration cells for domain (global id for element -> list of volume integration cells)
    std::map<int,GEO::DomainIntCells >        elementalDomainIntCells_;

    //! integration cells for boundary
    std::map<int,GEO::BoundaryIntCells >      elementalBoundaryIntCells_;

    //! for each xfem condition label there is a set of boundary elements
    std::map<int,std::set<int> >              boundaryElementsByLabel_;

     //! each boundary element belongs to one label (inverted information of elementsByLabel_)
    std::map<int,int>                         labelPerBoundaryElementId_;

    //! SearchTree octtree for timestep n+1
    const Teuchos::RCP<GEO::SearchTree>       octTreenp_;

    //! SearchTree octtree for timestep n
    const Teuchos::RCP<GEO::SearchTree>       octTreen_;


  private:

    //! don't want copy constructor
    explicit InterfaceHandle(const InterfaceHandle&);

    //! disabled assignment operator
    InterfaceHandle operator = (const InterfaceHandle& old);

  }; // class InterfaceHandle

} // namespace XFEM

#endif  // #ifndef INTERFACE_H
#endif  // #ifdef CCADISCRET
