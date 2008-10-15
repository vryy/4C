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

  class DRT::Discretization;

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
            const Teuchos::RCP<DRT::Discretization>  xfemdis   ///< discretization that will be intersected
            );

    //! destructor
    virtual ~InterfaceHandle();

    //! return string representation of the interface handle
    std::string toString() const;

    //! bullshit comment: return string representation of the interface handle
    virtual void toGmsh(const int step) const;

    //! return pointer to the intersected discretization
    Teuchos::RCP<DRT::Discretization> xfemdis() const {return xfemdis_;};

    //! get pointer to map with all domain integration cells on this proc
    const std::map<int,GEO::DomainIntCells>* elementalDomainIntCells() const {return &elementalDomainIntCells_;};
    //! get pointer to map with all boundary integration cells on this proc
    const std::map<int,GEO::BoundaryIntCells>* elementalBoundaryIntCells() const {return &elementalBoundaryIntCells_;};

    //! get pointer to oct tree
    Teuchos::RCP<GEO::SearchTree> getSearchTree() const {return octTreenp_;};

    /// Print Information about the contained search tree
    void PrintTreeInformation(const int step) const
    {
        octTreenp_->evaluateTreeMetrics(step);
        return;
    }

    //! get pointer to elementsByLabel_
    const std::map<int,std::set<int> >* elementsByLabel() const {return &elementsByLabel_;};

    /*! tell me, whether the element is intersected or not
     * \return true, if element is intersected
     */
    bool ElementIntersected(
            const int element_gid   ///< unique global element id
            ) const;

    //! return list of domain integrationcells for a given global element id, returned list can be empty -> always terminates
    GEO::DomainIntCells GetDomainIntCells(
            const int                              gid,     ///< unique global element id
            const DRT::Element::DiscretizationType distype  ///< shape of thefinite element
            ) const;

    //! return list of boundary integrationcells for a given global element id, returned list can be empty -> always terminates
    GEO::BoundaryIntCells GetBoundaryIntCells(
            const int                              gid      ///< unique global element id
            ) const;

    //! get domain information at x_in at timestep n+1
    virtual int PositionWithinConditionNP(const BlitzVec3& x_in) const;

    //! get domain information at x_in at timestep n
    virtual int PositionWithinConditionN(const BlitzVec3& x_in) const;

    //! get domain information at x_in at timestep n+1
    virtual int PositionWithinConditionNP(const BlitzVec3&     x_in,
                                          GEO::NearestObject&  nearestobject) const;

    //! get domain information at x_in at timestep n
    virtual int PositionWithinConditionN(const BlitzVec3&     x_in,
                                         GEO::NearestObject&  nearestobject) const;

  protected:

    //! some checks to ensure properly filled integrationcell maps
    void SanityChecks() const;

    //! intersected domain
    const Teuchos::RCP<DRT::Discretization> xfemdis_;

    //! integration cells for domain (global id for element -> list of volume integration cells)
    std::map<int,GEO::DomainIntCells >      elementalDomainIntCells_;

    //! integration cells for boundary
    std::map<int,GEO::BoundaryIntCells >    elementalBoundaryIntCells_;

    //! for each xfem condition label there is a set of boundary elements
    std::map<int,std::set<int> >            elementsByLabel_;

    //! SearchTree octtree for timestep n+1
    const Teuchos::RCP<GEO::SearchTree>           octTreenp_;

    //! SearchTree octtree for timestep n
    const Teuchos::RCP<GEO::SearchTree>           octTreen_;

  }; // class InterfaceHandle
} // namespace XFEM

#endif  // #ifndef INTERFACE_H
#endif  // #ifdef CCADISCRET
