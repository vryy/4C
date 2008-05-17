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

#include <vector>
#include <set>
#include <string>

#include "integrationcell.H"
#include "../drt_lib/drt_discret.H"
#include "intersection.H"
#include "../drt_io/io_gmsh.H"

namespace XFEM
{

/*!
 \brief It knows everything about the XFEM enrichments for one particular interface setup

 */
class InterfaceHandle
{
public:

    //! constructor
    explicit InterfaceHandle(
            const RCP<DRT::Discretization>        xfemdis,      ///< discretization that will be intersected
            const RCP<DRT::Discretization>        cutterdis,    ///< discretization with boundary elements
            const RCP<Epetra_Vector>              idispcol      ///< interface displacements for current position of interface
            );
    //! destructor
    ~InterfaceHandle();

    //! return string representation of the interface handle
    std::string toString() const;
    
    //! return string representation of the interface handle
    void toGmsh(const int step) const;

    //! return pointer to the intersected discretization
    RCP<DRT::Discretization> xfemdis() const {return xfemdis_;};
    //! return pointer to the intersecting discretization
    RCP<DRT::Discretization> cutterdis() const {return cutterdis_;};
    
    const std::map<int,BlitzVec3>* currentcutterpositions() const {return &currentcutterpositions_;};

    //! get pointer to map with all domain integration cells on this proc
    const std::map<int,DomainIntCells>* elementalDomainIntCells() const {return &elementalDomainIntCells_;};
    //! get pointer to map with all boundary integration cells on this proc
    const std::map<int,BoundaryIntCells>* elementalBoundaryIntCells() const {return &elementalBoundaryIntCells_;};
    
    //! get pointer to map with all elements gids sorted by their interface labels
    const std::map<int,std::set<int> >* elementsByLabel() const {return &elementsByLabel_;};

    //! return list of domain integrationcells for a given global element id, returned list can be empty -> always terminates
    DomainIntCells GetDomainIntCells(
            const int                              gid,     ///< unique global element id
            const DRT::Element::DiscretizationType distype  ///< shape of thefinite element
            ) const;

    //! return list of boundary integrationcells for a given global element id, returned list can be empty -> always terminates
    BoundaryIntCells GetBoundaryIntCells(
            const int                              gid      ///< unique global element id
            ) const;

//    set<int> getSetOfIntersectedElementIds() const
//    {
//        set<int> eleidset;
//        map<int,DomainIntCells>::const_iterator tmp;
//        for (tmp = elementalDomainIntCells_.begin(); tmp != elementalDomainIntCells_.end(); ++tmp) {
//            const int eleid = tmp->first;
//            eleidset.insert(eleid);
//        }
//        return eleidset;
//    }

    /*! return
     * \return true, if element is intersected
     */
    bool ElementIntersected(
            const int element_gid   ///< unique global element id
            ) const;

    /*!
     * get pointer to boundary element
     * \return pointer to boundary element
     */
    const DRT::Element* GetBoundaryEle(
            const int element_gid   ///< unique global element id of a boundary element
            ) const
    {
        return cutterdis_->gElement(element_gid);
    }

private:
    //! intersected domain
    const RCP<DRT::Discretization> xfemdis_;

    //! discretization that intersects the xfemdis_ -> should be a boundary discretization (surface of structure)
    const RCP<DRT::Discretization> cutterdis_;
    
    //! ///< current positions of all cutter nodes
    std::map<int,BlitzVec3>        currentcutterpositions_;

    //! integration cells for domain (global id for element -> list of volume integration cells)
    map<int, DomainIntCells > elementalDomainIntCells_;

    //! integration cells for boundary
    map<int, BoundaryIntCells > elementalBoundaryIntCells_;
    
    //! elements grouped by the xfem condition label
    map<int,set<int> >        elementsByLabel_;

};

}

#endif  // #ifndef INTERFACE_H
#endif  // #ifdef CCADISCRET
