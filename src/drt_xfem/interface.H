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
            const RCP<DRT::Discretization>        submerseddis  ///< discretization with structure elements
            );
    //! destructor
    ~InterfaceHandle();

    //! return string representation of the interface handle
    string toString() const;

    //! return pointer to the intersected discretization
    RCP<DRT::Discretization> xfemdis() const {return xfemdis_;};
    //! return pointer to the intersecting discretization
    RCP<DRT::Discretization> cutterdis() const {return cutterdis_;};

    //! return pointer to the sumbersed discretization (only used for enrichment calulation -> will go away
    RCP<DRT::Discretization> submerseddis() const {return submerseddis_;};

    //! get copy of map with all domain integration cells on this proc
    std::map<int,DomainIntCells> elementalDomainIntCells() const {return elementalDomainIntCells_;};
    //! get copy if map with all boundary integration cells on this proc
    std::map<int,BoundaryIntCells> elementalBoundaryIntCells() const {return elementalBoundaryIntCells_;};

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
    DRT::Element* GetBoundaryEle(
            const int element_gid   ///< unique global element id of a boundary element
            ) const
    {
        map<int, DRT::Element* >::const_iterator tmp = boundaryElements_.find(element_gid);
        if (tmp == boundaryElements_.end()){
            cout << element_gid << endl;
            dserror("GetBoundaryEle: no copy of boundary element with id = %d found on this proc!", element_gid);
            return NULL;
        }
        return tmp->second;
    }

private:
    //! intersected domain
    const RCP<DRT::Discretization> xfemdis_;

    //! domain that intersects the xfemdis_, should be a boundary discretization (surface of structure)
    RCP<DRT::Discretization> cutterdis_;

    //! domain (structure) that is submersed in the xfemdis_, used for enrichment evaluation -> will go away
    const RCP<DRT::Discretization> submerseddis_;

    //! integration cells for domain (global id for element -> list of volume integration cells)
    map<int, DomainIntCells > elementalDomainIntCells_;

    //! integration cells for boundary
    map<int, BoundaryIntCells > elementalBoundaryIntCells_;

    /*! list of boundary elements that intersect volume elements on this processor
     *  first integer is the global surface ele id to allow easy access
     */
    map<int, DRT::Element* > boundaryElements_;

    //! store all additional nodes for the cutter elements that are not already on this proc (int = gid of node)
    map<int, RCP<DRT::Node> > cutterNodeMap_;
    
    //! elements grouped by the xfem condition label
    map<int,set<DRT::Element*> >  elementsByLabel_;

};

}

#endif  // #ifndef INTERFACE_H
#endif  // #ifdef CCADISCRET
