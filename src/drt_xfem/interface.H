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
#include "../io/gmsh.H"

namespace XFEM
{

/*!
 \class InterfaceHandle

 \brief He knows everything about the XFEM enrichments for one particular interface setup

 */
class InterfaceHandle
{
public:

	explicit InterfaceHandle(
			const RCP<DRT::Discretization>        xfemdis, 
			const RCP<DRT::Discretization>        cutterdis);
    ~InterfaceHandle();

    string toString() const;
    
    RCP<DRT::Discretization> xfemdis() const {return xfemdis_;};
    RCP<DRT::Discretization> cutterdis() const {return cutterdis_;};
    
    // ok, we copy the whole thing, this can be more elegant later on or not used at all
    map<int,DomainIntCells> elementalDomainIntCells() const {return elementalDomainIntCells_;};
    map<int,BoundaryIntCells> elementalBoundaryIntCells() const {return elementalBoundaryIntCells_;};

    // return list of integrationcells, however the list can be empty -> always terminates
    DomainIntCells GetDomainIntCells(
            const int gid,
            const DRT::Element::DiscretizationType distype
            ) const;
    
    BoundaryIntCells GetBoundaryIntCells(
            const int gid
            ) const;
    
//    inline set<int> getSetOfIntersectedElementIds() const
//    {
//        set<int> eleidset;
//        map<int,DomainIntCells>::const_iterator tmp;
//        for (tmp = elementalDomainIntCells_.begin(); tmp != elementalDomainIntCells_.end(); ++tmp) {
//            const int eleid = tmp->first;
//            eleidset.insert(eleid);
//        }
//        return eleidset;
//    }

    bool PositionWithinSpecificClosedRegion(
            const blitz::Array<double,1>& actpos,
            const int xfemcondition_label
            ) const;
    
    bool ElementIntersected(
            const int element_gid
            ) const;
    
private:
	// intersected domain
	RefCountPtr<DRT::Discretization> xfemdis_;
	
	// domain that intersects the xfemdis_
	RefCountPtr<DRT::Discretization> cutterdis_;
	
	// integration cells for domain
	map<int, DomainIntCells > elementalDomainIntCells_;
	
	// integration cells for boundary
	map<int, BoundaryIntCells > elementalBoundaryIntCells_;
};

}

#endif  // #ifndef INTERFACE_H
#endif  // #ifdef CCADISCRET
