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
    
    DRT::Element* GetBoundaryEle(
            const int element_gid
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
    RCP<DRT::Discretization> xfemdis_;
	
	//! domain that intersects the xfemdis_
	RCP<DRT::Discretization> cutterdis_;
	
	//! integration cells for domain
	map<int, DomainIntCells > elementalDomainIntCells_;
	
	//! integration cells for boundary
	map<int, BoundaryIntCells > elementalBoundaryIntCells_;
	
	//! list of boundary elements that intersect volume elements on this processor
	//! first integer is the global surface ele id
	map<int, DRT::Element* > boundaryElements_;
	
	//! store all additional nodes for the cutter elements that are not already on this proc
	map< int, RCP<DRT::Node> > cutterNodeMap_;

};

}

#endif  // #ifndef INTERFACE_H
#endif  // #ifdef CCADISCRET
