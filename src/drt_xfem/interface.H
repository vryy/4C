/*!
 \file integrationcell.H

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
    const map<int,DomainIntCells> elementalDomainIntCells() const {return elementalDomainIntCells_;};
    const map<int,BoundaryIntCells> elementalBoundaryIntCells() const {return elementalBoundaryIntCells_;};

    // return set of integrationcells, however the set can be empty -> always terminates
    inline const DomainIntCells domainIntCells(const int gid) const
    {
    	map<int,DomainIntCells>::const_iterator tmp = elementalDomainIntCells_.find(gid);
    	if (tmp == elementalDomainIntCells_.end())
    	{	
    		const DomainIntCells x;
    		return x;
    	}
        return tmp->second;
    };
    
    inline const BoundaryIntCells boundaryIntCells(const int gid) const
    {
    	map<int,BoundaryIntCells>::const_iterator tmp = elementalBoundaryIntCells_.find(gid);
    	if (tmp == elementalBoundaryIntCells_.end())
    	{	
    		const BoundaryIntCells x;
    		return x;
    	}
        return tmp->second;
    };

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

};

#endif  // #ifndef INTERFACE_H
#endif  // #ifdef CCADISCRET
