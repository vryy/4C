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

#include "../drt_geometry/intersection.H"
#include "../drt_geometry/searchtree.H"
#include "geometry_service.H"

namespace XFEM
{

class DRT::Discretization;
  
/*!
 \brief It knows everything about the XFEM enrichments for one particular interface setup

 */
class InterfaceHandle
{
public:

    //! constructor
    explicit InterfaceHandle(
            const Teuchos::RCP<DRT::Discretization>        xfemdis,      ///< discretization that will be intersected
            const Teuchos::RCP<DRT::Discretization>        cutterdis,    ///< discretization with boundary elements
            const Epetra_Vector&                           idispcol      ///< interface displacements for current position of interface
            );
    //! destructor
    ~InterfaceHandle();

    //! return string representation of the interface handle
    std::string toString() const;
    
    //! return string representation of the interface handle
    void toGmsh(const int step) const;

    //! return pointer to the intersected discretization
    Teuchos::RCP<DRT::Discretization> xfemdis() const {return xfemdis_;};
    //! return pointer to the intersecting discretization
    Teuchos::RCP<DRT::Discretization> cutterdis() const {return cutterdis_;};
    
    //! returns pointer to map with (current) nodal positions per node of the cutter discretization 
    const std::map<int,BlitzVec3>* currentcutterpositions() const {return &currentcutterpositions_;};

    //! get pointer to map with all domain integration cells on this proc
    const std::map<int,GEO::DomainIntCells>* elementalDomainIntCells() const {return &elementalDomainIntCells_;};
    //! get pointer to map with all boundary integration cells on this proc
    const std::map<int,GEO::BoundaryIntCells>* elementalBoundaryIntCells() const {return &elementalBoundaryIntCells_;};
    
    //! get pointer to elementsByLabel_
    const std::map<int,std::set<int> >* elementsByLabel() const {return &elementsByLabel_;};

    //! get pointer to oct tree
    Teuchos::RCP<GEO::SearchTree> getSearchTree() const {return octTree_;};
 
    
    //! return list of domain integrationcells for a given global element id, returned list can be empty -> always terminates
    GEO::DomainIntCells GetDomainIntCells(
            const int                              gid,     ///< unique global element id
            const DRT::Element::DiscretizationType distype  ///< shape of thefinite element
            ) const;

    //! return list of boundary integrationcells for a given global element id, returned list can be empty -> always terminates
    GEO::BoundaryIntCells GetBoundaryIntCells(
            const int                              gid      ///< unique global element id
            ) const;

    /*! tell me, whether the element is intersected or not
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

    /// Print Information about the contained search tree
    void PrintTreeInformation(const int step) const
    {
        printf("tree metrics is not yet implemented");
        octTree_->printTree("TREE", step);
        return;
    }
    
private:
  
    //! don't want default constructor
    explicit InterfaceHandle();
    
    //! don't want copy constructor
    explicit InterfaceHandle(const InterfaceHandle&);
    
    //! put a position into the map for each nodal id
    void FillCurrentCutterPositionMap(
        const Teuchos::RCP<DRT::Discretization>  cutterdis,              ///< dis containing nodes
        const Epetra_Vector&                     idispcol,               ///< vector with displacements
        std::map<int,BlitzVec3>&                 currentcutterpositions  ///< map to be filled with positions
        ) const;
    
    //! for difficult cases during intersection, do sanity checks and create list of redundant integration cells 
    std::set<int> FindDoubleCountedIntersectedElements() const;
  
    //! intersected domain
    const Teuchos::RCP<DRT::Discretization> xfemdis_;

    //! discretization that intersects the xfemdis_ -> should be a boundary discretization (surface of structure)
    const Teuchos::RCP<DRT::Discretization> cutterdis_;
    
    //! current positions of all cutter nodes
    std::map<int,BlitzVec3>                 currentcutterpositions_;

    //! integration cells for domain (global id for element -> list of volume integration cells)
    std::map<int,GEO::DomainIntCells >          elementalDomainIntCells_;

    //! integration cells for boundary
    std::map<int,GEO::BoundaryIntCells >        elementalBoundaryIntCells_;
    
    //! for each xfem condition label there is a set of boundary elements
    std::map<int,std::set<int> >            elementsByLabel_;
    
    //! SearchTree octtree
    Teuchos::RCP<GEO::SearchTree>              octTree_;

};

/*!
\brief checks if a position in current coordinates lies within a certain condition region           
*/  
int PositionWithinCondition(
        const BlitzVec3&                  x,              ///< node in physical coordinates
        const XFEM::InterfaceHandle&      ih              ///< the ubiquitous interface handle with all geometry information
    );


/*!
\brief checks if a position in current coordinates lies within a certain condition region using an octree
*/
int PositionWithinConditionSearchTree(
    const BlitzVec3&                  x_in, ///< node in physical coordinates
    const XFEM::InterfaceHandle&      ih    ///< the ubiquitious interface handle with all geometry information
    );


/*!
\brief checks if a position in current coordinates lies within a certain condition region           
*/  
bool PositionWithinAnyInfluencingCondition(
        const BlitzVec3&                  x,              ///< node in physical coordinates
        const XFEM::InterfaceHandle&      ih,             ///< the ubiquitous interface handle with all geometry information
        const std::set<int>&              xlabelset       ///< influencing enrichment labels
    );


/*!
\brief checks if a position in current coordinates lies within a certain condition region: example brute force implementation
\warning does worrk only for convex bodies in serial and it is the worst case algorithm O(n*m) - use for debuging only
*/
int PositionWithinConditionBruteForce(
        const BlitzVec3&                  x,              ///< node in physical coordinates
        const XFEM::InterfaceHandle&      ih              ///< the ubiquitous interface handle with all geometry information
    );

}

#endif  // #ifndef INTERFACE_H
#endif  // #ifdef CCADISCRET
