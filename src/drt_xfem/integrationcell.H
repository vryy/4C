/*!
 \file integrationcell.H

 \brief integration cell classes for domain and boundary integration

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */

#ifdef CCADISCRET
#ifndef INTEGRATIONCELL_H
#define INTEGRATIONCELL_H

#include <vector>
#include <set>
#include <string>
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_integration.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "xfem.H"


namespace XFEM
{

//! little helper function
inline BlitzMat ConvertPosArrayToBlitz(
        const vector<vector<double> >&         pos_array,
        const DRT::Element::DiscretizationType distype,
        const int                              dim
        )
{
    const int numnode = DRT::UTILS::getNumberOfElementNodes(distype);
    BlitzMat pos_array_blitz(dim,numnode,blitz::ColumnMajorArray<2>());
    for (int inode=0; inode<numnode; ++inode)
    {
        for (int isd=0; isd<dim; ++isd)
        {
            pos_array_blitz(isd,inode) = pos_array[inode][isd];
        }
    }    
    return pos_array_blitz;
}


template<class Cell>
BlitzMat ComputePhysicalCoordinates(
        const DRT::Element&  ele,
        const Cell& cell
        )
{
    const BlitzMat eleCoord(DRT::UTILS::PositionArrayBlitz(&ele));
    const BlitzMat nodalPosXiDomain(cell.NodalPosXiDomainBlitz());
    
    const int nen_cell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // return value
    BlitzMat physicalCoordinates(3, nen_cell);
    physicalCoordinates = 0.0;
    // for each cell node, compute physical position
    for (int inen = 0; inen < nen_cell; ++inen)
    {
        // shape functions
        const BlitzVec funct(DRT::UTILS::shape_function_3D(
                nodalPosXiDomain(0, inen),
                nodalPosXiDomain(1, inen),
                nodalPosXiDomain(2, inen),
                ele.Shape()));

        blitz::firstIndex i;    // Placeholder for the first index
        blitz::secondIndex j;   // Placeholder for the second index
        const blitz::Range _  = blitz::Range::all();
        physicalCoordinates(_,inen) = blitz::sum(eleCoord(i, j) * funct(j),j);
    };
    return physicalCoordinates;
}



/*
 * \brief An integration cell is used for specialized element integration routines 
 */
class IntCell
{

public:

    /*!
     \brief Standard Constructor
     */
    explicit IntCell(
        const DRT::Element::DiscretizationType  distype  ///< distype of the integration cell
        );

    /*!
     \brief Copy Constructor

     Makes a deep copy of a Element
     */
    explicit IntCell(const IntCell& old);
    
    /*!
    \brief Destructor

    */
    ~IntCell();

    /*!  = operator
     * create true copy -> allows to have constant instance members
     */
    IntCell& operator = (const IntCell& old)
    {
        IntCell intcell(old.distype_);
        IntCell& ref = intcell;
        return ref;
    };
    
    /*!
     \brief returns the shape of the integration cell

     \return distype_
     */
    inline const DRT::Element::DiscretizationType Shape() const {return distype_;};
    
    inline const std::string toString() const;

protected:

    //! shape
    const DRT::Element::DiscretizationType distype_;
    
private:
    
    explicit IntCell();
};


/*
 * \brief An domain integration cell is used for integrating the XFEM elements domain 
 */
class DomainIntCell : IntCell
{

public:

    /*!
     \brief Standard Constructor
     */
    explicit DomainIntCell(
            const DRT::Element::DiscretizationType  distype,            ///< distype of the integration cell
            const vector<vector<double> >&          domainCoordinates   ///< coordinates in parent element coordinates
            );

    /*!
     \brief Standard Constructor
     */
    explicit DomainIntCell(
            const DRT::Element::DiscretizationType  distype,            ///< distype of the integration cell
            const BlitzMat&                         domainCoordinates   ///< coordinates in parent element coordinates
            );
    
    /*!
     \brief Standard Constructor for dummy integration cell

     This is the default integration cell, if no intersections are within the XFEM element
     
     the user just gives the distype of the XFEM element and gets an integration cell with the same distype as the parent
     */
    explicit DomainIntCell(
            const DRT::Element::DiscretizationType distype  ///< distype of the parent elementintegration cell
            );
    
    /*!  = operator
     * create true copy -> allows to have constant instance members
     */
    DomainIntCell& operator = (const DomainIntCell& old)
    {
        DomainIntCell cell(old.distype_, old.nodalpos_xi_domain_blitz_);
        DomainIntCell& ref = cell;
        return ref;
    };
    
    /*!
     \brief Copy Constructor

     Makes a deep copy of a Element
     */
    DomainIntCell(
            const DomainIntCell& old);

    /*!
     \brief returns the shape of the integration cell

     \return distype_
     */
    inline const DRT::Element::DiscretizationType Shape() const {return distype_;};
    
    /*!
     \brief returns number of nodes

     \return number of nodes
     */
    inline const int NumNode() const {return DRT::UTILS::getNumberOfElementNodes(distype_);};
    
    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    inline const BlitzMat NodalPosXiDomainBlitz() const {return nodalpos_xi_domain_blitz_;};
    
    /*!
     \brief returns the coordinates of the integration cell in physical coordinates

     \return cells nodal coordinates  
     */
    inline const BlitzMat NodalPosXYZ(const DRT::Element& ele) const
    {
        return ComputePhysicalCoordinates<DomainIntCell>(ele, (*this));
    };
    
    const std::string toString() const;
    
    const BlitzVec GetPhysicalCenterPosition(const DRT::Element& ele) const;

private:
    
    //! coordinates of the nodes in the parent domain element coordinates xi_domain
    const BlitzMat nodalpos_xi_domain_blitz_;
    
    explicit DomainIntCell();
    
    std::vector<std::vector<double> > GetDefaultCoordinates(
            const DRT::Element::DiscretizationType distype) const;
};

typedef std::vector<XFEM::DomainIntCell > DomainIntCells;
typedef std::set<XFEM::DomainIntCell > DomainIntCellSet;


/*
 * \brief An domain integration cell is used for integrating the XFEM elements domain 
 */
class BoundaryIntCell : IntCell
{

public:

    /*!
     \brief Standard Constructor

     \param coordinates      : coordinates of the integrationcell
     \param distype      : coordinates of the integrationcell
     */
    explicit BoundaryIntCell(
            const DRT::Element::DiscretizationType    distype,
            const int                                 surface_ele_gid,
            const vector<vector<double> >&            domainCoordinates,
            const vector<vector<double> >&            boundaryCoordinates
            );

    /*!
     \brief Standard Constructor

     \param coordinates      : coordinates of the integrationcell
     \param distype      : coordinates of the integrationcell
     */
    explicit BoundaryIntCell(
            const DRT::Element::DiscretizationType    distype,
            const int                                 surface_ele_gid,
            const BlitzMat&                           domainCoordinates,
            const BlitzMat&                           boundaryCoordinates
            );
    
    /*!
     \brief Copy Constructor

     Makes a deep copy of a Element
     */
    BoundaryIntCell(
            const BoundaryIntCell& old);

    /*!  = operator
     * create true copy -> allows to have constant instance members
     */
    BoundaryIntCell& operator = (const BoundaryIntCell& old)
    {
        BoundaryIntCell cell(
                old.distype_,
                old.surface_ele_gid_,
                old.nodalpos_xi_domain_blitz_,
                old.nodalpos_xi_boundary_blitz_);
        BoundaryIntCell& ref = cell;
        return ref;
    };
    
    /*!
     \brief returns the shape of the integration cell

     \return distype_
     */
    const DRT::Element::DiscretizationType Shape() const {return distype_;};
    
    /*!
     \brief returns number of nodes

     \return number of nodes
     */
    inline const int NumNode() const {return DRT::UTILS::getNumberOfElementNodes(distype_);};
    
    /*!
     \brief returns the coordinates of the integration cell in physical coordinates

     \return cells nodal coordinates  
     */
    inline const BlitzMat NodalPosXYZ(const DRT::Element& ele) const
    {
        return ComputePhysicalCoordinates<BoundaryIntCell>(ele, (*this));
    };
    
    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    inline const BlitzMat NodalPosXiDomainBlitz() const {return nodalpos_xi_domain_blitz_;};
    
    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    inline const BlitzMat NodalPosXiBoundaryBlitz() const {return nodalpos_xi_boundary_blitz_;};
    
    inline const int GetSurfaceEleGid() const {return surface_ele_gid_;}
    
    const std::string toString() const;

private:
    
    explicit BoundaryIntCell();
    
    //! the boundaryIntCell should know, to which cutterElement it belongs!?
    const int surface_ele_gid_;
    
    //! coordinates of the nodes
    const BlitzMat nodalpos_xi_domain_blitz_;
    
    //! boundary coordinates of the nodes
    const BlitzMat nodalpos_xi_boundary_blitz_;
    
};

typedef std::vector<XFEM::BoundaryIntCell > BoundaryIntCells;
typedef std::set<XFEM::BoundaryIntCell > BoundaryIntCellSet;

}

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef CCADISCRET
