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
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "xfem.H"


namespace XFEM
{

/*!
 * \brief create array with physical coordinates based an local coordinates of a parent element
 */
template<class Cell>
static void ComputePhysicalCoordinates(
        const DRT::Element&  ele,  ///< parent element
        const Cell&          cell, ///< integration cell whose coordinates we'd like to transform
        BlitzMat&            physicalCoordinates
        )
{
    const BlitzMat eleCoord(DRT::UTILS::InitialPositionArrayBlitz(&ele));
    //DRT::UTILS::fillInitialPositionArray(&ele, eleCoord);
    const BlitzMat* nodalPosXiDomain = cell.NodalPosXiDomainBlitz();
    
    const int nen_cell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // return value
    //BlitzMat physicalCoordinates(3, nen_cell);
    physicalCoordinates = 0.0;
    // for each cell node, compute physical position
    const int nen_ele = ele.NumNode();
    BlitzVec funct(nen_ele);
    for (int inen = 0; inen < nen_cell; ++inen)
    {
        // shape functions
        DRT::UTILS::shape_function_3D(funct,
                (*nodalPosXiDomain)(0, inen),
                (*nodalPosXiDomain)(1, inen),
                (*nodalPosXiDomain)(2, inen),
                ele.Shape());

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < nen_ele; ++j)
            {
              physicalCoordinates(i,inen) += eleCoord(i, j) * funct(j);
            }
        }
    };
    return;
}



/*!
 * \brief An integration cell is used for specialized element integration routines 
 */
class IntCell
{

public:

    //! Standard Constructor
    explicit IntCell(
        const DRT::Element::DiscretizationType  distype  ///< distype of the integration cell
        );

    //! Copy Constructor
    explicit IntCell(const IntCell& old);
    
    //! Destructor
    virtual ~IntCell();

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
    DRT::Element::DiscretizationType Shape() const {return distype_;};
    
    /*!
     \brief returns number of nodes

     \return number of nodes
     */
    int NumNode() const {return DRT::UTILS::getNumberOfElementNodes(Shape());};
    
    //! returns a string representation of the IntCell
    std::string toString() const;

private:

    //! shape
    const DRT::Element::DiscretizationType distype_;

    //! hidden default constructor
    explicit IntCell();
};


/*!
 * \brief An domain integration cell is used for integrating the XFEM elements domain 
 */
class DomainIntCell : public IntCell
{

public:

    //! Standard Constructor
    explicit DomainIntCell(
            const DRT::Element::DiscretizationType  distype,            ///< distype of the integration cell
            const vector<vector<double> >&          domainCoordinates   ///< coordinates in parent element coordinates
            );

    //! Standard Constructor
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
        DomainIntCell cell(old.Shape(), old.nodalpos_xi_domain_blitz_);
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
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    const BlitzMat* NodalPosXiDomainBlitz() const {return &nodalpos_xi_domain_blitz_;};
    
    /*!
     \brief returns the coordinates of the integration cell in physical coordinates

     \return cells nodal coordinates  
     */
    void NodalPosXYZ(const DRT::Element& ele, BlitzMat& xyz_cell) const
    {
        ComputePhysicalCoordinates(ele, (*this), xyz_cell);
        return;
    };
    
    //! return string representation of the DomainIntCell
    std::string toString() const;
    
    //! get geometric center of the cell in physical coordinates
    BlitzVec3 GetPhysicalCenterPosition(const DRT::Element& ele) const;

private:
    
    //! coordinates of the nodes in the parent domain element coordinates xi_domain
    const BlitzMat nodalpos_xi_domain_blitz_;
    
    //! hide default constructor
    explicit DomainIntCell();
    
    //! this function returns the coordinate array for all nodes of an element in parameter coordinates \f$ \xi \f$
    BlitzMat GetDefaultCoordinates(
            const DRT::Element::DiscretizationType distype  ///< shape of element
            ) const;
};

//! shortcut for a vector of DomainIntCells
typedef std::vector<XFEM::DomainIntCell > DomainIntCells;
//! shortcut for a set of DomainIntCells
typedef std::set<XFEM::DomainIntCell > DomainIntCellSet;


/*!
 * \brief An boundary integration cell is used for integrating along a discontinuity caused by XFEM enrichments 
 */
class BoundaryIntCell : public IntCell
{

public:

    //! Standard Constructor
    explicit BoundaryIntCell(
            const DRT::Element::DiscretizationType  distype,            ///< shape of the integration cell
            const int                               surface_ele_gid,    ///< global id of the boundary element (cutter)
            const vector<vector<double> >&          domainCoordinates,  ///< coordinates of the integrationcell in xsi_domain
            const vector<vector<double> >&          boundaryCoordinates ///< coordinates of the integrationcell in xsi_boundary
            );

    //! Standard Constructor
    explicit BoundaryIntCell(
            const DRT::Element::DiscretizationType  distype,             ///< shape of the integration cell
            const int                               surface_ele_gid,     ///< global id of the boundary element (cutter)
            const BlitzMat&                         domainCoordinates,   ///< coordinates of the integrationcell in xsi_domain
            const BlitzMat&                         boundaryCoordinates  ///< coordinates of the integrationcell in xsi_boundary
            );
    
    //! Copy Constructor
    BoundaryIntCell(
            const BoundaryIntCell& old);

    /*!  = operator
     * create true copy -> allows to have constant instance members
     */
    BoundaryIntCell& operator = (const BoundaryIntCell& old)
    {
        BoundaryIntCell cell(
                old.Shape(),
                old.surface_ele_gid_,
                old.nodalpos_xi_domain_blitz_,
                old.nodalpos_xi_boundary_blitz_);
        BoundaryIntCell& ref = cell;
        return ref;
    };
    
    /*!
     \brief returns the coordinates of the integration cell in physical coordinates

     \return cells nodal coordinates  
     */
    void NodalPosXYZ(
            const DRT::Element& ele,  ///< parent element, where the integration cell is located in
            BlitzMat&           xyz_cell
            ) const
    {
        ComputePhysicalCoordinates(ele, (*this), xyz_cell);
        return;
    };
    
    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    const BlitzMat* NodalPosXiDomainBlitz() const {return &nodalpos_xi_domain_blitz_;};
    
    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    const BlitzMat* NodalPosXiBoundaryBlitz() const {return &nodalpos_xi_boundary_blitz_;};
    
    //! return "parent" cutter element id (global id)
    int GetSurfaceEleGid() const {return surface_ele_gid_;}
    
    //! return string representation of the BoundaryIntCell
    std::string toString() const;
    
    //! get geometric center of the cell in physical coordinates
    BlitzVec3 GetPhysicalCenterPosition(const DRT::Element& ele) const;

private:
    
    //! hide default constructor
    explicit BoundaryIntCell();
    
    //! the boundaryIntCell should know, to which cutterElement it belongs!?
    const int surface_ele_gid_;
    
    //! coordinates of the nodes
    const BlitzMat nodalpos_xi_domain_blitz_;
    
    //! boundary coordinates of the nodes
    const BlitzMat nodalpos_xi_boundary_blitz_;
    
};

//! shortcut for a vector of BoundaryIntCells
typedef std::vector<XFEM::BoundaryIntCell > BoundaryIntCells;
//! shortcut for a set of BoundaryIntCells
typedef std::set<XFEM::BoundaryIntCell > BoundaryIntCellSet;

}

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef CCADISCRET
