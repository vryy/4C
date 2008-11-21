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

#include <set>
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


namespace GEO
{

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
    DRT::Element::DiscretizationType distype_;

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
            const DRT::Element::DiscretizationType    distype,            ///< distype of the integration cell
            const std::vector<std::vector<double> >&  domainCoordinates   ///< coordinates in parent element coordinates
            );

    //! Standard Constructor
    explicit DomainIntCell(
            const DRT::Element::DiscretizationType  distype,            ///< distype of the integration cell
            const LINALG::SerialDenseMatrix&        domainCoordinates   ///< coordinates in parent element coordinates
            );
    
    /*!
     \brief Standard Constructor for dummy integration cell

     This is the default integration cell, if no intersections are within the XFEM element
     
     the user just gives the distype of the XFEM element and gets an integration cell with the same distype as the parent
     */
    explicit DomainIntCell(
            const DRT::Element::DiscretizationType distype  ///< distype of the parent elementintegration cell
            );
    
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
    const LINALG::SerialDenseMatrix* NodalPosXiDomain() const {return &nodalpos_xi_domain_;};
    
    //! fills an array with the coordinates of the integration cell in physical coordinates
    void NodalPosXYZ(
        const DRT::Element&         ele,     ///< parent element
        LINALG::SerialDenseMatrix&  xyz_cell ///< fill with nodal positions
        ) const;
    
    //! return string representation of the DomainIntCell
    std::string toString() const;
    
    //! get geometric center of the cell in physical coordinates
    LINALG::Matrix<3,1> GetPhysicalCenterPosition(const DRT::Element& ele) const;

private:
    
    //! coordinates of the nodes in the parent domain element coordinates xi_domain
    LINALG::SerialDenseMatrix nodalpos_xi_domain_;
    
    //! hide default constructor
    explicit DomainIntCell();
    
    //! this function returns the coordinate array for all nodes of an element in parameter coordinates \f$ \xi \f$
    LINALG::SerialDenseMatrix GetDefaultCoordinates(
            const DRT::Element::DiscretizationType distype  ///< shape of element
            ) const;
};

//! shortcut for a vector of DomainIntCells
typedef std::vector<GEO::DomainIntCell > DomainIntCells;
//! shortcut for a set of DomainIntCells
typedef std::set<GEO::DomainIntCell > DomainIntCellSet;


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
            const std::vector<std::vector<double> >&          domainCoordinates,  ///< coordinates of the integrationcell in xsi_domain
            const std::vector<std::vector<double> >&          boundaryCoordinates ///< coordinates of the integrationcell in xsi_boundary
            );

    //! Standard Constructor
    explicit BoundaryIntCell(
            const DRT::Element::DiscretizationType  distype,             ///< shape of the integration cell
            const int                               surface_ele_gid,     ///< global id of the boundary element (cutter)
            const LINALG::SerialDenseMatrix&        domainCoordinates,   ///< coordinates of the integrationcell in xsi_domain
            const LINALG::SerialDenseMatrix&        boundaryCoordinates  ///< coordinates of the integrationcell in xsi_boundary
            );
    
    //! Copy Constructor
    BoundaryIntCell(
            const BoundaryIntCell& old);
    
    //! fills an array with the coordinates of the integration cell in physical coordinates
    void NodalPosXYZ(
        const DRT::Element&         ele,     ///< parent element
        LINALG::SerialDenseMatrix&  xyz_cell ///< fill with nodal positions
        ) const;
    
    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    const LINALG::SerialDenseMatrix* NodalPosXiDomain() const {return &nodalpos_xi_domain_;};
    
    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    const LINALG::SerialDenseMatrix* NodalPosXiBoundary() const {return &nodalpos_xi_boundary_;};
    
    //! return "parent" cutter element id (global id)
    int GetSurfaceEleGid() const {return surface_ele_gid_;}
    
    //! return string representation of the BoundaryIntCell
    std::string toString() const;
    
    //! get geometric center of the cell in physical coordinates
    LINALG::Matrix<3,1> GetPhysicalCenterPosition(const DRT::Element& ele) const;

private:
    
    //! hide default constructor
    explicit BoundaryIntCell();
    
    //! the boundaryIntCell should know, to which cutterElement it belongs!?
    int surface_ele_gid_;
    
    //! coordinates of the nodes
    LINALG::SerialDenseMatrix nodalpos_xi_domain_;
    
    //! boundary coordinates of the nodes
    LINALG::SerialDenseMatrix nodalpos_xi_boundary_;
    
};

//! shortcut for a vector of BoundaryIntCells
typedef std::vector<GEO::BoundaryIntCell > BoundaryIntCells;
//! shortcut for a set of BoundaryIntCells
typedef std::set<GEO::BoundaryIntCell > BoundaryIntCellSet;

}

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef CCADISCRET
