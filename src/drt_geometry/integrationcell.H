/*!----------------------------------------------------------------------
\file integrationcell.H

\brief integration cell classes for domain and boundary integration

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef INTEGRATIONCELL_H
#define INTEGRATIONCELL_H

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_geometry/element_coordtrafo.H"


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
        const DRT::Element::DiscretizationType&  distype  ///< distype of the integration cell
        );

    //! Copy Constructor
    explicit IntCell(const IntCell& old);
    
    //! virtual destructor
    virtual ~IntCell();

    //! brief returns the shape of the integration cell
    DRT::Element::DiscretizationType Shape() const {return distype_;};
    
    //! brief returns number of nodes
    int NumNode() const {return DRT::UTILS::getNumberOfElementNodes(Shape());};
    
    //! returns a string representation of the IntCell
    virtual std::string toString() const;


  private:

    //! hidden default constructor
    explicit IntCell();
    
    //! shape
    DRT::Element::DiscretizationType distype_;   
};




/*!
 * \brief An domain integration cell is used for integrating the XFEM elements domain 
 */
class DomainIntCell : public IntCell
{
  public:

    //! Standard Constructor
    explicit DomainIntCell(
        const DRT::Element::DiscretizationType&     distype,                    ///< distype of the integration cell
        const LINALG::SerialDenseMatrix&            xfemEleDomainCoordinates,   ///< coordinates in xfem element coordinates
        const LINALG::SerialDenseMatrix&            physDomainCoordinates       ///< coordinates in physical coordinates
        );
    
    /*!
     \brief Standard Constructor for dummy integration cell
            This is the default integration cell, if no intersections 
            are within the XFEM element
            the user just gives the distype of the XFEM element and 
            gets an integration cell with the same distype as the parent
     */
    explicit DomainIntCell(
        const DRT::Element::DiscretizationType&     distype,    ///< distype of the parent elementintegration cell
        const LINALG::SerialDenseMatrix&            xyz_ele     ///< coordinates in physical coordinates
        );
    
    //! brief Copy Constructor
    DomainIntCell(const DomainIntCell& old);
         
    //! virtual destructor
    virtual ~DomainIntCell();
    
     //! returns the coordinates of the integration cell in domain element coordinates
    const LINALG::SerialDenseMatrix&  CellNodalPosXiDomain() const {return nodalpos_xi_domain_;};
    
    //! returns an array with the coordinates of the integration cell in physical coordinates
    const LINALG::SerialDenseMatrix&  CellNodalPosXYZ() const {return nodalpos_xyz_domain_;}
    
    //! returns an array with the coordinates of the integration cell in physical coordinates
    const LINALG::Matrix<3,1>&  GetPhysicalCenterPosition() const{ return phys_center_;}
    
    double SizeXiDomain(const DRT::Element& ele) const;
       
    //! return string representation of the DomainIntCell
    virtual std::string toString() const;
    
    //! return label of domain intcell
    int getLabel() const {return label_; }
        
    //! set label of domain intcell
    void setLabel(const int   label);
    

  private:
     
    //! hide default constructor
    explicit DomainIntCell();
    
    //! get geometric center of the cell in physical coordinates
    LINALG::Matrix<3,1> ComputePhysicalCenterPosition(
        const DRT::Element::DiscretizationType&   distype,
        const LINALG::SerialDenseMatrix&          xyze
        ) const;
      
    //! coordinates of the nodes in the parent element coordinates system ( xi_domain )
    LINALG::SerialDenseMatrix                       nodalpos_xi_domain_;
    
    //! coordinates of the nodes in the physical coordinates system ( xyz_domain )
    LINALG::SerialDenseMatrix                       nodalpos_xyz_domain_;
    
    //! get geometric center of the cell in physical coordinates
    LINALG::Matrix<3,1>                             phys_center_;
    
    // xfem label if fluid label = 0; otherwise label = solid id
    int label_;
};

//! shortcut for a vector of DomainIntCells
typedef std::vector<GEO::DomainIntCell > DomainIntCells;




/*!
 * \brief An boundary integration cell is used for 
 *        integrating at a discontinuity caused by XFEM enrichments 
 */
class BoundaryIntCell : public IntCell
{
  public:
  
    //! Standard Constructor
    explicit BoundaryIntCell(
      const DRT::Element::DiscretizationType&     distype,                    ///< shape of the integration cell
      const int                                   surface_ele_gid,            ///< global id of the boundary element (cutter)
      const LINALG::SerialDenseMatrix&            xfemEleDomainCoordinates,   ///< coordinates in element parameter space xsi
      const LINALG::SerialDenseMatrix&            eleBoundaryCoordinates,     ///< coordinates in boundary parameter space eta 
      const LINALG::SerialDenseMatrix&            physDomainCoordinates       ///< coordinates of the integrationcell in physical domain
      );
    
    //! Copy Constructor
    BoundaryIntCell( const BoundaryIntCell& old);
      
    //! virtual destructor
    virtual ~BoundaryIntCell();
    
    //! returns the coordinates of the integration cell in parent element coordinates xsi
    const LINALG::SerialDenseMatrix&  CellNodalPosXiDomain() const {return nodalpos_xi_domain_;};
    
    //! returns the coordinates of the integration cell in boundary parent space eta 
    const LINALG::SerialDenseMatrix&  CellNodalPosXiBoundary() const {return nodalpos_xi_boundary_;};
    
    //! returns an array with the coordinates of the integration cell in physical coordinates
    const LINALG::SerialDenseMatrix&  CellNodalPosXYZ() const{ return nodalpos_xyz_domain_;}
    
    //! return "parent" cutter element id (global id)
    int GetSurfaceEleGid() const {return surface_ele_gid_;}
        
    //! return string representation of the BoundaryIntCell
    virtual std::string toString() const;


  private:
    
    //! hide default constructor
    explicit BoundaryIntCell();
    
    //! the boundaryIntCell should know, to which cutterElement it belongs!?
    int                                   surface_ele_gid_;
    
    //! coordinates of the nodes of the integration cell in parent element coordinates xsi
    LINALG::SerialDenseMatrix             nodalpos_xi_domain_;
    
    //! boundary coordinates of the nodes of the integration cell in boundary element coordinates eta
    LINALG::SerialDenseMatrix             nodalpos_xi_boundary_;
     
    //! coordinates of the nodes of the integration cell in physical coordinates
    LINALG::SerialDenseMatrix             nodalpos_xyz_domain_;
};

//! shortcut for a vector of BoundaryIntCells
typedef std::vector<GEO::BoundaryIntCell > BoundaryIntCells;


}

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef CCADISCRET
