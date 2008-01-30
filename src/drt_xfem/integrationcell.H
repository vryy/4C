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
#include "../drt_lib/drt_utils_integration.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "xfem.H"


namespace XFEM
{

/*
 * \brief An integration cell is used for specialized element integration routines 
 */
class IntCell
{

public:

    /*!
     \brief Standard Constructor

     \param distype      : distype of the integrationcell
     */
    IntCell(
        const DRT::Element::DiscretizationType distype);

    /*!
     \brief Copy Constructor

     Makes a deep copy of a Element
     */
    IntCell(const IntCell& old);
    
    /*!
    \brief Destructor

    */
    virtual ~IntCell();

    /*!
     \brief returns the shape of the integration cell

     \return distype_
     */
    inline DRT::Element::DiscretizationType Shape() const {return distype_;};
    
    /*!
     \brief returns the coordinates of the integration cell

     \return cells nodal coordinates  
     */
    virtual std::vector< std::vector<double> > GetDomainCoord() const;
    
    /*!
     \brief returns the coordinates of the integration cell in physical coordinates

     \return cells nodal coordinates  
     */
    virtual std::vector< std::vector<double> > GetPhysicalCoord(DRT::Element& ele) const;
    
    virtual std::string Print() const;

protected:

    //! shape
    DRT::Element::DiscretizationType distype_;

    vector<vector<double> > ComputePhysicalCoordinates(
            DRT::Element&  ele) const;
};


/*
 * \brief An domain integration cell is used for integrating the XFEM elements domain 
 */
class DomainIntCell : IntCell
{

public:

    /*!
     \brief Standard Constructor

     \param distype          : distype of the integration cell
     \param coordinates      : coordinates of the integrationcell in parent element coordinates
     */
    DomainIntCell(
            const DRT::Element::DiscretizationType distype,
            const std::vector< std::vector<double> >& domainCoordinates);

    /*!
     \brief Standard Constructor for dummy integration cell

     This is the default integration cell, if no intersections are within the XFEM element
     
     the user just gives the distype of the XFEM element and gets an integration cell with the same distype as the parent

     \param distype          : distype of the parent elementintegration cell
     */
    DomainIntCell(
            const DRT::Element::DiscretizationType distype);
    
    
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
    inline DRT::Element::DiscretizationType Shape() const {return distype_;};
    
    /*!
     \brief returns number of nodes

     \return number of nodes
     */
    inline int NumNode() const {return DRT::UTILS::getNumberOfElementNodes(distype_);};
    
    /*!
     \brief returns the coordinates of the integration cell

     \return cells nodal coordinates  
     */
    inline std::vector< std::vector<double> > GetDomainCoord() const {return nodalpos_xi_domain_;};

    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    inline BlitzMat GetDomainCoordBlitz() const
    {
        const int numnode = this->NumNode();
        BlitzMat positions_xi_domain(3,numnode,blitz::ColumnMajorArray<2>());
        for (int inode=0; inode<numnode; inode++)
        {
            positions_xi_domain(0,inode) = this->GetDomainCoord()[inode][0];
            positions_xi_domain(1,inode) = this->GetDomainCoord()[inode][1];
            positions_xi_domain(2,inode) = this->GetDomainCoord()[inode][2];
        }    
        return positions_xi_domain;
    }
    
    /*!
     \brief returns the coordinates of the integration cell in physical coordinates

     \return cells nodal coordinates  
     */
    inline std::vector< std::vector<double> > GetPhysicalCoord(DRT::Element& ele) const
    {
        return ComputePhysicalCoordinates(ele);
    };
    
    vector<double> modifyGaussRule3D(
            const bool standard_integration,
            const double& cell_e0,
            const double& cell_e1,
            const double& cell_e2) const;
    
    virtual std::string Print() const;
    
    blitz::Array<double,1> GetPhysicalCenterPosition(DRT::Element& ele) const;

protected:
    
    //! coordinates of the nodes in the parent domain element coordinates xi_domain
    std::vector< std::vector<double> > nodalpos_xi_domain_;
    
private:
    
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
    BoundaryIntCell(
            const DRT::Element::DiscretizationType distype,
            const int xfemConditionLabel,
            const std::vector< std::vector<double> >& domainCoordinates,
            const std::vector< std::vector<double> >& boundaryCoordinates);

    /*!
     \brief Copy Constructor

     Makes a deep copy of a Element
     */
    BoundaryIntCell(
            const BoundaryIntCell& old);

    /*!
     \brief returns the shape of the integration cell

     \return distype_
     */
    inline DRT::Element::DiscretizationType Shape() const {return distype_;};
    
    /*!
     \brief returns number of nodes

     \return number of nodes
     */
    inline int NumNode() const {return DRT::UTILS::getNumberOfElementNodes(distype_);};
    
    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    inline std::vector< std::vector<double> > GetDomainCoord() const {return nodalpos_xi_domain_;};
    
    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    inline BlitzMat GetDomainCoordBlitz() const
    {
        const int numnode = this->NumNode();
        BlitzMat positions_xi_domain(3,numnode,blitz::ColumnMajorArray<2>());
        for (int inode=0; inode<numnode; inode++)
        {
            positions_xi_domain(0,inode) = this->GetDomainCoord()[inode][0];
            positions_xi_domain(1,inode) = this->GetDomainCoord()[inode][1];
            positions_xi_domain(2,inode) = this->GetDomainCoord()[inode][2];
        }    
        return positions_xi_domain;
    }
    
    /*!
     \brief returns the coordinates of the integration cell in physical coordinates

     \return cells nodal coordinates  
     */
    inline std::vector< std::vector<double> > GetPhysicalCoord(DRT::Element& ele) const
    {
        return ComputePhysicalCoordinates(ele);
    };
    
    /*!
     \brief returns the coordinates of the integration cell in surface element coordinates

     \return cells nodal coordinates  
     */
    inline std::vector< std::vector<double> > GetBoundaryCoord() const {return nodalpos_xi_boundary_;};
    
    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    inline BlitzMat GetBoundaryCoordBlitz() const
    {
        const int numnode = this->NumNode();
        BlitzMat positions_xi_domain(2,numnode,blitz::ColumnMajorArray<2>());
        for (int inode=0; inode<numnode; inode++)
        {
            positions_xi_domain(0,inode) = this->GetBoundaryCoord()[inode][0];
            positions_xi_domain(1,inode) = this->GetBoundaryCoord()[inode][1];
        }    
        return positions_xi_domain;
    }
    
    inline int GetSurfaceEleGid() const {return surface_ele_gid_;}
    
    std::string Print() const;

private:
    
    //! the boundaryIntCell must know to which xfemCondition it belongs
    int xfemConditionLabel_;
    
    //! the boundaryIntCell should know, to which cutterElement it belongs!?
    int surface_ele_gid_;
    
    //! coordinates of the nodes
    std::vector< std::vector<double> > nodalpos_xi_domain_;
    
    //! boundary coordinates of the nodes
    std::vector< std::vector<double> > nodalpos_xi_boundary_;
    
};

typedef std::vector<XFEM::BoundaryIntCell > BoundaryIntCells;
typedef std::set<XFEM::BoundaryIntCell > BoundaryIntCellSet;

}

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef CCADISCRET
