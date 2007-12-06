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
    inline int NumNode() const {return DRT::Utils::getNumberOfElementNodes(distype_);};
    
    /*!
     \brief returns the coordinates of the integration cell

     \return cells nodal coordinates  
     */
    inline std::vector< std::vector<double> > GetDomainCoord() const {return domainCoordinates_;};

    /*!
     \brief returns the coordinates of the integration cell in physical coordinates

     \return cells nodal coordinates  
     */
    inline std::vector< std::vector<double> > GetPhysicalCoord(DRT::Element& ele) const
    {
        return ComputePhysicalCoordinates(ele);
    };
    
    const vector<double> modifyGaussRule3D(
            const bool standard_integration,
            const double& cell_e0,
            const double& cell_e1,
            const double& cell_e2) const;
    
    virtual std::string Print() const;
    
    blitz::Array<double,1> GetCenterPosition(DRT::Element& ele) const;

protected:
    
    //! coordinates of the nodes in the parent volume element coordinates xi
    std::vector< std::vector<double> > domainCoordinates_;
    
private:
    
    void SetDefaultCoordinates(
            const DRT::Element::DiscretizationType distype);
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
            const std::vector< std::vector<double> > domainCoordinates,
            const std::vector< std::vector<double> > boundaryCoordinates);

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
    inline int NumNode() const {return DRT::Utils::getNumberOfElementNodes(distype_);};
    
    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    inline std::vector< std::vector<double> > GetDomainCoord() const {return domainCoordinates_;};
    
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
    inline std::vector< std::vector<double> > GetBoundaryCoord() const {return boundaryCoordinates_;};
    
    std::string Print() const;

private:
    
    //! coordinates of the nodes
    std::vector< std::vector<double> > domainCoordinates_;
    
    //! boundary coordinates of the nodes
    std::vector< std::vector<double> > boundaryCoordinates_;
    
};

typedef std::vector<XFEM::BoundaryIntCell > BoundaryIntCells;
typedef std::set<XFEM::BoundaryIntCell > BoundaryIntCellSet;

};

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef CCADISCRET
