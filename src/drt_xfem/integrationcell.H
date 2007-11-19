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
    DRT::Element::DiscretizationType Shape() const {return distype_;};
    
    /*!
     \brief returns the coordinates of the integration cell

     \return cells nodal coordinates  
     */
    std::vector< std::vector<double> > GetDomainCoord() const;
    
    virtual std::string Print() const;

protected:

    //! shape
    DRT::Element::DiscretizationType distype_;

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
            const std::vector< std::vector<double> > domainCoordinates);

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
    DRT::Element::DiscretizationType Shape() const {return distype_;};
    
    /*!
     \brief returns the coordinates of the integration cell

     \return cells nodal coordinates  
     */
    std::vector< std::vector<double> > GetDomainCoord() const;
    
    std::string Print() const;
    
    double VolumeRatio(
            const DRT::Element::DiscretizationType  parentdistype) const;
    
    double Volume() const {return volume_;};

protected:
    
    //! coordinates of the nodes
    std::vector< std::vector<double> > domainCoordinates_;
    
private:
    
    // volume with respect to the parent element coordinates
    double  volume_;
    
    void SetDefaultCoordinates(
            const DRT::Element::DiscretizationType distype);
    
    double GetVolume() const;
    
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
    DRT::Element::DiscretizationType Shape() const {return distype_;};
    
    /*!
     \brief returns the coordinates of the integration cell in domain element coordinates

     \return cells nodal coordinates  
     */
    std::vector< std::vector<double> > GetDomainCoord() const;
    
    /*!
     \brief returns the coordinates of the integration cell in surface element coordinates

     \return cells nodal coordinates  
     */
    std::vector< std::vector<double> > GetBoundaryCoord() const;
    
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
