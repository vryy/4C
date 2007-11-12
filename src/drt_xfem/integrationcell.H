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

     \param coordinates      : coordinates of the integrationcell
     */
    IntCell(
        const std::vector< std::vector<double> > domainCoordinates);

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
     \brief returns the coordinates of the integration cell

     \return cells nodal coordinates  
     */
    std::vector< std::vector<double> > GetDomainCoord() const;
    
    virtual std::string Print() const;

protected:

    //! coordinates of the nodes
    std::vector< std::vector<double> > domainCoordinates_;

};

typedef std::vector<XFEM::IntCell > IntCells;
typedef std::set<XFEM::IntCell > IntCellSet;



/*
 * \brief An domain integration cell is used for integrating the XFEM elements domain 
 */
class DomainIntCell : IntCell
{

public:

    /*!
     \brief Standard Constructor

     \param coordinates      : coordinates of the integrationcell
     */
    DomainIntCell(
        const std::vector< std::vector<double> > coordinates);

    /*!
     \brief Copy Constructor

     Makes a deep copy of a Element
     */
    DomainIntCell(const DomainIntCell& old);

    /*!
     \brief returns the coordinates of the integration cell

     \return cells nodal coordinates  
     */
    std::vector< std::vector<double> > GetDomainCoord() const;
    
    std::string Print() const;

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
     */
    BoundaryIntCell(
        const std::vector< std::vector<double> > coordinates);

    /*!
     \brief Copy Constructor

     Makes a deep copy of a Element
     */
    BoundaryIntCell(const BoundaryIntCell& old);

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
    
    //! boundary coordinates of the nodes
    std::vector< std::vector<double> > boundaryCoordinates_;
    
};

typedef std::vector<XFEM::BoundaryIntCell > BoundaryIntCells;
typedef std::set<XFEM::BoundaryIntCell > BoundaryIntCellSet;

};

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef CCADISCRET
