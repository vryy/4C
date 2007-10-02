/*!
 \file integrationcell.H

 \brief integration cell

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
#include <string>

namespace XFEM
{

/*
 * \brief An integrationcell is used for specialized element integration routines 
 */
class Integrationcell
{

public:

    /*!
     \brief Standard Constructor

     \param id               : A unique global id
     \param coordinates      : coordinates of the integrationcell
     */
    Integrationcell(
        const int id,
        const std::vector< std::vector<double> > coordinates);

    /*!
     \brief Copy Constructor

     Makes a deep copy of a Element
     */
    Integrationcell(const Integrationcell& old);

    /*!
     \brief returns the coordinates of the integration cell

     \return cells nodal coordinates  
     */
    std::vector< std::vector<double> > GetCoord() const;
    
    std::string Print() const;

private:

    //! The global integration cell id
    int id_;

    //! coordinates of the nodes
    std::vector< std::vector<double> > coordinates_;

};

typedef std::vector<XFEM::Integrationcell > IntCells;

};

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef CCADISCRET
