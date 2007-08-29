/*!
 \file integrationcell.H

 \brief integration cell

 \author Axel Gerstenberger,
 gerstenberger@lnm.mw.tum.de,
 http://www.lnm.mw.tum.de,
 089 - 289-15236
 */

#ifdef XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef INTEGRATIONCELL_H
#define INTEGRATIONCELL_H

#include <vector>

namespace Xfem
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
    std::vector< std::vector<double> > GetCoord();

private:

    //! The global integration cell id
    int id_;

    //! coordinates of the nodes
    std::vector< std::vector<double> > coordinates_;

};
};

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef XFEM
