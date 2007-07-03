/*!----------------------------------------------------------------------
\file integrationcell.H

\brief integration cell

<pre>
Maintainer: 
</pre>

*----------------------------------------------------------------------*/


#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef INTEGRATIONCELL_H
#define INTEGRATIONCELL_H

#include <vector>

class Integrationcell
{
        
    public:  
    
    /*!
    \brief Standard Constructor

    \param id               : A unique global id
    \param coordinates      : coordinates of the integrationcell
    */
    Integrationcell(int id, std::vector< std::vector<double> > coordinates);
    
    
    /*!
    \brief Copy Constructor

    Makes a deep copy of a Element

    */
    Integrationcell(const Integrationcell& old);

    
    /*!
    \brief returns the coordinates of the integration cell

    \return returns the coordinates 

    */
    std::vector< std::vector<double> > GetCoord();
    
     
    private:
    
    //! The global integration cell id
    int                                             id_;
    
     //! coordinates of the nodes
    std::vector< std::vector<double> >                       coordinates_;
        
};

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM



