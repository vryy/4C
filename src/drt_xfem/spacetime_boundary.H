/*!
 \file spacetime_boundary.H

 \brief tools for proper space-time treatment of moving XFEM interfaces

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */

#ifdef CCADISCRET
#ifndef SPACETIME_BOUNDARY_H
#define SPACETIME_BOUNDARY_H

#include <set>
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


namespace XFEM
{

/*!
 * \brief A space-time boundary element (whatever this may mean) 
 */
class SpaceTimeBoundaryCell
{

public:

    //! Standard Constructor
    explicit SpaceTimeBoundaryCell(
            const DRT::Element* boundaryele,            ///< boundary element we want to track
            const BlitzMat&     posnp,                  ///< nodal positions at n+1
            const BlitzMat&     posn                    ///< nodal positions at n
            );

    //! Copy Constructor
    SpaceTimeBoundaryCell(
            const SpaceTimeBoundaryCell& old);
    
    /*!  = operator
     * create true copy
     */
    SpaceTimeBoundaryCell& operator = (const SpaceTimeBoundaryCell& old)
    {
      SpaceTimeBoundaryCell cell(old);
      SpaceTimeBoundaryCell& ref = cell;
      return ref;
    };
    
    //! return string representation
    std::string toString() const;

private:
  
  const DRT::Element* boundaryele_;           ///< boundary element we want to track
  const BlitzMat      posnp_;                 ///< nodal positions at n+1
  const BlitzMat      posn_;                  ///< nodal positions at n
  
    
  //! hide default constructor
  explicit SpaceTimeBoundaryCell();
    
};
}

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef CCADISCRET
