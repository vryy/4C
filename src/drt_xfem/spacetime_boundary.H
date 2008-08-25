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
            const int           bele_id,                ///< boundary element id we want to track
            const BlitzMat&     posnp,                  ///< nodal positions at n+1
            const BlitzMat&     posn                    ///< nodal positions at n
            );

    //! Copy Constructor
    SpaceTimeBoundaryCell(
            const SpaceTimeBoundaryCell& old);
    
    explicit SpaceTimeBoundaryCell();
    
    /*!  = operator
     * create true copy
     */
    SpaceTimeBoundaryCell& operator = (const SpaceTimeBoundaryCell& old)
    {
      SpaceTimeBoundaryCell cell(old);
      SpaceTimeBoundaryCell& ref = cell;
      return ref;
    };
    
    BlitzMat get_xyzt() const {return xyzt_;};
    
    //! return string representation
    std::string toString() const;

private:
  
  const int           bele_id_;           ///< boundary element id we want to track
  const BlitzMat      posnp_;             ///< nodal positions at n+1
  const BlitzMat      posn_;              ///< nodal positions at n
  const BlitzMat      xyzt_;              ///< nodal positions at (x,y,z,tnp,tn)
  
  BlitzMat getLinearPositionArray(
      const BlitzMat&      posnp,                 ///< nodal positions at n+1
      const BlitzMat&      posn                   ///< nodal positions at n
      ) const;
    
};
}

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef CCADISCRET
