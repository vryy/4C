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
            const int                            bele_id,  ///< boundary element id we want to track
            const LINALG::SerialDenseMatrix&     posnp,    ///< nodal positions at n+1
            const LINALG::SerialDenseMatrix&     posn      ///< nodal positions at n
            );

//    //! Copy Constructor
//    SpaceTimeBoundaryCell(
//            const SpaceTimeBoundaryCell& old);
    
    //! Default Constructor
    explicit SpaceTimeBoundaryCell();
    
//    /*!  = operator
//     * create true copy
//     */
//    SpaceTimeBoundaryCell& operator = (const SpaceTimeBoundaryCell& old)
//    {
//      SpaceTimeBoundaryCell cell(old.bele_id_, posnp_, posn_);
//      SpaceTimeBoundaryCell& ref = cell;
//      return ref;
//    };
    
    //! return position array
    LINALG::SerialDenseMatrix get_xyzt() const {return xyzt_;};
    
    //! get global id of boundary element for which this space-time cell is constructed
    int getBeleId() const {return bele_id_;};
    
    //! return string representation
    std::string toString() const;

private:
  
  int           bele_id_;           ///< boundary element id we want to track
  LINALG::SerialDenseMatrix      posnp_;             ///< nodal positions at n+1
  LINALG::SerialDenseMatrix      posn_;              ///< nodal positions at n
  LINALG::SerialDenseMatrix      xyzt_;              ///< nodal positions at (x,y,z,tnp,tn)
  
  //! from surface positions at two time steps, create prism element (e.g. hex)
  LINALG::SerialDenseMatrix getLinearPositionArray(
      const LINALG::SerialDenseMatrix&      posnp,                 ///< nodal positions at n+1
      const LINALG::SerialDenseMatrix&      posn                   ///< nodal positions at n
      ) const;
    
};
}

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef CCADISCRET
