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
  * \brief A space-time boundary element
  *
  * If a boundary travels through time, we can for each time step connect
  * the boundary at its new position with one or several of its old positions.
  * This leads to a volumetric element. For example, for one time step,
  * a quad4 surface element generates a hex8 space-time element,
  * while a tri3 generates a wedge element.
  */
  class SpaceTimeBoundaryCell
  {

  public:

    //! Standard Constructor
    explicit SpaceTimeBoundaryCell(
        const int                              bele_id,  ///< id of boundary element we want to track
        const DRT::Element::DiscretizationType surf_distype,  ///< shape of the boundary element
        const LINALG::SerialDenseMatrix&       posnp,    ///< nodal positions at n+1
        const LINALG::SerialDenseMatrix&       posn      ///< nodal positions at n
    );

    //! Copy Constructor
    SpaceTimeBoundaryCell(
            const SpaceTimeBoundaryCell& old);
    
    //! Default Constructor
    explicit SpaceTimeBoundaryCell();
    
    //! operator
    SpaceTimeBoundaryCell& operator = (const SpaceTimeBoundaryCell& old);
    
    //! return position array
    const LINALG::SerialDenseMatrix& get_xyzt() const {return xyzt_;};
    
    //! get global id of boundary element for which this space-time cell is constructed
    int getBeleId() const {return bele_id_;};
    
    //! get global id of boundary element for which this space-time cell is constructed
    DRT::Element::DiscretizationType Shape() const {return DRT::Element::hex8;};
    
    //! return string representation
    std::string toString() const;
    
  private:
     
    /// id of boundary element we want to track
    int bele_id_;
    
    /// shape of the boundary element
    DRT::Element::DiscretizationType surf_distype_;
    
    /// number of timesteps this cell covers
    int num_timestep_;
    
    /// nodal positions at (x,y,z,tnp,tn)
    LINALG::SerialDenseMatrix xyzt_;
    
    //! from surface positions at two time steps, create prism element (e.g. hex)
    LINALG::SerialDenseMatrix getLinearPositionArray(
        const LINALG::SerialDenseMatrix&      posnp, ///< nodal positions at n+1
        const LINALG::SerialDenseMatrix&      posn   ///< nodal positions at n
    ) const;
    
  };
}

#endif  // #ifndef INTEGRATIONCELL_H
#endif  // #ifdef CCADISCRET
