/*----------------------------------------------------------------------*/
/*!
\file fluid3_stationary.H

\brief Internal implementation of Fluid3 element -- stationary formulation

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef FLUID3_STATIONARY_H
#define FLUID3_STATIONARY_H

#include "fluid3.H"
#include <blitz/array.h>


namespace DRT
{
namespace ELEMENTS
{

  /// Internal Fluid3 implementation -- stationary formulation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Fluid3 element in stationary formulation. 
    Additionally the method Sysmat() provides a clean and fast 
    element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Fluid3)
    from the mathematical contents (this class). Of course there are
    different implementations of the Fluid3 element, this is just one
    such implementation.

    The Fluid3 will allocate exactly one object of this class for all
    fluid elements with the same number of nodes in the mesh. This
    allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    For matrix calculation the blitz library is used.

    <h3>History</h3>
    Right now we do not read any stabilization parameters from the
    input file but have a fixed version with precompiled values.

    \author gjb
    \date 11/07
  */
  class Fluid3Stationary
  {
  public:
    /// Constructor with number of nodes
    explicit Fluid3Stationary(int iel);

    /// calculating element matrix and rhs
    /*!
      Calculate matrix and rhs for stationary problem formulation
      Here the whole thing is hidden.

      \param ele     		(i) the element those matrix is calculated
      \param evelnp  		(i) nodal velocities at n+1
      \param eprenp  		(i) nodal pressure
      \param estif   		(o) element matrix to calculate
      \param esv                (o) fine-scale subgrid-viscosity element matrix
      \param eforce  		(o) element rhs to calculate
      \param material 		(i) fluid material
      \param pseudotime     (i) current simulation pseudotime
      \param newton     	(i) full Newton or fixed-point-like
      \param pspg     		(i) boolean flag for stabilisation 
      \param supg     		(i) boolean flag for stabilisation
      \param vstab     		(i) boolean flag for stabilisation
      \param cstab     		(i) boolean flag for stabilisation
      \param cross     		(i) boolean flag for stabilisation
      \param reynolds     	(i) boolean flag for stabilisation
      \param Cs_fs  (i) Smagorinsky model param. for fine-scale subgrid viscosity
     */
    void Sysmat(Fluid3* ele,
                const blitz::Array<double,2>&            evelnp,
                const blitz::Array<double,1>&            eprenp,
                blitz::Array<double,2>&                  estif,
                blitz::Array<double,2>&                  esv,
                blitz::Array<double,1>&                  eforce,
                blitz::Array<double,1>&                  sugrvisc,
                struct _MATERIAL*                        material,
                double                                   pseudotime,
                bool                                     newton,
                int                                      fssgv,
                const enum Fluid3::StabilisationAction   pspg,
                const enum Fluid3::StabilisationAction   supg,
                const enum Fluid3::StabilisationAction   vstab,
                const enum Fluid3::StabilisationAction   cstab,
                const enum Fluid3::StabilisationAction   cross,
                const enum Fluid3::StabilisationAction   reynolds,
                double                                   Cs_fs);

  private:

    /*!
      \brief calculate stabilization parameter for stationary formulation
    */
    void CalTauStationary(
      Fluid3* ele,
      const blitz::Array<double,2>&           evelnp,
      const DRT::Element::DiscretizationType  distype,
      const double                            visc,
      const int                               fssgv,
      const double                            Cs_fs);

    /*!
     * \brief calculate the body force from the nodes condition
     */
    void BodyForce(Fluid3* ele, const double pseudotime);

    /*!
     * \brief calculate second global derivatives w.r.t x/y/z at point r,s,t
     */
    void gder2(Fluid3* ele);

  private:
    int iel_;                          /// number of nodes
    double vart_;                      /// fine-scale subgrid viscosity
    blitz::Array<double,2> xyze_;      /// nodal position array
    blitz::Array<double,2> edeadng_;   /// nodal dead load 
    blitz::Array<double,1> funct_;     /// nodal shape functions
    blitz::Array<double,2> deriv_;     /// shape derivatives with respect to \xi
    blitz::Array<double,2> deriv2_;    /// second shape derivatives with respect to \xi
    blitz::Array<double,2> xjm_;       /// transposed of the jacobian matrix d x / d \xi
    blitz::Array<double,2> xji_;       /// inverse of the jacobian matrix
    blitz::Array<double,2> vderxy_;    /// velocity vector derivative with respect to x
    blitz::Array<double,1> pderxy_;    /// pressure derivative with respect to x
    blitz::Array<double,2> vderxy2_;   /// second derivatives of the velocity vector
    blitz::Array<double,2> derxy_;     /// shape derivative with respect to x
    blitz::Array<double,2> derxy2_;    /// second shape derivative with respect to x
    blitz::Array<double,1> bodyforce_; /// volumetric load at gausspoint
    blitz::Array<double,1> velino_;
    blitz::Array<double,1> velint_;    /// velocity
    blitz::Array<double,1> gradp_;     /// pressure gradient
    blitz::Array<double,1> tau_;       /// tau_Mu, tau_Mp, tau_C
    blitz::Array<double,3> viscs2_;    /// viscous term including 2nd derivatives
    blitz::Array<double,1> conv_c_;    /// linearisation of convect, convective part
    blitz::Array<double,1> conv_g_;    /// linearisation of convect, grid part
    blitz::Array<double,3> conv_r_;    /// linearisation of convect, reactive part
    blitz::Array<double,1> rhsint_;    /// total right hand side terms at int.-point
    blitz::Array<double,1> conv_old_;
    blitz::Array<double,1> visc_old_;
    blitz::Array<double,1> res_old_;
    blitz::Array<double,1> conv_resM_;
    blitz::Array<double,2> xder2_;
    blitz::Array<int,1> numepn_;
  };

}
}

#endif

#endif
#endif
