/*----------------------------------------------------------------------*/
/*!
\file fluid2_stationary.H

\brief Internal implementation of Fluid2 element (stationary formulation)

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID2
#ifdef CCADISCRET

#ifndef FLUID2_STATIONARY_H
#define FLUID2_STATIONARY_H

#include "fluid2.H"
#include <blitz/array.h>


namespace DRT
{
namespace ELEMENTS
{

  /// Internal Fluid2 implementation -- stationary formulation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Fluid2 element in stationary formulation.
    Additionally the method Sysmat() provides a clean and fast
    element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Fluid2)
    from the mathematical contents (this class). Of course there are
    different implementations of the Fluid2 element, this is just one
    such implementation.

    The Fluid2 will allocate exactly one object of this class for all
    fluid elements with the same number of nodes in the mesh. This
    allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    For matrix calculation the blitz library is used.

    \author vg
    \date 08/08
  */
  class Fluid2Stationary
  {
  public:
    /// Constructor with number of nodes
    explicit Fluid2Stationary(int iel);

    /// calculating element matrix and rhs
    /*!
      Calculate matrix and rhs for stationary problem formulation
      Here the whole thing is hidden.

      \param ele     		  (i) the element those matrix is calculated
      \param evelnp  		  (i) nodal velocities at n+1
      \param fsevelnp         (i) fine-scale nodal velocities at n+1
      \param eprenp  		  (i) nodal pressure
      \param edensnp          (i) nodal density
      \param estif   		  (o) element matrix to calculate
      \param eforce  		  (o) element rhs to calculate
      \param material 		  (i) fluid material
      \param pseudotime       (i) current simulation pseudotime
      \param newton     	  (i) full Newton or fixed-point-like
      \param loma             (i) boolean flag for potential low-Mach-number solver
      \param higher_order_ele (i) keep or drop second derivatives
      \param fssgv            (i) flag for type of fine-scale subgrid viscosity
      \param pspg     		  (i) boolean flag for stabilisation
      \param supg     		  (i) boolean flag for stabilisation
      \param vstab     		  (i) boolean flag for stabilisation
      \param cstab     		  (i) boolean flag for stabilisation
      \param cross     		  (i) boolean flag for stabilisation
      \param reynolds     	  (i) boolean flag for stabilisation
      \param Cs               (i) Smagorinsky model parameter
     */
    void Sysmat(Fluid2*                                  ele,
                const blitz::Array<double,2>&            evelnp,
                const blitz::Array<double,2>&            fsevelnp,
                const blitz::Array<double,1>&            eprenp,
                const blitz::Array<double,1>&            edensnp,
                blitz::Array<double,2>&                  estif,
                blitz::Array<double,1>&                  eforce,
                struct _MATERIAL*                        material,
                double                                   pseudotime,
                bool                                     newton,
                bool                                     loma,
                const bool                               higher_order_ele,
                const enum Fluid2::StabilisationAction   fssgv,
                const enum Fluid2::StabilisationAction   pspg,
                const enum Fluid2::StabilisationAction   supg,
                const enum Fluid2::StabilisationAction   vstab,
                const enum Fluid2::StabilisationAction   cstab,
                const enum Fluid2::StabilisationAction   cross,
                const enum Fluid2::StabilisationAction   reynolds,
                const double                             Cs);

    /// Internal implementation class for specific fluid element
    static Fluid2Stationary* StationaryImpl(DRT::ELEMENTS::Fluid2* f2);

  private:

    /*!
      \brief calculate stabilization parameter for stationary formulation
    */
    void CalTauStationary(Fluid2*                                 ele,
                          const blitz::Array<double,2>&           evelnp,
                          const blitz::Array<double,2>&           fsevelnp,
                          const blitz::Array<double,1>&           edensnp,
                          const DRT::Element::DiscretizationType  distype,
                          const double                            visc,
                          const enum Fluid2::StabilisationAction  fssgv,
                          const double                            Cs);

    /*!
     * \brief calculate the body force from the nodes condition
     */
    void BodyForce(Fluid2*      ele,
                   const double pseudotime);

  private:
    int iel_;                          /// number of nodes
    double vart_;                      /// fine-scale subgrid viscosity
    blitz::Array<double,2> xyze_;      /// nodal position array
    blitz::Array<double,2> edeadng_;   /// nodal dead load
    blitz::Array<double,1> funct_;     /// nodal shape functions
    blitz::Array<double,1> densfunct_; /// density-weighted nodal shape functions
    blitz::Array<double,1> functdens_; /// nodal shape functions divided by density
    blitz::Array<double,2> deriv_;     /// shape derivatives with respect to \xi
    blitz::Array<double,2> deriv2_;    /// second shape derivatives with respect to \xi
    blitz::Array<double,2> xjm_;       /// transposed of the jacobian matrix d x / d \xi
    blitz::Array<double,2> xji_;       /// inverse of the jacobian matrix
    blitz::Array<double,2> vderxy_;    /// velocity vector derivative
    blitz::Array<double,2> mderxy_;    /// momentum derivative
    blitz::Array<double,2> fsvderxy_;  /// fine-scale velocity vector derivative
    blitz::Array<double,2> vderxy2_;   /// second derivatives of the velocity vector
    blitz::Array<double,2> derxy_;     /// shape derivative with respect to x
    blitz::Array<double,2> densderxy_; /// density-weighted shape derivative w.r.t. x
    blitz::Array<double,2> derxy2_;    /// second shape derivative with respect to x
    blitz::Array<double,1> bodyforce_; /// volumetric load at gausspoint
    //blitz::Array<double,1> velino_;
    blitz::Array<double,1> velint_;    /// velocity
    blitz::Array<double,1> fsvelint_;  /// fine-scale velocity
    blitz::Array<double,1> gradp_;     /// pressure gradient
    blitz::Array<double,1> tau_;       /// tau_Mu, tau_Mp, tau_C
    blitz::Array<double,3> viscs2_;    /// viscous term including 2nd derivatives
    blitz::Array<double,1> conv_c_;    /// linearisation of convect, convective part
    double mdiv_;                      /// momentum divergence
    double vdiv_;                      /// velocity divergence
    blitz::Array<double,1> rhsmom_;    /// momentum equation right hand side
    blitz::Array<double,1> conv_old_;
    blitz::Array<double,1> visc_old_;
    blitz::Array<double,1> res_old_;
    blitz::Array<double,1> conv_resM_;
    blitz::Array<double,2> xder2_;
  };

}
}

#endif

#endif
#endif
