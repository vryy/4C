/*----------------------------------------------------------------------*/
/*!
\file fluid3_impl.H

\brief Internal implementation of Fluid3 element

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef FLUID3_IMPL_H
#define FLUID3_IMPL_H

#include "fluid3.H"
#include <blitz/array.h>

namespace DRT
{
namespace ELEMENTS
{

  /// Internal Fluid3 implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Fluid3 element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

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

    The implementation here is based on the work by Christiane
    Foerster. The core element code has been linearized using MuPAD.
    The reaction terms have been removed from the lhs.

    Right now we do not read any stabilization parameters from the
    input file but have a fixed version with precompiled values.

    \author u.kue
    \date 07/07
  */
  class Fluid3Impl
  {
  public:
    /// Constructor with number of nodes
    explicit Fluid3Impl(int iel);

    /// calculating element matrix and rhs
    /*!
      Calculate matrix and rhs. Here the whole thing is hidden.

      \note Stationary calculations are not (yet) supported.

      \param ele              (i) the element those matrix is calculated
      \param evelnp           (i) nodal velocities at n+1
      \param csevelnp         (i) ???????????????????????
      \param fsevelnp         (i) fine-scale nodal velocities at n+1
      \param cseconvnp        (i) ???????????????????????
      \param eprenp           (i) nodal pressure
      \param evhist           (i) rhs from beginning of time step
      \param edispnp          (i) nodal displacements (on moving mesh)
      \param egridv           (i) grid velocity (on moving mesh)
      \param estif            (o) element matrix to calculate
      \param emesh            (o) linearization wrt mesh motion
      \param eforce           (o) element rhs to calculate
      \param material         (i) fluid material
      \param time             (i) current simulation time
      \param dt               (i) timestep
      \param timefac          (i) time discretization factor
      \param newton           (i) boolean flag for linearisation
      \param fssgv            (i) flag for type of fine-scale subgrid viscosity
      \param higher_order_ele (i) keep or drop second derivatives
      \param fssgv            (i) ?????????????????????????
      \param pspg             (i) boolean flag for stabilisation
      \param supg             (i) boolean flag for stabilisation
      \param vstab            (i) boolean flag for stabilisation
      \param cstab            (i) boolean flag for stabilisation
      \param cross            (i) boolean flag for stabilisation
      \param reynolds         (i) boolean flag for stabilisation
      \param whichtau         (i) selecting available tau definition
      \param turb_mod_action  (i) selecting turbulence model (none, Smagorisky,
                                  dynamic Smagorinsky, Smagorinsky with van Driest
                                  damping for channel flows)
      \param Cs               (i) Smagorinsky model parameter
      \param Cs_delta_sq      (i) Model parameter computed by dynamic Smagorinsky
                                  approach (Cs*h*h)
      \param visceff          (i) nu+nu_turbulent
      \param l_tau            (i) viscous length scale, required for van driest
                                  damping function and defined on input
      */
    void Sysmat(Fluid3*                                  ele,
                const blitz::Array<double,2>&            evelnp,
                const blitz::Array<double,2>&            csevelnp,
                const blitz::Array<double,2>&            fsevelnp,
                const blitz::Array<double,2>&            cseconvnp,
                const blitz::Array<double,1>&            eprenp,
                const blitz::Array<double,2>&            evhist,
                const blitz::Array<double,2>&            edispnp,
                const blitz::Array<double,2>&            egridv,
                blitz::Array<double,2>&                  estif,
                blitz::Array<double,2>&                  emesh,
                blitz::Array<double,1>&                  eforce,
                struct _MATERIAL*                        material,
                double                                   time,
                double                                   dt,
                double                                   timefac,
                bool                                     newton,
                const bool                               higher_order_ele,
                const enum Fluid3::StabilisationAction   fssgv,
                const enum Fluid3::StabilisationAction   pspg,
                const enum Fluid3::StabilisationAction   supg,
                const enum Fluid3::StabilisationAction   vstab,
                const enum Fluid3::StabilisationAction   cstab,
                const enum Fluid3::StabilisationAction   cross,
                const enum Fluid3::StabilisationAction   reynolds,
                const enum Fluid3::TauType               whichtau,
                const enum Fluid3::TurbModelAction       turb_mod_action,
                double&                                  Cs,
                double&                                  Cs_delta_sq,
                double&                                  visceff,
                double&                                  l_tau);


    /// Internal implementation class for specific fluid element
    static Fluid3Impl* Impl(DRT::ELEMENTS::Fluid3* f3);

  private:

    /*!
      \brief calculate stabilization parameter

      o Standard:
           tauMp: Barrenechea, Valentin
           tauMu: Franca, Valentin
           tauC : Wall

      o Bazilevs: bazilevs
           tauM : Bazilevs, Whiting
           tauC : Bazilevs

      o Codina: codina
           tauMp: Barrenechea, Valentin
           tauMu: Franca, Valentin
           tauC : Codina

      \param evelnp           (i) nodal velocities at n+1
      \param fsevelnp         (i) fine-scale nodal velocities at n+1
      \param distype          (i) hex8, tet4 etc.
      \param whichtau         (i) selecting available tau definition
      \param visc             (i) material viscosity
      \param timefac          (i) time discretization factor
      \param dt               (i) timestep
      \param turb_mod_action  (i) selecting turbulence model (none, Smagorisky,
                                  dynamic Smagorinsky, Smagorinsky with van Driest
                                  damping for channel flows)
      \param Cs               (i) Smagorinsky model parameter
      \param Cs_delta_sq      (i) Model parameter computed by dynamic Smagorinsky
                                  approach (Cs*h*h)
      \param viscturb         (o) kinematic viscosity turbulence
      \param visceff          (i) nu+nu_turbulent
      \param l_tau            (i) viscous length scale, required for van driest
                                  damping function and defined on input
      \param fssgv            (i) ?????????????????????????

    */
    void Caltau(
      Fluid3* ele,
      const blitz::Array<double,2>&           evelnp,
      const blitz::Array<double,2>&           fsevelnp,
      const DRT::Element::DiscretizationType  distype,
      const enum Fluid3::TauType              whichtau,
      struct _MATERIAL*                       material,
      double&                                 visc,
      const double                            timefac,
      const double                            dt,
      const enum Fluid3::TurbModelAction      turb_mod_action,
      double&                                 Cs,
      double&                                 Cs_delta_sq,
      double&								  viscturb,
      double&                                 visceff,
      double&                                 l_tau,
      const enum Fluid3::StabilisationAction  fssgv);


    /*!
         \brief calculates viscosity at every gaussian point

  		 \param material		 (i) material
         \param visc             (o) kinematic material viscosity
    */
    void CalVisc_AtGaussianPoint(
      const struct _MATERIAL*                 material,
      double&                           	  visc);

    /*!
     * \brief calculate the body force from the nodes condition
     */
    void BodyForce(Fluid3* ele, const double time, struct _MATERIAL* material);

    /*!
     * \brief calculate second global derivatives w.r.t x/y/z at point r,s,t
     */
    void gder2(Fluid3* ele);

  private:
    //! number of nodes
    int iel_;
    //! ????
    double vart_;
    //! node coordinates
    blitz::Array<double,2> xyze_;
    //! dead load in nodes
    blitz::Array<double,2> edeadng_;
    //! array for shape functions
    blitz::Array<double,1> funct_;
    //! array for shape function derivatives w.r.t r,s,t
    blitz::Array<double,2> deriv_;
    //! array for second derivatives of shape function w.r.t r,s,t
    blitz::Array<double,2> deriv2_;
    //! transposed jacobian "dx/ds"
    blitz::Array<double,2> xjm_;
    //! inverse of transposed jacobian "ds/dx"
    blitz::Array<double,2> xji_;
    //! global velocity derivatives in gausspoint w.r.t x,y,z
    blitz::Array<double,2> vderxy_;
    //! ???????
    blitz::Array<double,2> csvderxy_;
    //! ???????
    blitz::Array<double,2> fsvderxy_;
    //! global second velocity derivatives in gausspoint w.r.t x,y,z
    blitz::Array<double,2> vderxy2_;
    //! global derivatives of shape functions w.r.t x,y,z
    blitz::Array<double,2> derxy_;
    //! global second derivatives of shape functions w.r.t x,y,z
    blitz::Array<double,2> derxy2_;
    //! bodyforce in gausspoint
    blitz::Array<double,1> bodyforce_;
    //! vector containing all values from previous timelevel n
    blitz::Array<double,1> histvec_;
    //! direction of flow (normed velocity vector)
    blitz::Array<double,1> velino_;
    //! velocity vector in gausspoint
    blitz::Array<double,1> velint_;
    //! ??????????
    blitz::Array<double,1> csvelint_;
    //! ??????????
    blitz::Array<double,1> fsvelint_;
    //! ??????????
    blitz::Array<double,1> csconvint_;
    //! ale convective velocity in gausspoint, c=u-u_G
    blitz::Array<double,1> convvelint_;
    //! pressure gradient in gausspoint
    blitz::Array<double,1> gradp_;
    //! the stabilisation parameters
    blitz::Array<double,1> tau_;
    //! viscous term including 2nd derivatives
    blitz::Array<double,3> viscs2_;
    //! linearisation of convect, convective part
    blitz::Array<double,1> conv_c_;
    //! total right hand side terms at int.-point
    blitz::Array<double,1> rhsint_;
    //! (u_old*nabla)u_old
    blitz::Array<double,1> conv_old_;
    //! convective part of scale-similarity term
    blitz::Array<double,1> conv_s_;
    //! div epsilon(u_old)
    blitz::Array<double,1> visc_old_;
    //! the old residual
    blitz::Array<double,1> res_old_;
    //! convection by subgrid-scale velocity
    blitz::Array<double,1> conv_resM_;
    //! 2nd derivatives of coord.-functions w.r.t r,s,t
    blitz::Array<double,2> xder2_;

    blitz::Array<double,2> vderiv_;
    blitz::Array<double,4> derxjm_;
  };

}
}

#endif

#endif
#endif
