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
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/linalg_serialdensematrix.H"

namespace DRT
{
namespace ELEMENTS
{

  /// Interface base class for Fluid3Impl
  /*!
    This class exists to provide a common interface for all template
    versions of Fluid3Impl. The only function
    this class actually defines is Impl, which returns a pointer to
    the appropriate version of Fluid3Impl.
   */
  class Fluid3ImplInterface
  {
  public:
    /// Empty constructor
    Fluid3ImplInterface() {}
    /// Empty destructor
    virtual ~Fluid3ImplInterface() {}
    /// Evaluate the element
    /*!
      This class does not provide a definition for this function, it
      must be defined in Fluid3Impl.
     */
    virtual int Evaluate(Fluid3*                   ele,
                         ParameterList&            params,
                         DRT::Discretization&      discretization,
                         vector<int>&              lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         RefCountPtr<MAT::Material> mat) = 0;

    /// Internal implementation class for fluid element
    static Fluid3ImplInterface* Impl(DRT::ELEMENTS::Fluid3* f3);

  };

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

    <h3>History</h3>

    The implementation here is based on the work by Christiane
    Foerster. The core element code has been linearized using MuPAD.
    The reaction terms have been removed from the lhs.

    Right now we do not read any stabilization parameters from the
    input file but have a fixed version with precompiled values.

    \author u.kue
    \date 07/07
  */
  template<DRT::Element::DiscretizationType distype>
  class Fluid3Impl: public Fluid3ImplInterface
  {
  public:
    /// Constructor
    Fluid3Impl();

    //! number of nodes
    static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    /*!
      \brief evaluate function for Fluid3 element
     */
    virtual int Evaluate(Fluid3*                    ele,
                         ParameterList&             params,
                         DRT::Discretization&       discretization,
                         vector<int>&               lm,
                         Epetra_SerialDenseMatrix&  elemat1_epetra,
                         Epetra_SerialDenseMatrix&  elemat2_epetra,
                         Epetra_SerialDenseVector&  elevec1_epetra,
                         Epetra_SerialDenseVector&  elevec2_epetra,
                         Epetra_SerialDenseVector&  elevec3_epetra,
                         RefCountPtr<MAT::Material> mat);

    /*!
      \brief calculate element matrix and rhs

      \param ele              (i) the element those matrix is calculated
      \param evelaf           (i) nodal velocities at n+alpha_F/n+1
      \param eveln            (i) nodal velocities at n
      \param fsevelaf         (i) fine-scale nodal velocities at n+alpha_F/n+1
      \param epreaf           (i) nodal pressure at n+alpha_F/n+1
      \param eaccam           (i) nodal accelerations at n+alpha_M
      \param escaaf           (i) nodal scalar at n+alpha_F/n+1
      \param escaam           (i) nodal scalar at n+alpha_M/n
      \param escadtam         (i) nodal scalar derivatives at n+alpha_M/n+1
      \param emhist           (i) time rhs for momentum equation
      \param edispnp          (i) nodal displacements (on moving mesh)
      \param egridv           (i) grid velocity (on moving mesh)
      \param estif            (o) element matrix to calculate
      \param emesh            (o) linearization wrt mesh motion
      \param eforce           (o) element rhs to calculate
      \param sgvelvisc        (o) subgrid-velocity/viscosity vector
      \param material         (i) fluid material
      \param time             (i) current simulation time
      \param timefac          (i) time discretization factor
      \param newton           (i) boolean flag for linearisation
      \param loma             (i) boolean flag for potential low-Mach-number solver
      \param conservative     (i) boolean flag for conservative form
      \param sgvel            (i) boolean flag for subgrid-scale velocity
      \param is_genalpha      (i) boolean flag for generalized-alpha time integration
      \param higher_order_ele (i) keep or drop second derivatives
      \param fssgv            (i) flag for type of fine-scale subgrid viscosity
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
      \param l_tau            (i) viscous length scale, required for van driest
                                  damping function and defined on input
      */
    void Sysmat(Fluid3*                                  ele,
                const LINALG::Matrix<3,iel>&             evelaf,
                const LINALG::Matrix<3,iel>&             eveln,
                const LINALG::Matrix<3,iel>&             fsevelaf,
                const LINALG::Matrix<iel,1>&             epreaf,
                const LINALG::Matrix<3,iel>&             eaccam,
                const LINALG::Matrix<iel,1>&             escaaf,
                const LINALG::Matrix<iel,1>&             escaam,
                const LINALG::Matrix<iel,1>&             escadtam,
                const LINALG::Matrix<3,iel>&             emhist,
                const LINALG::Matrix<3,iel>&             edispnp,
                const LINALG::Matrix<3,iel>&             egridv,
                LINALG::Matrix<4*iel,4*iel>&             estif,
                LINALG::Matrix<4*iel,4*iel>&             emesh,
                LINALG::Matrix<4*iel,    1>&             eforce,
                LINALG::Matrix<4*iel,    1>&             sgvelvisc,
                Teuchos::RCP<const MAT::Material>        material,
                double                                   time,
                double                                   timefac,
                const bool                               newton,
                const bool                               loma,
                const bool                               conservative,
                const bool                               sgvel,
                const bool                               is_genalpha,
                const bool                               higher_order_ele,
                const bool                               tau_gp,
                const bool                               mat_gp,
                const enum Fluid3::FineSubgridVisc       fssgv,
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
                double&                                  l_tau);


  private:

    /*!
     * \brief calculate body force from node condition
     */
    void BodyForce(Fluid3*      ele,
                   const double time);

    /*!
     * \brief evaluate shape functions and their derivatives at element center
     */
    void EvalShapeFuncAndDerivsAtEleCenter(const int  eleid,
                                           const bool higher_order_ele);

    /*!
     * \brief evaluate shape functions and their derivatives at integration point
     */
    void EvalShapeFuncAndDerivsAtIntPoint(
        const DRT::UTILS::IntegrationPoints3D& intpoints,
        const int                              iquad,
        const int                              eleid,
        const bool                             higher_order_ele);

    /*!
       \brief get material parameters
     */
    void GetMaterialParams(Teuchos::RCP<const MAT::Material>  material,
                           const LINALG::Matrix<3,iel>&       evelaf,
                           const LINALG::Matrix<iel,1>&       escaaf,
                           const LINALG::Matrix<iel,1>&       escaam,
                           const bool                         is_genalpha);

    /*!
      \brief calculate (all-scale) subgrid viscosity
    */
    void CalcSubgrVisc(Fluid3*                             ele,
                       const LINALG::Matrix<3,iel>&        evelaf,
                       LINALG::Matrix<4*iel,1>&            sgvelvisc,
                       const double                        vol,
                       const enum Fluid3::TurbModelAction  turb_mod_action,
                       double&                             Cs,
                       double&                             Cs_delta_sq,
                       double&                             l_tau,
                       const bool                          transfer);

    /*!
      \brief calculate fine-scale subgrid viscosity
    */
    void CalcFineScaleSubgrVisc(Fluid3*                             ele,
                                const LINALG::Matrix<3,iel>&        evelaf,
                                const LINALG::Matrix<3,iel>&        fsevelaf,
                                LINALG::Matrix<4*iel,1>&            sgvelvisc,
                                const double                        vol,
                                const enum Fluid3::FineSubgridVisc  fssgv,
                                double&                             Cs,
                                const bool                          transfer);

    /*!
      \brief calculate stabilization parameter

      Options for tau:

      o Franca-Barrenechea-Valentin-Wall:
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
    */
    void CalcStabParameter(const double               timefac,
                           const double               vol,
                           const enum Fluid3::TauType whichtau);

    /*!
      \brief calculate subgrid-scale velocity at element center
             for transfer to coupled scalar transport solver
    */
    void CalcSubgrVelocity(Fluid3*                      ele,
                           const LINALG::Matrix<3,iel>& evelaf,
                           const LINALG::Matrix<iel,1>& epreaf,
                           const LINALG::Matrix<3,iel>& eaccam,
                           const LINALG::Matrix<3,iel>& emhist,
                           LINALG::Matrix<4*iel,1>&     sgvelvisc,
                           const double                 timefac,
                           const bool                   loma,
                           const bool                   is_genalpha,
                           const bool                   higher_order_ele);

    /*!
      \brief calculate rate of strain of (fine-scale) velocity

      \param evel       (i) nodal velocity values
      \param derxy      (i) shape function derivatives
      \param velderxy   (o) velocity derivatives

      \return computed rate of strain
     */
    double GetStrainRate(const LINALG::Matrix<3,iel>& evel,
                         const LINALG::Matrix<3,iel>& derxy,
                         LINALG::Matrix<3,3>        & velderxy)
    {
      double rateofstrain=0;

      // get velocity derivatives at integration point
      //
      //              +-----  dN (x)
      //   dvel (x)    \        k
      //   -------- =   +     ------ * vel
      //      dx       /        dx        k
      //        j     +-----      j
      //              node k
      //
      // j : direction of derivative x/y/z
      //
      velderxy.MultiplyNT(evel,derxy);

      // compute (resolved) rate of strain
      //
      //          +-                                 -+ 1
      //          |          /   \           /   \    | -
      //          | 2 * eps | vel |   * eps | vel |   | 2
      //          |          \   / ij        \   / ij |
      //          +-                                 -+
      //
      LINALG::Matrix<3,3> two_epsilon;
      for(int rr=0;rr<3;++rr)
      {
        for(int mm=0;mm<3;++mm)
        {
          two_epsilon(rr,mm) = velderxy(rr,mm) + velderxy(mm,rr);
        }
      }

      for(int rr=0;rr<3;rr++)
      {
        for(int mm=0;mm<3;mm++)
        {
          rateofstrain += two_epsilon(rr,mm)*two_epsilon(mm,rr);
        }
      }

      // sqrt(two_epsilon(rr,mm)*two_epsilon(mm,rr)/4.0*2.0)

      return(sqrt(rateofstrain/2.0));
    }


  private:
    //! node coordinates
    LINALG::Matrix<3,iel> xyze_;
    //! dead load in nodes at n+alpha_F/n+1
    LINALG::Matrix<3,iel> edeadaf_;
    //! array for shape functions
    LINALG::Matrix<iel,1> funct_;
    //! array for shape function derivatives w.r.t r,s,t
    LINALG::Matrix<3,iel> deriv_;
    //! array for second derivatives of shape function w.r.t r,s,t
    LINALG::Matrix<6,iel> deriv2_;
    //! transposed jacobian "dx/ds"
    LINALG::Matrix<3,3> xjm_;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<3,3> xji_;
    //! global velocity derivatives in gausspoint w.r.t x,y,z
    LINALG::Matrix<3,3> vderxy_;
    //! fine-scale global velocity derivatives in gausspoint w.r.t x,y,z
    LINALG::Matrix<3,3> fsvderxy_;
    //! global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<3,iel> derxy_;
    //! global second derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<6,iel> derxy2_;
    //! bodyforce in gausspoint
    LINALG::Matrix<3,1> bodyforce_;
    //! vector containing all values from previous timelevel n for momentum equation
    LINALG::Matrix<3,1> histmom_;
    //! direction of flow (normed velocity vector)
    LINALG::Matrix<3,1> velino_;
    //! density-weighted velocity vector in gausspoint
    LINALG::Matrix<3,1> velint_;
    //! fine-scale velocity vector in gausspoint
    LINALG::Matrix<3,1> fsvelint_;
    //! subgrid-scale velocity vector in gausspoint
    LINALG::Matrix<3,1> sgvelint_;
    //! density-weighted ale convective velocity in gausspoint, c=rho*u-rho*u_G
    LINALG::Matrix<3,1> convvelint_;
    //! acceleration vector in gausspoint
    LINALG::Matrix<3,1> accint_;
    //! pressure gradient in gausspoint
    LINALG::Matrix<3,1> gradp_;
    //! the stabilisation parameters
    LINALG::Matrix<3,1> tau_;
    //! viscous term including 2nd derivatives
    //! (This array once had three dimensions, now the first two are combined to one.)
    LINALG::Matrix<3*3,iel> viscs2_;
    //! linearisation of convection, convective part
    LINALG::Matrix<iel,1> conv_c_;
    //! linearisation of subgrid-scale convection, convective part
    LINALG::Matrix<iel,1> sgconv_c_;
    //! velocity divergence
    double vdiv_;
    //! total right hand side terms at int.-point for momentum equation
    LINALG::Matrix<3,1> rhsmom_;
    //! (u_old*nabla)u_old
    LINALG::Matrix<3,1> conv_old_;
    //! div epsilon(u_old)
    LINALG::Matrix<3,1> visc_old_;
    //! old residual of momentum equation
    LINALG::Matrix<3,1> momres_old_;
    //! old residual of continuity equation
    double conres_old_;
    //! 2nd derivatives of coord.-functions w.r.t r,s,t
    LINALG::Matrix<6,3> xder2_;

    LINALG::Matrix<3,3> vderiv_;

    //! flag for material evaluation at Gaussian integration points
    bool mat_gp_;
    //! flag for stabilization parameter evaluation at Gaussian integration points
    bool tau_gp_;
    //! Jacobian determinant
    double det_;
    //! integration factor
    double fac_;
    //! physical viscosity
    double visc_;
    //! (all-scale) subgrid viscosity
    double sgvisc_;
    //! effective viscosity = physical viscosity + (all-scale) subgrid viscosity
    double visceff_;
    //! fine-scale subgrid viscosity
    double fssgvisc_;
    //! time-step length
    double dt_;
    //! factor for right-hand side due to time-integration scheme
    double omtheta_;

    //! LOMA-specific variables:
    //! right-hand-side term at int.-point for continuity equation
    double rhscon_;
    //! density at t_(n+alpha_F) or t_(n+1)
    double densaf_;
    //! density at t_(n+alpha_M)
    double densam_;
    //! factor for scalar time derivative
    double scadtfac_;
    //! factor for convective scalar term at t_(n+alpha_F) or t_(n+1)
    double scaconvfacaf_;
    //! factor for convective scalar term at t_(n)
    double scaconvfacn_;
    //! addition to continuity equation due to thermodynamic pressure
    double thermpressadd_;
    //! global velocity derivatives in gausspoint w.r.t x,y,z at t_(n)
    LINALG::Matrix<3,3> vderxyn_;
    //! scalar gradient at t_(n+alpha_F) or t_(n+1)
    LINALG::Matrix<3,1> grad_scaaf_;
    //! scalar gradient at t_(n)
    LINALG::Matrix<3,1> grad_scan_;
    //! convective scalar term at t_(n+alpha_F) or t_(n+1)
    double conv_scaaf_;
    //! convective scalar term at t_(n)
    double conv_scan_;
    //! thermodynamic pressure at t_(n+alpha_F) or t_(n+1)
    double thermpressaf_;
    //! thermodynamic pressure at t_(n+alpha_M) or t_(n)
    double thermpressam_;
    //! time derivative of thermodynamic pressure at t_(n+alpha_M) or t_(n+1)
    double thermpressdtam_;
  };

}
}

#endif

#endif
#endif
