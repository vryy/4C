/*----------------------------------------------------------------------*/
/*!
\file fluid2_impl.H

\brief Internal implementation of Fluid2 element (one-step-theta/BDF2)

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

#ifndef FLUID2_IMPL_H
#define FLUID2_IMPL_H

#include "fluid2.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/linalg_fixedsizematrix.H"
//#include <blitz/array.h>

namespace DRT
{
namespace ELEMENTS
{

  /// Interface base class for Fluid2Impl
  /*!
    This class exists to provide a common interface for all template
    versions of Fluid2Impl. The only function
    this class actually defines is Impl, which returns a pointer to
    the appropriate version of Fluid2Impl.
   */
  class Fluid2ImplInterface
  {
  public:
    /// Empty constructor
    Fluid2ImplInterface() {}
    /// Empty destructor
    virtual ~Fluid2ImplInterface() {}
    /// Evaluate the element
    /*!
      This class does not provide a definition for this function, it
      must be defined in Fluid3Impl.
     */
    virtual int Evaluate(Fluid2*                   ele,
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
    static Fluid2ImplInterface* Impl(DRT::ELEMENTS::Fluid2* f2);

  };


  /// Internal Fluid2 implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Fluid2 element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

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
  template <DRT::Element::DiscretizationType distype>
  class Fluid2Impl: public Fluid2ImplInterface
  {
  public:
    /// Constructor
    explicit Fluid2Impl();

    //! number of nodes
    static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    /// Evaluate
    /*!
      The evaluate function for the general fluid case.
      Why doesn't someone who actually knows what this function does
      write a useful comment?
     */
    virtual int Evaluate(Fluid2*                   ele,
                         ParameterList&            params,
                         DRT::Discretization&      discretization,
                         vector<int>&              lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         RefCountPtr<MAT::Material> mat);

    /*!
      \brief calculate element matrix and rhs

      \param ele              (i) the element those matrix is calculated
      \param evelnp           (i) nodal velocities at n+1
      \param fsevelnp         (i) fine-scale nodal velocities at n+1
      \param eprenp           (i) nodal pressure
      \param edensnp          (i) nodal density
      \param eaccam           (i) nodal accelerations at n+alpha_M
      \param ededtam          (i) nodal density derivatives at n+alpha_M
      \param edensam          (i) nodal densities at n+alpha_M
      \param emhist           (i) time rhs for momentum equation
      \param echist           (i) time rhs for continuity equation
      \param edispnp          (i) nodal displacements (on moving mesh)
      \param egridv           (i) grid velocity (on moving mesh)
      \param estif            (o) element matrix to calculate
      \param emesh            (o) linearization wrt mesh motion
      \param eforce           (o) element rhs to calculate
      \param sgvelvisc        (o) subgrid-velocity/viscosity vector
      \param material         (i) fluid material
      \param time             (i) current simulation time
      \param dt               (i) timestep
      \param timefac          (i) time discretization factor
      \param eosfac           (i) equation-of-state factor
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
      \param visceff          (i) viscosity + subgrid viscosity
      */
    void Sysmat(Fluid2*                                  ele,
                const LINALG::Matrix<2,iel>&             evelnp,
                const LINALG::Matrix<2,iel>&             fsevelnp,
                const LINALG::Matrix<iel,1>&             eprenp,
                const LINALG::Matrix<iel,1>&             edensnp,
                const LINALG::Matrix<2,iel>&             eaccam,
                const LINALG::Matrix<iel,1>&             ededtam,
                const LINALG::Matrix<iel,1>&             edensam,
                const LINALG::Matrix<2,iel>&             emhist,
                const LINALG::Matrix<iel,1>&             echist,
                const LINALG::Matrix<2,iel>&             edispnp,
                const LINALG::Matrix<2,iel>&             egridv,
                LINALG::Matrix<3*iel,3*iel>&             estif,
                LINALG::Matrix<3*iel,3*iel>&             emesh,
                LINALG::Matrix<3*iel,    1>&             eforce,
                LINALG::Matrix<3*iel,    1>&             sgvelvisc,
                Teuchos::RCP<const MAT::Material>        material,
                double                                   time,
                double                                   dt,
                double                                   timefac,
                const double                             eosfac,
                const bool                               newton,
                const bool                               loma,
                const bool                               conservative,
                const bool                               sgvel,
                const bool                               is_genalpha,
                const bool                               higher_order_ele,
                const enum Fluid2::FineSubgridVisc       fssgv,
                const enum Fluid2::StabilisationAction   pspg,
                const enum Fluid2::StabilisationAction   supg,
                const enum Fluid2::StabilisationAction   vstab,
                const enum Fluid2::StabilisationAction   cstab,
                const enum Fluid2::StabilisationAction   cross,
                const enum Fluid2::StabilisationAction   reynolds,
                const enum Fluid2::TauType               whichtau,
                const enum Fluid2::TurbModelAction       turb_mod_action,
                double&                                  Cs,
                double&                                  visceff);


  private:

    /*!
      \brief calculate various values at element center: stabilization parameter,
             (non-linear) viscosity, subgrid viscosity, subgrid velocity

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

      \param ele              (i) the element those matrix is calculated
      \param evelnp           (i) nodal velocities at n+1
      \param fsevelnp         (i) fine-scale nodal velocities at n+1
      \param eprenp           (i) nodal pressure
      \param edensnp          (i) nodal density
      \param eaccam           (i) nodal accelerations at n+alpha_M
      \param edensam          (i) nodal densities at n+alpha_M
      \param emhist           (i) time rhs for momentum equation
      \param sgvelvisc        (o) subgrid-velocity/viscosity vector
      \param material         (i) fluid material
      \param dt               (i) timestep
      \param timefac          (i) time discretization factor
      \param loma             (i) boolean flag for potential low-Mach-number solver
      \param conservative     (i) boolean flag for conservative form
      \param sgvel            (i) boolean flag for subgrid-scale velocity
      \param is_genalpha      (i) boolean flag for generalized-alpha time integration
      \param higher_order_ele (i) keep or drop second derivatives
      \param fssgv            (i) flag for type of fine-scale subgrid viscosity
      \param whichtau         (i) selecting available tau definition
      \param turb_mod_action  (i) selecting turbulence model (none, Smagorisky,
                                  dynamic Smagorinsky, Smagorinsky with van Driest
                                  damping for channel flows)
      \param Cs               (i) Smagorinsky model parameter
      \param Cs_delta_sq      (i) Model parameter computed by dynamic Smagorinsky
                                  approach (Cs*h*h)
      \param visc             (i) viscosity
      \param visceff          (i) viscosity + subgrid viscosity
    */
    void Caltau(Fluid2*                                 ele,
                const LINALG::Matrix<2,iel>&            evelnp,
                const LINALG::Matrix<2,iel>&            fsevelnp,
                const LINALG::Matrix<iel,1>&            eprenp,
                const LINALG::Matrix<iel,1>&            edensnp,
                const LINALG::Matrix<2,iel>&            eaccam,
                const LINALG::Matrix<iel,1>&            edensam,
                const LINALG::Matrix<2,iel>&            emhist,
                LINALG::Matrix<3*iel,1>&                sgvelvisc,
                Teuchos::RCP<const MAT::Material>       material,
                const double                            dt,
                const double                            timefac,
                const double                            eosfac,
                const bool                              loma,
                const bool                              conservative,
                const bool                              sgvel,
                const bool                              is_genalpha,
                const bool                              higher_order_ele,
                const enum Fluid2::FineSubgridVisc      fssgv,
                const enum Fluid2::TauType              whichtau,
                const enum Fluid2::TurbModelAction      turb_mod_action,
                double&                                 Cs,
                double&                                 visc,
                double&                                 visceff);


    /*!
       \brief calculate viscosity at element center

       \param material     (i) material
       \param visc         (o) kinematic material viscosity
       \param rateofshear  (i) shear rate
       \param dens         (i) density
       \param eosfac       (i) equation-of-state factor
     */
    void CalVisc(Teuchos::RCP<const MAT::Material>  material,
                 double&                            visc,
                 const double &                     rateofshear,
                 const double &                     dens,
                 const double &                     eosfac);

    /*!
     * \brief calculate body force from the node condition
     */
    void BodyForce(Fluid2*      ele,
                   const double time);

    /*!
      \brief calculate rate of strain of (fine-scale) velocity

      \param evel       (i) nodal velocity values
      \param derxy      (i) shape function derivatives
      \param velderxy   (o) velocity derivatives

      \return computed rate of strain
     */
    double GetStrainRate(const LINALG::Matrix<2,iel>& evel,
                         const LINALG::Matrix<2,iel>& derxy,
                         LINALG::Matrix<2,2>        & velderxy)
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
      // j : direction of derivative x/y
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
      LINALG::Matrix<2,2> two_epsilon;
      for(int rr=0;rr<2;++rr)
      {
        for(int mm=0;mm<2;++mm)
        {
          two_epsilon(rr,mm) = velderxy(rr,mm) + velderxy(mm,rr);
        }
      }

      for(int rr=0;rr<2;rr++)
      {
        for(int mm=0;mm<2;mm++)
        {
          rateofstrain += two_epsilon(rr,mm)*two_epsilon(mm,rr);
        }
      }

      // sqrt(two_epsilon(rr,mm)*two_epsilon(mm,rr)/4.0*2.0)

      return(sqrt(rateofstrain/2.0));
    }

  private:
    //! artificial viscosity
    double vart_;
    //! node coordinates
    LINALG::Matrix<2,iel> xyze_;
    //! dead load in nodes
    LINALG::Matrix<2,iel> edeadng_;
    //! array for shape functions
    LINALG::Matrix<iel,1> funct_;
    //! array for density-weighted shape functions
    LINALG::Matrix<iel,1> densfunct_;
    //! array for density-weighted shape functions at n+alphaM
    LINALG::Matrix<iel,1> densamfunct_;
    //! array for shape functions divided by density
    LINALG::Matrix<iel,1> functdens_;
    //! array for shape function derivatives w.r.t r,s,t
    LINALG::Matrix<2,iel> deriv_;
    //! array for second derivatives of shape function w.r.t r,s,t
    LINALG::Matrix<3,iel> deriv2_;
    //! transposed jacobian "dx/ds"
    LINALG::Matrix<2,2> xjm_;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<2,2> xji_;
    //! global velocity derivatives in gausspoint w.r.t x,y,z
    LINALG::Matrix<2,2> vderxy_;
    //! global momentum derivatives in gausspoint w.r.t x,y,z
    LINALG::Matrix<2,2> mderxy_;
    //! fine-scale global velocity derivatives in gausspoint w.r.t x,y,z
    LINALG::Matrix<2,2> fsvderxy_;
    //! global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<2,iel> derxy_;
    //! density-weighted global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<2,iel> densderxy_;
    //! global second derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<3,iel> derxy2_;
    //! bodyforce in gausspoint
    LINALG::Matrix<2,1> bodyforce_;
    //! vector containing all values from previous timelevel n for momentum equation
    LINALG::Matrix<2,1> histmom_;
    //! value from previous timelevel n for continuity equation
    double histcon_;
    //! direction of flow (normed velocity vector, currently not used)
    //blitz::Array<double,1> velino_;
    //! density-weighted velocity vector in gausspoint
    LINALG::Matrix<2,1> velint_;
    //! non-density-weighted velocity vector in gausspoint
    LINALG::Matrix<2,1> ndwvelint_;
    //! fine-scale velocity vector in gausspoint
    LINALG::Matrix<2,1> fsvelint_;
    //! subgrid-scale velocity vector in gausspoint
    LINALG::Matrix<2,1> sgvelint_;
    //! density-weighted ale convective velocity in gausspoint, c=rho*u-rho*u_G
    LINALG::Matrix<2,1> convvelint_;
    //! acceleration vector at n+alphaM in gausspoint
    LINALG::Matrix<2,1> accintam_;
    //! pressure gradient in gausspoint
    LINALG::Matrix<2,1> gradp_;
    //! the stabilisation parameters
    LINALG::Matrix<3,1> tau_;
    //! viscous term including 2nd derivatives
    LINALG::Matrix<4,iel> viscs2_;
    //! linearisation of convection, convective part
    LINALG::Matrix<iel,1> conv_c_;
    //! linearisation of subgrid-scale convection, convective part
    LINALG::Matrix<iel,1> sgconv_c_;
    //! momentum divergence
    double mdiv_;
    //! velocity divergence
    double vdiv_;
    //! total right hand side terms at int.-point for momentum equation
    LINALG::Matrix<2,1> rhsmom_;
    //! total right hand side term at int.-point for continuity equation
    double rhscon_;
    //! (u_old*nabla)u_old
    LINALG::Matrix<2,1> conv_old_;
    //! div epsilon(u_old)
    LINALG::Matrix<2,1> visc_old_;
    //! the old residual
    LINALG::Matrix<2,1> res_old_;
    //! 2nd derivatives of coord.-functions w.r.t r,s,t
    LINALG::Matrix<3,2> xder2_;

    LINALG::Matrix<2,2> vderiv_;
  };

}
}

#endif

#endif
#endif
