
/*----------------------------------------------------------------------*/
/*!
\file artery_impl.H

\brief Internal implementation of Artery element

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_ARTNET
#ifdef CCADISCRET

#ifndef ARTERY_IMPL_H
#define ARTERY_IMPL_H

#include "artery.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/linalg_fixedsizematrix.H"
//#include <blitz/array.h>

namespace DRT
{
namespace ELEMENTS
{

  /// Interface base class for ArteryImpl
  /*!
    This class exists to provide a common interface for all template
    versions of AretryImpl. The only function
    this class actually defines is Impl, which returns a pointer to
    the appropriate version of ArteryImpl.
   */
  class AretryImplInterface
  {
  public:
    /// Empty constructor
    ArteryImplInterface() {}
    /// Empty destructor
    virtual ~ArteryImplInterface() {}
    /// Evaluate the element
    /*!
      This class does not provide a definition for this function, it
      must be defined in ArteryImpl.
     */
    virtual int Evaluate(Artery*                   ele,
                         ParameterList&            params,
                         DRT::Discretization&      discretization,
                         vector<int>&              lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         RefCountPtr<MAT::Material> mat,
                         _MATERIAL* actmat) = 0;

    /// Internal implementation class for artery element
    static ArteryImplInterface* Impl(DRT::ELEMENTS::Artery* art);

  };


  /// Internal artery implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the artery element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class artery)
    from the mathematical contents (this class). Of course there are
    different implementations of the artery element, this is just one
    such implementation.

    The Artery will allocate exactly one object of this class for all
    artery elements with the same number of nodes in the mesh. This
    allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    For matrix calculation the blitz library is used.

    \author ism
    \date 01/09
  */
  template <DRT::Element::DiscretizationType distype>
  class ArteryImpl: public ArteryImplInterface
  {
  public:
    /// Constructor
    explicit ArteryImpl();

    //! number of nodes
    static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    /// Evaluate
    /*!
      The evaluate function for the general artery case.
      Why doesn't someone who actually knows what this function does
      write a useful comment?
     */
    virtual int Evaluate(Artery*                   ele,
                         ParameterList&            params,
                         DRT::Discretization&      discretization,
                         vector<int>&              lm,
                         Epetra_SerialDenseMatrix& elemat1_epetra,
                         Epetra_SerialDenseMatrix& elemat2_epetra,
                         Epetra_SerialDenseVector& elevec1_epetra,
                         Epetra_SerialDenseVector& elevec2_epetra,
                         Epetra_SerialDenseVector& elevec3_epetra,
                         RefCountPtr<MAT::Material> mat,
                         _MATERIAL* actmat);
    /// calculating element matrix and rhs
    /*!
      Calculate matrix and rhs. Here the whole thing is hidden.

      \note Stationary calculations are not (yet) supported.

      \param ele              (i) the element whose matrix is calculated
      \param evelnp           (i) nodal velocities at n+1
      \param eveln            (i) nodal velocities at n
      \param eareanp          (i) nodal area at n+1
      \param earean           (i) nodal area at n
      \param emhist           (i) time rhs for momentum equation
      \param echist           (i) time rhs for continuity equation
      \param estif            (o) element stiffness matrix to calculate
      \param emass            (o) element mass matrix to calculate
      \param eforce           (o) element rhs to calculate
      \param material         (i) artery material
      \param time             (i) current simulation time
      \param dt               (i) timestep
      \param lin              (i) boolean flag for linearisation
      */
    void Sysmat(Artery*                               ele,
                const LINALG::Matrix<2,1>&            evelnp,
                const LINALG::Matrix<2,1>&            eveln,
                const LINALG::Matrix<2,1>&            eareanp,
                const LINALG::Matrix<2,1>&            earean,
                LINALG::Matrix<2,1>&                  estif,
                LINALG::Matrix<2,2>&                  emesh,
                struct _MATERIAL*                     material,
                double                                time,
                double                                dt,
                bool                                  lin);


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
      \param visceff          (i) nu+nu_turbulent
      \param fssgv            (i) flag for type of fine-scale subgrid viscosity

    */
    void Caltau(Fluid2*                                 ele,
                const LINALG::Matrix<2,iel>&           evelnp,
                const LINALG::Matrix<2,iel>&           fsevelnp,
                const LINALG::Matrix<iel,1>&           edensnp,
                const enum Fluid2::TauType              whichtau,
                struct _MATERIAL*                       material,
                double&                                 visc,
                const double                            timefac,
                const double                            dt,
                const enum Fluid2::TurbModelAction      turb_mod_action,
                double&                                 Cs,
                double&                                 visceff,
                const enum Fluid2::StabilisationAction  fssgv);


    /*!
       \brief calculates viscosity at a given point (element center or Gaussian point)

       \param material  (i) material
       \param visc      (o) kinematic material viscosity
     */
    void CalVisc(const struct _MATERIAL*   material,
                 double&                   visc);

    /*!
     * \brief calculate the body force from the nodes condition
     */
    void BodyForce(Fluid2*      ele,
                   const double time);

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
    //! velocity vector in gausspoint
    LINALG::Matrix<2,1> velint_;
    //! fine-scale velocity vector in gausspoint
    LINALG::Matrix<2,1> fsvelint_;
    //! ale convective velocity in gausspoint, c=u-u_G
    LINALG::Matrix<2,1> convvelint_;
    //! pressure gradient in gausspoint
    LINALG::Matrix<2,1> gradp_;
    //! the stabilisation parameters
    LINALG::Matrix<3,1> tau_;
    //! viscous term including 2nd derivatives
    LINALG::Matrix<4,iel> viscs2_;
    //! linearisation of convect, convective part
    LINALG::Matrix<iel,1> conv_c_;
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
    //! convection by subgrid-scale velocity
    LINALG::Matrix<iel,1> conv_resM_;
    //! 2nd derivatives of coord.-functions w.r.t r,s,t
    LINALG::Matrix<3,2> xder2_;

    LINALG::Matrix<2,2> vderiv_;
  };

}
}

#endif

#endif
#endif
