/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjoint3_impl.H

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef TOPOPT_FLUIDADJOINT3_IMPL_H_
#define TOPOPT_FLUIDADJOINT3_IMPL_H_


#include "topopt_fluidAdjoint3_interface.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"


namespace FLD
{
template<DRT::Element::DiscretizationType distype> class RotationallySymmetricPeriodicBC;
}

namespace DRT
{
namespace ELEMENTS
{
class FluidAdjoint3ImplParameter;


/// Internal FluidAdjoint3 implementation
/*!
    This internal class keeps all the working arrays needed to calculate
    the Adjoint of the Fluid3 element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Fluid3) from the
    mathematical contents (this class). There are different implementations of
    the Fluid3 element, this is just one such implementation.

    The Fluid3 will allocate exactly one object of this class for all
    fluid elements with the same number of nodes in the mesh. This
    allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    <h3>Usability</h3>

    The calculations are done by the Evaluate() method. The non-virtual
    Evaluate() method must be callable without an actual Fluid3 object.

    <h3>History</h3>

    The implementation here is based on the work by Christiane
    Foerster. The core element code has been linearized using MuPAD.
    The reaction terms have been removed from the lhs.

    \author winklmaier
    \date 02/12
 */
template<DRT::Element::DiscretizationType distype>
class FluidAdjoint3Impl : public FluidAdjoint3ImplInterface
{
  friend class FluidAdjoint3ImplParameter;

  /// private Constructor since we are a Singleton.
  FluidAdjoint3Impl();

public:

  virtual ~FluidAdjoint3Impl() {}

  /// Singleton access method
  static FluidAdjoint3Impl<distype> * Instance( bool create=true );

  //! nen_: number of element nodes (P. Hughes: The Finite Element Method)
  static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions
  static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  static const int numdofpernode_ = nsd_+1;

  /// called upon destruction
  virtual void Done();

  /// Evaluate element ERROR
  /*!
        general function to compute the error (analytical solution) for particular problem type
   */
  virtual int ComputeError(
      DRT::ELEMENTS::Fluid3*          ele,
      ParameterList&                  params,
      Teuchos::RCP<MAT::Material>&    mat,
      DRT::Discretization&            discretization,
      vector<int>&                    lm,
      Epetra_SerialDenseVector&       elevec1);

  /// Evaluate the element
  /*!
      Generic virtual interface function. Called via base pointer.
   */
  virtual int Evaluate(
      DRT::ELEMENTS::Fluid3*        ele,
      DRT::Discretization &         discretization,
      const std::vector<int> &      lm,
      Teuchos::ParameterList&       params,
      Teuchos::RCP<MAT::Material> & mat,
      Epetra_SerialDenseMatrix&     elesysmat,
      Epetra_SerialDenseMatrix&     eledummymat,
      Epetra_SerialDenseVector&     elerhs,
      Epetra_SerialDenseVector&     eledummyrhs1,
      Epetra_SerialDenseVector&     eledummyrhs2 );

  /// Evaluate the element at specified gauss points
  virtual int Evaluate(
      DRT::ELEMENTS::Fluid3*               ele,
      DRT::Discretization &                discretization,
      const std::vector<int> &             lm,
      Teuchos::ParameterList&              params,
      Teuchos::RCP<MAT::Material> &        mat,
      Epetra_SerialDenseMatrix&            elesysmat,
      Epetra_SerialDenseVector&            elerhs,
      const DRT::UTILS::GaussIntegration & intpoints );

  /*!
      \brief evaluate function for Fluid3 element

      Specific evaluate function without any knowledge about DRT objects. This
      way the element evaluation is independent of the specific mesh storage.
   */
  int Evaluate(
      int                                           eid,
      Teuchos::ParameterList&                       params,
      LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> & elesysmat,
      LINALG::Matrix<(nsd_+1)*nen_,            1> & elerhs,
      const LINALG::Matrix<nsd_,nen_> &             eveln,
      const LINALG::Matrix<nsd_,nen_> &             evelnp,
      const LINALG::Matrix<nen_,1>    &             epren,
      const LINALG::Matrix<nen_,1>    &             eprep,
      const LINALG::Matrix<nsd_,nen_> &             efluidveln,
      const LINALG::Matrix<nsd_,nen_> &             efluidvelnp,
      const LINALG::Matrix<nen_,1> &                eporo,
      Teuchos::RCP<MAT::Material>                   mat,
      const DRT::UTILS::GaussIntegration &          intpoints);

  /*!
      \brief calculate element matrix and rhs

      \param eid              (i) element id
      \param ebofoaf          (i) body force at n+alpha_F/n+1
      \param evelaf           (i) nodal velocities at n+alpha_F/n+1
      \param evelnp           (i) nodal velocities at n+1 (np_genalpha)
      \param epreaf           (i) nodal pressure at n+alpha_F/n+1
      \param eaccam           (i) nodal accelerations at n+alpha_M
      \param emhist           (i) time rhs for momentum equation
      \param estif            (o) element matrix to calculate
      \param eporo            (i) nodal porosities
      \param eforce           (o) element rhs to calculate
      \param material         (i) fluid material
      \param intpoints        (i) Gaussian integration points

   */
  void Sysmat(
      int                                            eid,
      const LINALG::Matrix<nsd_,nen_>&               eveln,
      const LINALG::Matrix<nsd_,nen_>&               evelnp,
      const LINALG::Matrix<nen_,1>&                  epren,
      const LINALG::Matrix<nen_,1>&                  eprenp,
      const LINALG::Matrix<nsd_,nen_> &              efluidveln,
      const LINALG::Matrix<nsd_,nen_> &              efluidvelnp,
      LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&   estif,
      LINALG::Matrix<(nsd_+1)*nen_,    1>&           eforce,
      const LINALG::Matrix<nen_,1> &                 eporo,
      Teuchos::RCP<const MAT::Material>              material,
      const DRT::UTILS::GaussIntegration &           intpoints);

  LINALG::Matrix<nsd_,nen_> & xyz() { return xyze_; }

private:

  //! number of components necessary to store second derivatives
  /*!
     1 component  for nsd=1:  (N,xx)

     3 components for nsd=2:  (N,xx ; N,yy ; N,xy)

     6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
   */
  static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

  //! evaluate shape functions and their derivatives at element center
  void EvalShapeFuncAndDerivsAtEleCenter(
      const int  eleid  ///< element ID
  );

  //! brief evaluate shape functions and their derivatives at integration point
  void EvalShapeFuncAndDerivsAtIntPoint(
      DRT::UTILS::GaussIntegration::iterator & iquad,             ///< actual integration point
      const int                                eleid              ///< element ID
  );


  //! get material parameters
  void GetMaterialParams(Teuchos::RCP<const MAT::Material>       material,  ///< reference pointer to material
      const LINALG::Matrix<nsd_,nen_>&        evelaf    ///< velocity at time n+alpha_f / n+1
  );

  //! calculate stabilization parameter
  void CalcStabParameter(const double  vol);  ///< volume

  //! calculate characteristic element length
  void CalcCharEleLength(const double  vol,       ///< volume
      const double  fluidvel_norm,  ///< norm of velocity vector
      double&       strle,     ///< streamlength (for tau_Mu)
      double&       hk);       ///< length for tau_Mp

  void MassReactionGalPart(
      LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_u,       ///< block (test-vel x shp-vel)
      LINALG::Matrix<nsd_,nen_> &             velforce,      ///< rhs forces velocity
      const double &                          timefacfac,    ///< = timefac x fac
      const double &                          timefacfacrhs  ///< = timefacrhs * fac
  );

  //! Compute element matrix and rhs entries: convective terms of the Galerkin part
  void ConvectionGalPart(
      LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_u,       ///< block (test-vel x shp-vel)
      LINALG::Matrix<nsd_,nen_> &             velforce,      ///< rhs forces velocity
      const double &                          timefacfac,    ///< = timefac x fac
      const double &                          timefacfacrhs  ///< = timefacrhs * fac
  );

  //! Compute element matrix entries: for the viscous terms of the Galerkin part
  void ViscousGalPart(
      LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_u,        ///< block (test-vel x shp-vel)
      LINALG::Matrix<nsd_,nen_> &             velforce,       ///< rhs forces velocity
      const double &                          timefacfac,     ///< = timefac x fac
      const double &                          timefacfacrhs   ///< = timefacrhs * fac
  );

  //! Compute element matrix entries: pressure terms of the Garlerkin part and rhs
  void PressureGalPart(
      LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,        ///< block (test-vel x shp-vel)
      LINALG::Matrix<nsd_,nen_> &               velforce,         ///< rhs forces velocity
      const double &                            timefacfacpre,    ///< = timefacpre x fac
      const double &                            timefacfacrhspre, ///< = timefacrhspre x fac
      const double &                            press,            ///< pressure at integration point
      const double &                            press_old         ///< old pressure at integration point
  );

  //! Compute element matrix entries: continuity terms of the Garlerkin part and rhs
  void ContinuityGalPart(
      LINALG::Matrix<nen_, nen_*nsd_> &         estif_q_u,        ///< block (weighting function q x u)
      LINALG::Matrix<nen_,1> &                  preforce,         ///< rhs forces pressure
      const double &                            timefacfacdiv,    ///< = timefac x fac
      const double &                            timefacfacdivrhs  ///< right-hand-side factor
  );

  //! Compute element matrix entries: body force terms on rhs
  void BodyForceGalPart(
      LINALG::Matrix<nsd_,nen_> &               velforce,
      const double &                            timefacfac,
      const double &                            timefacfacrhs
  );

  void MomRes(
      LINALG::Matrix<nsd_*nsd_,nen_> &        lin_resM_Du,    ///< linearisation of the Garlerkin momentum residual
      LINALG::Matrix<nsd_,1> &                  StrongResMomScaled,
      const double &                          timefacfac,     ///< = timefac x fac
      const double &                          timefacfacrhs,
      const double &                          timefacfacpre,
      const double &                          timefacfacprerhs,
      const LINALG::Matrix<nsd_,nen_>&        eveln,
      const LINALG::Matrix<nsd_,nen_>&        evelnp,
      const LINALG::Matrix<nsd_,nen_>&        efluidveln,
      const LINALG::Matrix<nsd_,nen_>&        efluidvelnp
  );

  //! Compute element matrix and rhs entries: inertia and reactive terms of the Galerkin part
  void CalcDivEps(
      const LINALG::Matrix<nsd_,nen_>&      eveln,
      const LINALG::Matrix<nsd_,nen_>&      evelnp,
      LINALG::Matrix<nsd_,1>&               viscs,
      LINALG::Matrix<nsd_,1>&               viscs_old,
      LINALG::Matrix<nsd_*nsd_,nen_>&       visc_shp
  );

  /// calculate body force from nodal conditions
  void BodyForce(
      const LINALG::Matrix<nsd_,nen_>&           efluidveln,
      const LINALG::Matrix<nsd_,nen_>&           efluidvelnp,
      LINALG::Matrix<nsd_,1>&                    bodyforce,
      LINALG::Matrix<nsd_,1>&                    bodyforce_old
  );

  //! Compute element matrix entries: PSPG
  void PSPG(
      LINALG::Matrix<nen_, nen_*nsd_> &         estif_q_u,      ///< block (weighting function q x u)
      LINALG::Matrix<nen_,nen_> &               ppmat,          ///< block (weighting function q x p)
      LINALG::Matrix<nen_,1> &                  preforce,       ///< rhs forces pressure
      LINALG::Matrix<nsd_*nsd_,nen_> &          GalMomResnU,    ///< linearisation of the stabilization residual
      LINALG::Matrix<nsd_,1> &                  StrongResMomScaled,
      const double &                            timefacfac,     ///< = timefac x fac
      const double &                            timefacfacrhs,  ///< = timefacrhs x fac
      const double &                            timefacfacpre,  ///< = timefacpre x fac
      const double &                            timefacfacprerhs///< = timefacprerhs x fac
  );

  //! Compute element matrix entries: SUPG
  void SUPG(
      LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,        ///< block (weighting function v x u)
      LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,      ///< block (weighting function v x p)
      LINALG::Matrix<nsd_,nen_> &               velforce,       ///< rhs forces velocity
      LINALG::Matrix<nsd_*nsd_,nen_> &          GalMomResnU,    ///< linearisation of the stabilization residual
      LINALG::Matrix<nsd_,1> &                  StrongResMomScaled,
      const double &                            timefacfac,     ///< = timefac x fac
      const double &                            timefacfacrhs,  ///< = timefacrhs x fac
      const double &                            timefacfacpre,  ///< = timefacpre x fac
      const double &                            timefacfacprerhs///< = timefacprerhs x fac
  );

  void ContRes(
      double &                                  StrongResContScaled,
      const double &                            timefacfacdiv,  ///< = timefac x fac
      const double &                            timefacfacdivrhs///< = timefac x fac
  );
  void ContForce(
      double &                            contForce,  ///< = timefac x fac
      double &                            contForce_old///< = timefac x fac
  );

      //! Compute element matrix entries: div-grad stabilization and the rhs of the viscous term
  void ContStab(
      LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,        ///< block (weighting function v x u)
      LINALG::Matrix<nsd_,nen_> &               velforce,       ///< rhs forces velocity
      const double &                            timefacfacdiv,  ///< = timefac x fac
      const double &                            timefacfacdivrhs///< = timefac x fac
  );

  /*!
   * \brief fill elment matrix and vectors with the global values
   */
  void ExtractValuesFromGlobalVector( const DRT::Discretization&   discretization, ///< discretization
      const vector<int>&           lm,             ///<
      LINALG::Matrix<nsd_,nen_> *  matrixtofill,   ///< vector field
      LINALG::Matrix<nen_,1> *     vectortofill,   ///< scalar field
      const std::string            state          ///< state of the global vector
  );



  //! @name coordinate data

  //! node coordinates
  LINALG::Matrix<nsd_,nen_> xyze_;
  //@}


  //! @name shape functions

  //! array for shape functions
  LINALG::Matrix<nen_,1> funct_;
  //! array for shape function derivatives w.r.t r,s,t
  LINALG::Matrix<nsd_,nen_> deriv_;
  //! array for second derivatives of shape function w.r.t r,s,t
  LINALG::Matrix<numderiv2_,nen_> deriv2_;
  //! transposed jacobian "dx/ds"
  LINALG::Matrix<nsd_,nsd_> xjm_;
  //! inverse of transposed jacobian "ds/dx"
  LINALG::Matrix<nsd_,nsd_> xji_;
  //! global derivative of shape functions w.r.t x,y,z
  LINALG::Matrix<nsd_,nen_> derxy_;
  //! global second derivatives of shape functions w.r.t x,y,z
  LINALG::Matrix<numderiv2_,nen_> derxy2_;
  //@}


  //! @name variable values at last iteration (at time t^n)

  //! velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> velint_;
  //! global velocity derivatives in gausspoint w.r.t x,y,z
  LINALG::Matrix<nsd_,nsd_> vderxy_;
  //! pressure gradient in gausspoint
  LINALG::Matrix<nsd_,1> gradp_;
  //! fluid velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> fluidvelint_;
  //! global fluid derivative of velocity w.r.t x,y,z
  LINALG::Matrix<nsd_,nsd_> fluidvelxy_;
  //! velocity divergence
  double vdiv_;
  //! u_j v_i,j
  LINALG::Matrix<nsd_,1> conv1_;
  //! u_j,i * v_j
  LINALG::Matrix<nsd_,1> conv2_;
  //@}

  //! @name variable values at last time step (= t^n+1)

  //! velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> velint_old_;
  //! global velocity derivatives in gausspoint w.r.t x,y,z
  LINALG::Matrix<nsd_,nsd_> vderxy_old_;
  //! pressure gradient in gausspoint
  LINALG::Matrix<nsd_,1> gradp_old_;
  //! fluid velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> fluidvelint_old_;
  //! global fluid derivative of velocity w.r.t x,y,z
  LINALG::Matrix<nsd_,nsd_> fluidvelxy_old_;
  //! velocity divergence
  double vdiv_old_;
  //! u_j v_i,j
  LINALG::Matrix<nsd_,1> conv1_old_;
  //! u_j,i * v_j
  LINALG::Matrix<nsd_,1> conv2_old_;
  //@}


  //! the stabilisation parameters -> it is a (3,1) vector for 2D and 3D
  LINALG::Matrix<3,1> tau_;


  //! @name gauss point data

  // gauss points
  DRT::UTILS::GaussIntegration intpoints_;
  //! coordinates of current integration point in reference coordinates
  LINALG::Matrix<nsd_,1> xsi_;
  //! Jacobian determinant
  double det_;
  //! integration factor
  double fac_;
  //@}


  //! @name physical coefficients

  //! physical viscosity
  double visc_;
  //! reaction coefficient
  double reacoeff_;
  //! density
  double dens_;
  //@}

  //! Flag to (de)activate higher order elements
  //! elements with only mixed second order derivatives are not counted as higher order elements
  //! (see definition of higher order elements in fluid3_ele_impl_utils.cpp)
  bool is_higher_order_ele_;

  //! pointer to parameter list
  DRT::ELEMENTS::FluidAdjoint3ImplParameter* fluidAdjoint3Parameter_;
};
}
}

#endif /* TOPOPT_FLUIDADJOINT3_IMPL_H_ */
#endif  // #ifdef CCADISCRET
