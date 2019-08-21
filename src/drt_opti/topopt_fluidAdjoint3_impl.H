/*---------------------------------------------------------------------*/
/*! \file

\brief Interface for element level of adjoint fluid equations

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef TOPOPT_FLUIDADJOINT3_IMPL_H_
#define TOPOPT_FLUIDADJOINT3_IMPL_H_


#include "topopt_fluidAdjoint3_interface.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"


namespace FLD
{
  template <DRT::Element::DiscretizationType distype>
  class RotationallySymmetricPeriodicBC;
}

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidAdjoint3ImplParameter;


    /// Internal FluidAdjoint3 implementation
    /*!
        This internal class keeps all the working arrays needed to calculate
        the Adjoint of the Fluid element. Additionally the method Sysmat()
        provides a clean and fast element implementation.

        <h3>Purpose</h3>

        The idea is to separate the element maintenance (class Fluid) from the
        mathematical contents (this class). There are different implementations of
        the Fluid element, this is just one such implementation.

        The FluidAdjoint3 will allocate exactly one object of this class for all
        fluid elements with the same number of nodes in the mesh. This
        allows us to use exactly matching working arrays (and keep them
        around.)

        The instationary adjoint equations start at the endtime T and stop at time 0
        Nevertheless in order to remain consistent, the variables at the new time
        step end with "...np", the old ones with "...n"

        <h3>Names</h3>
        The adjoint equations contain one more variable compared to the primal equations:
        - adjoint velocity v
        - primal velocity u (out of the Navier-Stokes equations)
        - adjoint pressure q
        - primal pressure p (out of the Navier-Stokes equations, only used for
          objective evaluations on the right-hand side)
        - velocity test function w
        - pressure test function r

        \author winklmaier
        \date 02/12
     */
    template <DRT::Element::DiscretizationType distype>
    class FluidAdjoint3Impl : public FluidAdjoint3ImplInterface
    {
      friend class FluidAdjoint3ImplParameter;

      /// private Constructor since we are a Singleton.
      FluidAdjoint3Impl();

     public:
      virtual ~FluidAdjoint3Impl() {}

      /// Singleton access method
      static FluidAdjoint3Impl<distype>* Instance(bool create = true);

      //! nen_: number of element nodes (P. Hughes: The Finite Element Method)
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      /// number of degrees of freedom per node
      static const int numdofpernode_ = nsd_ + 1;

      /// called upon destruction
      virtual void Done();

      /// Evaluate the element
      /*!
          Generic virtual interface function. Called via base pointer.
       */
      virtual int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elesysmat,
          Epetra_SerialDenseVector& elerhs);

      /// Evaluate the element at specified gauss points
      virtual int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elesysmat,
          Epetra_SerialDenseVector& elerhs, const DRT::UTILS::GaussIntegration& intpoints);

      /*!
          \brief calculate element matrix and rhs

          Specific evaluate routine without any knowledge about DRT objects. This
          way the element evaluation is independent of the specific mesh storage.

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
      void Sysmat(const LINALG::Matrix<nsd_, nen_>& eveln, const LINALG::Matrix<nsd_, nen_>& evelnp,
          const LINALG::Matrix<nen_, 1>& epren, const LINALG::Matrix<nen_, 1>& eprenp,
          const LINALG::Matrix<nsd_, nen_>& efluidveln,
          const LINALG::Matrix<nsd_, nen_>& efluidvelnp,
          const LINALG::Matrix<nsd_, nen_>& efluidvelnpp, const LINALG::Matrix<nen_, 1>& efluidpren,
          const LINALG::Matrix<nen_, 1>& efluidprenp, const LINALG::Matrix<nen_, 1>& efluidprenpp,
          const LINALG::Matrix<nsd_, nen_>& efluidbodyforcenp,
          const LINALG::Matrix<nsd_, nen_>& efluidbodyforcenpp,
          LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
          LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce, const LINALG::Matrix<nen_, 1>& edens,
          Teuchos::RCP<const MAT::Material> material,
          const DRT::UTILS::GaussIntegration& intpoints);

     private:
      //! number of components necessary to store second derivatives
      /*!
         1 component  for nsd=1:  (N,xx)

         3 components for nsd=2:  (N,xx ; N,yy ; N,xy)

         6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
       */
      static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

      //! evaluate shape functions and their derivatives at element center
      void EvalShapeFuncAndDerivsAtEleCenter();

      //! brief evaluate shape functions and their derivatives at integration point
      void EvalShapeFuncAndDerivsAtIntPoint(
          DRT::UTILS::GaussIntegration::iterator& iquad  ///< actual integration point
      );

      //! calculate stabilization parameter
      void CalcStabParameter(const double vol,      ///< volume
          const LINALG::Matrix<nsd_, 1>& fluidvel,  ///< velocity vector
          LINALG::Matrix<3, 1>& tau                 ///< stabilization parameter
      );

      //! calculate characteristic element length
      void CalcCharEleLength(const double vol,     ///< volume
          const LINALG::Matrix<nsd_, 1> fluidvel,  ///< velocity vector
          const double fluidvel_norm,              ///< norm of velocity vector
          double& strle,                           ///< streamlength (for tau_Mu)
          double& hk                               ///< length for tau_Mp
          ) const;

      /// calculate body force from nodal conditions
      void BodyForce(const LINALG::Matrix<nsd_, nen_>& efluidveln,
          const LINALG::Matrix<nsd_, nen_>& efluidvelnp);

      /// calculate body force from nodal conditions
      void FluidBodyForce(const double time, LINALG::Matrix<nsd_, nen_>& efluidbodyforce);

      /// compute force term of continuity equation
      void ContForce();

      /// compute element matrix and rhs entries: mass matrix and reaction terms of the Galerkin
      /// part
      void MassReactionGalPart(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (test-vel x shp-vel)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          const double& timefacfac,                           ///< = timefac x fac
          const double& timefacfacrhs                         ///< = timefacrhs * fac
          ) const;

      //! Compute element matrix and rhs entries: convective terms of the Galerkin part
      void ConvectionGalPart(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (test-vel x shp-vel)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          const double& timefacfac,                           ///< = timefac x fac
          const double& timefacfacrhs                         ///< = timefacrhs * fac
          ) const;

      //! Compute element matrix entries: for the viscous terms of the Galerkin part
      void ViscousGalPart(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (test-vel x shp-vel)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          const double& timefacfac,                           ///< = timefac x fac
          const double& timefacfacrhs                         ///< = timefacrhs * fac
          ) const;

      //! Compute element matrix entries: pressure terms of the Garlerkin part and rhs
      void PressureGalPart(
          LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,  ///< block (test-vel x shp-vel)
          LINALG::Matrix<nsd_, nen_>& velforce,          ///< rhs forces velocity
          const double& timefacfacpre,                   ///< = timefacpre x fac
          const double& timefacfacrhspre                 ///< = timefacrhspre x fac
          ) const;

      //! Compute element matrix entries: continuity terms of the Garlerkin part and rhs
      void ContinuityGalPart(
          LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          LINALG::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          const double& timefacfacdiv,                   ///< = timefac x fac
          const double& timefacfacdivrhs                 ///< right-hand-side factor
          ) const;

      //! Compute element matrix entries: body force terms on rhs
      void BodyForceGalPart(LINALG::Matrix<nsd_, nen_>& velforce, const double& timefacfac,
          const double& timefacfacrhs) const;

      /// compute element matrix entries: continuous force terms on rhs
      void ContForceGalPart(LINALG::Matrix<nen_, 1>& preforce, const double& timefacfacdiv,
          const double& timefacfacdivrhs) const;

      /// compute momentum residuum
      void DiscreteGalMom(LINALG::Matrix<nsd_ * nsd_, nen_>& GalMomTestStat,
          const double& timefacfac, const double& timefacfacrhs, const double& timefacfacpre,
          const double& timefacfacprerhs) const;

      /// compute momentum residuum
      void MomRes(LINALG::Matrix<nsd_ * nsd_, nen_>&
                      lin_resM_Du,  ///< linearisation of the Garlerkin momentum residual
          LINALG::Matrix<nsd_, 1>& StrongResMomScaled,
          const double& timefacfac,  ///< = timefac x fac
          const double& timefacfacrhs, const double& timefacfacpre,
          const double& timefacfacprerhs) const;

      //! Compute element matrix and rhs entries: inertia and reactive terms of the Galerkin part
      void CalcDivEps(const LINALG::Matrix<nsd_, nen_>& eveln,
          const LINALG::Matrix<nsd_, nen_>& evelnp, const LINALG::Matrix<nsd_, nen_>& efluidvelnp,
          const LINALG::Matrix<nsd_, nen_>& efluidvelnpp);

      //! Compute element matrix entries: PSPG
      void DiscretePSPG(LINALG::Matrix<nen_ * nsd_, nen_>& estif_w_q,
          LINALG::Matrix<nen_, nen_>& estif_r_q, LINALG::Matrix<nsd_, nen_>& velforce,
          LINALG::Matrix<nen_, 1>& preforce,
          const LINALG::Matrix<nsd_ * nsd_, nen_>& GalMomTestStat, const double& timefacfac,
          const double& timefacfacrhs, const double& timefacfacpre,
          const double& timefacfacprerhs) const;

      //! Compute element matrix entries: SUPG
      void DiscreteSUPG(LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>&
                            estif_w_v,                   ///< block (weighting function v x u)
          LINALG::Matrix<nen_, nen_ * nsd_>& estif_r_v,  ///< block (weighting function v x p)
          LINALG::Matrix<nsd_, nen_>& velforce,          ///< rhs forces velocity
          LINALG::Matrix<nen_, 1>& preforce,
          const LINALG::Matrix<nsd_ * nsd_, nen_>&
              GalMomTestStat,             ///< linearisation of the stabilization residual
          const double& timefacfac,       ///< = timefac x fac
          const double& timefacfacrhs,    ///< = timefacrhs x fac
          const double& timefacfacpre,    ///< = timefacpre x fac
          const double& timefacfacprerhs  ///< = timefacprerhs x fac
          ) const;

      void DiscreteContStab(LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_w_v,
          LINALG::Matrix<nsd_, nen_>& velforce, const double& timefacfacdiv,
          const double& timefacfacdivrhs) const;

      //! Compute element matrix entries: PSPG
      void PSPG(LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          LINALG::Matrix<nen_, nen_>& ppmat,                   ///< block (weighting function q x p)
          LINALG::Matrix<nen_, 1>& preforce,                   ///< rhs forces pressure
          const LINALG::Matrix<nsd_ * nsd_, nen_>&
              GalMomResnU,  ///< linearisation of the stabilization residual
          const LINALG::Matrix<nsd_, 1>& StrongResMomScaled,
          const double& timefacfac,       ///< = timefac x fac
          const double& timefacfacrhs,    ///< = timefacrhs x fac
          const double& timefacfacpre,    ///< = timefacpre x fac
          const double& timefacfacprerhs  ///< = timefacprerhs x fac
          ) const;

      //! Compute element matrix entries: SUPG
      void SUPG(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,       ///< block (weighting function v x p)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          const LINALG::Matrix<nsd_ * nsd_, nen_>&
              GalMomResnU,  ///< linearisation of the stabilization residual
          const LINALG::Matrix<nsd_, 1>& StrongResMomScaled,
          const double& timefacfac,       ///< = timefac x fac
          const double& timefacfacrhs,    ///< = timefacrhs x fac
          const double& timefacfacpre,    ///< = timefacpre x fac
          const double& timefacfacprerhs  ///< = timefacprerhs x fac
          ) const;

      /// compute residuum of continuity equation
      void ContRes(double& StrongResContScaled,
          const double& timefacfacdiv,    ///< = timefac x fac
          const double& timefacfacdivrhs  ///< = timefac x fac
          ) const;

      //! Compute element matrix entries: div-grad stabilization
      void ContStab(
          LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<nsd_, nen_>& velforce,               ///< rhs forces velocity
          const double& timefacfacdiv,                        ///< = timefac x fac
          const double& timefacfacdivrhs                      ///< = timefac x fac
          ) const;

      /*!
       * \brief fill elment matrix and vectors with the global values
       */
      void ExtractValuesFromGlobalVector(
          const DRT::Discretization& discretization,  ///< discretization
          const std::vector<int>& lm,                 ///<
          LINALG::Matrix<nsd_, nen_>* matrixtofill,   ///< vector field
          LINALG::Matrix<nen_, 1>* vectortofill,      ///< scalar field
          const std::string state                     ///< state of the global vector
          ) const;



      //! @name coordinate data

      //! node coordinates
      LINALG::Matrix<nsd_, nen_> xyze_;
      //@}


      //! @name shape functions

      //! array for shape functions
      LINALG::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives w.r.t r,s,t
      LINALG::Matrix<nsd_, nen_> deriv_;
      //! array for second derivatives of shape function w.r.t r,s,t
      LINALG::Matrix<numderiv2_, nen_> deriv2_;
      //! transposed jacobian "dx/ds"
      LINALG::Matrix<nsd_, nsd_> xjm_;
      //! inverse of transposed jacobian "ds/dx"
      LINALG::Matrix<nsd_, nsd_> xji_;
      //! global derivative of shape functions w.r.t x,y,z
      LINALG::Matrix<nsd_, nen_> derxy_;
      //! global second derivatives of shape functions w.r.t x,y,z
      LINALG::Matrix<numderiv2_, nen_> derxy2_;
      //! shape function for viscous term
      LINALG::Matrix<nsd_ * nsd_, nen_> visc_shp_;
      //@}


      //! @name variable values at last iteration (at time t^n)

      //! velocity vector in gausspoint
      LINALG::Matrix<nsd_, 1> velint_;
      //! global velocity derivatives in gausspoint w.r.t x,y,z
      LINALG::Matrix<nsd_, nsd_> vderxy_;
      //! pressure in gausspoint
      double pres_;
      //! pressure gradient in gausspoint
      LINALG::Matrix<nsd_, 1> gradp_;
      //! viscous velocity in gausspoint
      LINALG::Matrix<nsd_, 1> visc_;
      //! fluid velocity vector in gausspoint
      LINALG::Matrix<nsd_, 1> fluidvelint_;
      //! global fluid derivative of velocity w.r.t x,y,z
      LINALG::Matrix<nsd_, nsd_> fluidvelxy_;
      //! fluid pressure in gausspoint
      double fluidpres_;
      //! fluid pressure gradient in gausspoint
      LINALG::Matrix<nsd_, 1> fluidgradp_;
      //! viscous fluid velocity in gausspoint
      LINALG::Matrix<nsd_, 1> fluidvisc_;
      //! body force of momentum equation
      LINALG::Matrix<nsd_, 1> fluidbodyforce_;
      //! body force of adjoint momentum equation
      LINALG::Matrix<nsd_, 1> bodyforce_;
      //! body force of adjoint continuity equation
      double contforce_;
      //! velocity divergence
      double vdiv_;
      //! u_j v_i,j
      LINALG::Matrix<nsd_, 1> conv1_;
      //! u_j,i * v_j
      LINALG::Matrix<nsd_, 1> conv2_;

      //@}

      //! @name variable values at last time step (= t^n+1)

      //! velocity vector in gausspoint
      LINALG::Matrix<nsd_, 1> velint_old_;
      //! global velocity derivatives in gausspoint w.r.t x,y,z
      LINALG::Matrix<nsd_, nsd_> vderxy_old_;
      //! pressure in gausspoint
      double pres_old_;
      //! pressure gradient in gausspoint
      LINALG::Matrix<nsd_, 1> gradp_old_;
      //! viscous velocity in gausspoint
      LINALG::Matrix<nsd_, 1> visc_old_;
      //! fluid velocity vector in gausspoint
      LINALG::Matrix<nsd_, 1> fluidvelint_old_;
      //! global fluid derivative of velocity w.r.t x,y,z
      LINALG::Matrix<nsd_, nsd_> fluidvelxy_old_;
      //! body force of adjoint momentum equation
      LINALG::Matrix<nsd_, 1> bodyforce_old_;
      //! body force of adjoint continuity equation
      double contforce_old_;
      //! velocity divergence
      double vdiv_old_;
      //! u_j v_i,j
      LINALG::Matrix<nsd_, 1> conv1_old_;
      //! u_j,i * v_j
      LINALG::Matrix<nsd_, 1> conv2_old_;
      //@}


      //! @name variable values at next time step (= t^n-1)

      //! fluid velocity vector in gausspoint
      LINALG::Matrix<nsd_, 1> fluidvelint_new_;
      //! global fluid derivative of velocity w.r.t x,y,z
      LINALG::Matrix<nsd_, nsd_> fluidvelxy_new_;
      //! fluid pressure gradient in gausspoint
      LINALG::Matrix<nsd_, 1> fluidgradp_new_;
      //! viscous fluid velocity in gausspoint
      LINALG::Matrix<nsd_, 1> fluidvisc_new_;
      //! body force of momentum equation
      LINALG::Matrix<nsd_, 1> fluidbodyforce_new_;
      //@}


      //! the stabilisation parameters -> it is a (3,1) vector for 2D and 3D
      LINALG::Matrix<3, 1> tau_;
      LINALG::Matrix<3, 1> tau_old_;


      //! @name gauss point data

      // gauss points
      DRT::UTILS::GaussIntegration intpoints_;
      //! coordinates of current integration point in reference coordinates
      LINALG::Matrix<nsd_, 1> xsi_;
      //! Jacobian determinant
      double det_;
      //! integration factor
      double fac_;
      //@}


      //! @name physical coefficients

      //! reaction coefficient
      double reacoeff_;
      //@}

      //! element pointer
      DRT::ELEMENTS::Fluid* ele_;
      //! Flag to (de)activate higher order elements
      //! elements with only mixed second order derivatives are not counted as higher order elements
      //! (see definition of higher order elements in fluid3_ele_impl_utils.cpp)
      bool is_higher_order_ele_;

      //! pointer to parameter list
      Teuchos::RCP<DRT::ELEMENTS::FluidAdjoint3ImplParameter> fldAdPara_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* TOPOPT_FLUIDADJOINT3_IMPL_H_ */
