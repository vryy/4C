/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of poro Fluid element (p1, mixed poro fluid)


\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_CALC_PORO_P1_HPP
#define FOUR_C_FLUID_ELE_CALC_PORO_P1_HPP

#include "baci_config.hpp"

#include "baci_fluid_ele_calc_poro.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  class StructPoro;
}

namespace DRT
{
  namespace ELEMENTS
  {
    //! Class for Evaluating boundary integrals for porous media problems
    /*!
     This class is derived from the FluidEleCalc class, i.e. it is capable of
     evaluated all integrals implemented there. It will do so, if the evaluate action
     given by the control routine is not known (see method Evaluate).

     The main methods are the Evaluate and the EvaluateOD routines. Therein,
     the stiffness matrixes of a porous fluid problem are evaluated. OD means
     off diagonal, indicating linearizations with respect to structural degrees of freedom,
     that will be assembled into off diagonal entries in the global system matrix.
     The terms are eventually evaluated in the GaussPointLoop.. methods

     This a calculation class implemented as a singleton, like all calc classes in fluid
     (see comments on base classes for more details). In short this means that on instance
     exists for every discretization type of the boundary element (because of the template).

     This is the poro P1 implementation, i.e. meant to be coupled with a structure problem
     which solves for the porosity. For the fluid not much changes compared with the
     standard implementation. Only the porosity is evaluated in a different way.

     \author vuong 10/14
     */
    template <CORE::FE::CellType distype>
    class FluidEleCalcPoroP1 : public FluidEleCalcPoro<distype>
    {
      using Base = FluidEleCalcPoro<distype>;

     protected:
      using Base::nen_;
      using Base::nsd_;

      //! private Constructor since we are a Singleton.
      FluidEleCalcPoroP1();

     public:
      //! Singleton access method
      static FluidEleCalcPoroP1<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      /*!
      \brief calculate element matrix and rhs for porous flow (2)

      \param eid              (i) element id
      \param discretization   (i) fluid discretization the element belongs to
      \param lm               (i) location matrix of element
      \param params           (i) element parameter list
      \param mat              (i) material
      \param elemat1_epetra   (o) element matrix to calculate
      \param elemat2_epetra   (o) element matrix to calculate
      \param elevec1_epetra   (o) element vector to calculate
      \param elevec2_epetra   (o) element vector to calculate
      \param elevec3_epetra   (o) element vector to calculate
      \param intpoints        (i) Gaussian integration points

      */
      int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          const CORE::FE::GaussIntegration& intpoints) override;

      /*!
      \brief Evaluate coupling terms (off diagonal terms) of the element at specified gauss points
      for porous flow

      \param eid              (i) element id
      \param discretization   (i) fluid discretization the element belongs to
      \param lm               (i) location matrix of element
      \param params           (i) element parameter list
      \param mat              (i) material
      \param elemat1_epetra   (o) element matrix to calculate
      \param elemat2_epetra   (o) element matrix to calculate
      \param elevec1_epetra   (o) element vector to calculate
      \param elevec2_epetra   (o) element vector to calculate
      \param elevec3_epetra   (o) element vector to calculate
      \param intpoints        (i) Gaussian integration points

      */
      int EvaluateOD(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          const CORE::FE::GaussIntegration& intpoints) override;

     protected:
      /*!
        \brief evaluate function for Fluid element for porous flow

        Specific evaluate function without any knowledge about DRT objects. This
        way the element evaluation is independent of the specific mesh storage.

            \param params           (i) element parameter list
            \param elemat1          (o) element matrix to be filled
            \param elevec1          (o) element rhs vector to be filled
            \param evelaf           (i) nodal velocities at n+alpha_F/n+1
            \param epreaf           (i) nodal pressure at n+alpha_F/n+1
            \param evelnp           (i) nodal velocities at n+1 (np_genalpha)
            \param eprenp           (i) nodal pressure at n+alpha_F/n+1
            \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
            \param eaccam           (i) nodal accelerations at n+alpha_M
            \param edispnp          (i) nodal displacements at n+1 (on moving mesh)
            \param egridv           (i) grid velocity at n+1
            \param escaaf           (i) nodal scalar at n+alpha_F/n+1
            \param material         (i) fluid material
            \param isale            (i) ALE flag
            \param intpoints        (i) Gaussian integration points
       */
      int EvaluateOD(Teuchos::ParameterList& params,
          const CORE::LINALG::Matrix<nsd_, nen_>& ebofoaf,
          CORE::LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& elemat1,
          CORE::LINALG::Matrix<(nsd_ + 1) * nen_, 1>& elevec1,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelaf,
          const CORE::LINALG::Matrix<nen_, 1>& epreaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelnp,
          const CORE::LINALG::Matrix<nsd_, nen_>& eveln,
          const CORE::LINALG::Matrix<nen_, 1>& eprenp, const CORE::LINALG::Matrix<nen_, 1>& epren,
          const CORE::LINALG::Matrix<nsd_, nen_>& emhist,
          const CORE::LINALG::Matrix<nen_, 1>& echist,
          const CORE::LINALG::Matrix<nen_, 1>& epressnp_timederiv,
          const CORE::LINALG::Matrix<nen_, 1>& epressam_timederiv,
          const CORE::LINALG::Matrix<nen_, 1>& epressn_timederiv,
          const CORE::LINALG::Matrix<nsd_, nen_>& eaccam,
          const CORE::LINALG::Matrix<nsd_, nen_>& edispnp,
          const CORE::LINALG::Matrix<nsd_, nen_>& edispn,
          const CORE::LINALG::Matrix<nsd_, nen_>& egridv,
          const CORE::LINALG::Matrix<nsd_, nen_>& egridvn,
          const CORE::LINALG::Matrix<nen_, 1>& escaaf,
          const CORE::LINALG::Matrix<nen_, 1>* eporositynp, Teuchos::RCP<MAT::Material> mat,
          bool isale, const CORE::FE::GaussIntegration& intpoints);

      /*!
        \brief calculate off diagonal element matrix and rhs for porous flow

        \param params           (i) element parameter list
        \param evelaf           (i) nodal velocities at n+alpha_F/n+1
        \param evelnp           (i) nodal velocities at n+1 (np_genalpha)
        \param epreaf           (i) nodal pressure at n+alpha_F/n+1
        \param eprenp           (i) nodal pressure at n+alpha_F/n+1
        \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
        \param eaccam           (i) nodal accelerations at n+alpha_M
        \param edispnp          (i) nodal displacements at n+1 (on moving mesh)
        \param egridv           (i) grid velocity at n+1
        \param escaaf           (i) nodal scalar at n+alpha_F/n+1
        \param ecoupl           (o) element matrix to calculate
        \param eforce           (o) element rhs to calculate
        \param material         (i) fluid material
        \param isale            (i) ALE flag
        \param intpoints        (i) Gaussian integration points
      */
      void SysmatOD(Teuchos::ParameterList& params, const CORE::LINALG::Matrix<nsd_, nen_>& ebofoaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelnp,
          const CORE::LINALG::Matrix<nsd_, nen_>& eveln,
          const CORE::LINALG::Matrix<nen_, 1>& epreaf, const CORE::LINALG::Matrix<nen_, 1>& eprenp,
          const CORE::LINALG::Matrix<nen_, 1>& epren,
          const CORE::LINALG::Matrix<nsd_, nen_>& emhist,
          const CORE::LINALG::Matrix<nen_, 1>& echist,
          const CORE::LINALG::Matrix<nen_, 1>& epressnp_timederiv,
          const CORE::LINALG::Matrix<nen_, 1>& epressam_timederiv,
          const CORE::LINALG::Matrix<nen_, 1>& epressn_timederiv,
          const CORE::LINALG::Matrix<nsd_, nen_>& eaccam,
          const CORE::LINALG::Matrix<nsd_, nen_>& edispnp,
          const CORE::LINALG::Matrix<nsd_, nen_>& edispn,
          const CORE::LINALG::Matrix<nsd_, nen_>& egridv,
          const CORE::LINALG::Matrix<nsd_, nen_>& egridvn,
          const CORE::LINALG::Matrix<nen_, 1>& escaaf,
          const CORE::LINALG::Matrix<nen_, 1>* eporositynp,
          CORE::LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& ecoupl,
          CORE::LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce,
          Teuchos::RCP<const MAT::Material> material, bool isale,
          const CORE::FE::GaussIntegration& intpoints);

      /*!
        \brief Gauss point loop for evaluation of off-diagonal terms

        \param params           ( i) element parameter list
        \param evelaf             (i) nodal velocities at n+alpha_F/n+1
        \param evelnp             (i) nodal velocities at n+1 (np_genalpha)
        \param epreaf             (i) nodal pressure at n+alpha_F/n+1
        \param eprenp             (i) nodal pressure at n+alpha_F/n+1
        \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
        \param eaccam             (i) nodal accelerations at n+alpha_M
        \param edispnp            (i) nodal displacements at n+1 (on moving mesh)
        \param egridv             (i) grid velocity at n+1
        \param escaaf             (i) nodal scalar at n+alpha_F/n+1
        \param eporositynp        (i) nodal porosity at n+alpha_F/n+1
        \param eforce             (o) coupling rhs force vector
        \param ecoupl_u           (o) coupling element matrix for fluid velocity
        \param ecoupl_p           (o) coupling element matrix for fluid pressure
        \param ecouplp1_u         (o) coupling element matrix for fluid velocity (porosity dependent
        terms) \param ecouplp1_p         (o) coupling element matrix for fluid pressure (porosity
        dependent terms) \param material           (i) fluid material \param intpoints          (i)
        Gaussian integration points
      */
      virtual void GaussPointLoopP1OD(Teuchos::ParameterList& params,
          const CORE::LINALG::Matrix<nsd_, nen_>& ebofoaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelnp,
          const CORE::LINALG::Matrix<nsd_, nen_>& eveln,
          const CORE::LINALG::Matrix<nen_, 1>& epreaf, const CORE::LINALG::Matrix<nen_, 1>& eprenp,
          const CORE::LINALG::Matrix<nen_, 1>& epren,
          const CORE::LINALG::Matrix<nsd_, nen_>& emhist,
          const CORE::LINALG::Matrix<nen_, 1>& echist,
          const CORE::LINALG::Matrix<nen_, 1>& epressnp_timederiv,
          const CORE::LINALG::Matrix<nen_, 1>& epressam_timederiv,
          const CORE::LINALG::Matrix<nen_, 1>& epressn_timederiv,
          const CORE::LINALG::Matrix<nsd_, nen_>& eaccam,
          const CORE::LINALG::Matrix<nsd_, nen_>& edispnp,
          const CORE::LINALG::Matrix<nsd_, nen_>& edispn,
          const CORE::LINALG::Matrix<nsd_, nen_>& egridv,
          const CORE::LINALG::Matrix<nsd_, nen_>& egridvn,
          const CORE::LINALG::Matrix<nen_, 1>& escaaf,
          const CORE::LINALG::Matrix<nen_, 1>* eporositynp,
          CORE::LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce,
          CORE::LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& ecoupl_u,
          CORE::LINALG::Matrix<nen_, nen_ * nsd_>& ecoupl_p,
          CORE::LINALG::Matrix<nen_ * nsd_, nen_>& ecouplp1_u,
          CORE::LINALG::Matrix<nen_, nen_>& ecouplp1_p, Teuchos::RCP<const MAT::Material> material,
          const CORE::FE::GaussIntegration& intpoints);

      /*!
       \brief evaluate pressure equation (i.e. continuity equation for standard poro elements)
              This function is overwritten by the poro_p1 element

        \param params         (i) element parameter list
        \param timefacfacpre  (i) fluid pressure at gauss point
        \param rhsfac         (i) jacobian determinant at gauss point
        \param dphi_dp        (i) derivative of porosity gradient w.r.t. fluid pressure
        \param dphi_dJ        (i) derivative of porosity gradient w.r.t. jacobian determinant
        \param dphi_dJdp      (i) mixed derivative of porosity gradient w.r.t. fluid pressure and
       jacobian determinant \param dphi_dJJ       (i) second derivative of porosity gradient w.r.t.
       jacobian determinant \param dphi_dpp       (i) second derivative of porosity gradient w.r.t.
       fluid pressure \param eporositydot   (i) nodal values of porosity time derivative at time
       step n+1 \param eporositydotn  (i) nodal values of porosity time derivative at time step n
        \param echist         (i) nodal values of history values of continuity equation
        \param dgradphi_dp    (i) derivative of porosity gradient w.r.t. fluid pressure
        \param estif_q_u      (o) element matrix (pressure - fluid velocity weighting) to be filled
        \param ppmat          (o) element matrix (pressure - pressure weighting) to be filled
        \param preforce       (o) element rhs vector to be filled
       * */
      void EvaluatePressureEquation(Teuchos::ParameterList& params, const double& timefacfacpre,
          const double& rhsfac, const double& dphi_dp, const double& dphi_dJ,
          const double& dphi_dJdp, const double& dphi_dpp,
          const CORE::LINALG::Matrix<nen_, 1>* eporositydot,
          const CORE::LINALG::Matrix<nen_, 1>* eporositydotn,
          const CORE::LINALG::Matrix<nen_, 1>& echist,
          const CORE::LINALG::Matrix<nsd_, nen_>& dgradphi_dp,
          CORE::LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u,
          CORE::LINALG::Matrix<nen_, nen_>& ppmat,
          CORE::LINALG::Matrix<nen_, 1>& preforce) override;

      /*!
        \brief Compute porosity and derivatives

        \param params         (i) element parameter list
        \param press          (i) fluid pressure at gauss point
        \param J              (i) jacobian determinant at gauss point
        \param gp             (i) number of actual gauss point
        \param shapfct        (i) shape function values at gauss point
        \param myporosity     (i) nodal porosities
        \param porosity       (o) porosity at gauss point
        \param dphi_dp        (o) derivative of porosity gradient w.r.t. fluid pressure
        \param dphi_dJ        (o) derivative of porosity gradient w.r.t. jacobian determinant
        \param dphi_dJdp      (o) mixed derivative of porosity gradient w.r.t. fluid pressure and
        jacobian determinant \param dphi_dJJ       (o) second derivative of porosity gradient w.r.t.
        jacobian determinant \param dphi_dpp       (o) second derivative of porosity gradient w.r.t.
        fluid pressure \param save           (i) flag for saving porosity within structure material
      */
      void ComputePorosity(Teuchos::ParameterList& params, const double& press, const double& J,
          const int& gp, const CORE::LINALG::Matrix<nen_, 1>& shapfct,
          const CORE::LINALG::Matrix<nen_, 1>* myporosity, double& porosity, double* dphi_dp,
          double* dphi_dJ, double* dphi_dJdp, double* dphi_dJJ, double* dphi_dpp,
          bool save) override;

      /*!
        \brief Compute spatial gradient of porosity

        \param dphidp         (i) derivative of porosity w.r.t. fluid pressure
        \param dphidJ         (i) derivative of porosity w.r.t. jacobian determinant
        \param eporositynp    (i) nodal porosities at n+1
        \param gradJ          (i) spatial gradient of jacobian determinant
        \param grad_porosity  (o) spatial gradient of porosity
      */
      void ComputePorosityGradient(const double& dphidp, const double& dphidJ,
          const CORE::LINALG::Matrix<nsd_, 1>& gradJ, const CORE::LINALG::Matrix<nsd_, 1>& gradp,
          const CORE::LINALG::Matrix<nen_, 1>* eporositynp,
          CORE::LINALG::Matrix<nsd_, 1>& grad_porosity,
          CORE::LINALG::Matrix<nsd_, 1>& refgrad_porosity) override;

      /*!
        \brief Compute linearization needed for diagonal terms (lin. of porosity gradient w.r.t.
        fluid pressure)

        \param dphidp         (i) derivative of porosity w.r.t. fluid pressure
        \param dphi_dpp       (i) second derivative of porosity w.r.t. fluid pressure
        \param dphi_dJp       (i) mixed derivative of porosity w.r.t. fluid pressure and jacobian
        determinant \param gradJ          (i) spatial gradient of jacobian determinant \param
        dgradphi_dp    (o) derivate of spatial gradient of porosity w.r.t. fluid pressure
       */
      void ComputeLinearization(const double& dphi_dp, const double& dphi_dpp,
          const double& dphi_dJp, const CORE::LINALG::Matrix<nsd_, 1>& gradJ,
          CORE::LINALG::Matrix<nsd_, nen_>& dgradphi_dp) override;

      /*!
        \brief Compute linearization needed for off diagonal terms
          (lin. of jacobian, porosity and porosity gradient w.r.t. structure displacement)

        \param dphi_dJ        (i) derivative of porosity w.r.t. jacobian determinant
        \param dphi_dJJ       (i) second derivative of porosity w.r.t. jacobian determinant
        \param dphi_dJp       (i) mixed derivative of porosity w.r.t. fluid pressure and jacobian
        determinant \param defgrd_inv     (i) inverse deformation gradient \param defgrd_IT_vec  (o)
        inverse transposed deformation gradient in vector notation \param F_x            (i)
        derivative of deformation gradient w.r.t. current coordinates xyz \param F_X            (i)
        derivative of deformation gradient w.r.t. material coordinates XYZ \param gradJ          (i)
        spatial gradient of jacobian determinant \param dJ_dus         (o) derivative of jacobian
        determinant w.r.t. structure displacments \param dphi_dus       (o) derivative of porosity
        determinant w.r.t. structure displacments \param dgradphi_dus   (o) derivate of spatial
        gradient of porosity w.r.t. structure displacments
       */
      void ComputeLinearizationOD(const double& dphi_dJ, const double& dphi_dJJ,
          const double& dphi_dJp, const CORE::LINALG::Matrix<nsd_, nsd_>& defgrd_inv,
          const CORE::LINALG::Matrix<nsd_ * nsd_, 1>& defgrd_IT_vec,
          const CORE::LINALG::Matrix<nsd_ * nsd_, nsd_>& F_x,
          const CORE::LINALG::Matrix<nsd_ * nsd_, nsd_>& F_X,
          const CORE::LINALG::Matrix<nsd_, 1>& gradJ, CORE::LINALG::Matrix<1, nsd_ * nen_>& dJ_dus,
          CORE::LINALG::Matrix<1, nsd_ * nen_>& dphi_dus,
          CORE::LINALG::Matrix<nsd_, nen_ * nsd_>& dgradphi_dus) override;

      //! Compute element matrix entries: PSPG
      void PSPG(
          CORE::LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u,  //!< block (weighting function q x u)
          CORE::LINALG::Matrix<nen_, nen_>& ppmat,             //!< block (weighting function q x p)
          CORE::LINALG::Matrix<nen_, 1>& preforce,             //!< rhs forces pressure
          const CORE::LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,  //!< linearisation of the stabilization residual
          const CORE::LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resMRea_Du,
          const CORE::LINALG::Matrix<nsd_, nen_>&
              lin_resM_Dp,        //!< linearisation of the stabilization residual w.r.t. pressure
          const double& dphi_dp,  //!< linearisation of porosity w.r.t. pressure
          const double& fac3,     //!< factor for residual in current subgrid velocities
          const double& timefacfac,     //!< = timefac x fac
          const double& timefacfacpre,  //!< = timefacpre x fac
          const double& rhsfac          //!< right-hand-side factor for residuals
          ) override;

      //! Compute element matrix entries: reactive stabilization
      void ReacStab(CORE::LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>&
                        estif_u,                               //!< block (weighting function v x u)
          CORE::LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,  //!< block (weighting function v x p)
          CORE::LINALG::Matrix<nsd_, nen_>& velforce,          //!< rhs forces velocity
          CORE::LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,  //!< linearisation of the stabilization residual
          const CORE::LINALG::Matrix<nsd_, nen_>&
              lin_resM_Dp,        //!< linearisation of the stabilization residual w.r.t. pressure
          const double& dphi_dp,  //!< linearisation of porosity w.r.t. pressure
          const double& timefacfac,     //!< = timefac x fac
          const double& timefacfacpre,  //!< = timefacpre x fac
          const double& rhsfac,         //!< right-hand-side factor for residuals
          const double& fac3            //!< factor for residual in current subgrid velocities
          ) override;

      int ComputeVolume(Teuchos::ParameterList& params,
          DRT::ELEMENTS::Fluid* ele,                //!< current fluid element
          DRT::Discretization& discretization,      //!< fluid discretization
          std::vector<int>& lm,                     //!< location vector for DOF management
          CORE::LINALG::SerialDenseVector& elevec1  //!< reference to element vector to be filled
          ) override;
    };
  }  // namespace ELEMENTS
}  // namespace DRT


FOUR_C_NAMESPACE_CLOSE

#endif
