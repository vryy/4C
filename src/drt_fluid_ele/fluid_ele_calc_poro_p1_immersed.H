/*----------------------------------------------------------------------*/
/*! \file

\brief internal implementation of poro immersed fluid element (p1 poro fluid)

\maintainer Martin Kronbichler

\level 3

*/
/*----------------------------------------------------------------------*/


#ifndef FLUID_ELE_CALC_PORO_P1_IMMERSED_H_
#define FLUID_ELE_CALC_PORO_P1_IMMERSED_H_

#include "fluid_ele_calc_poro_p1.H"


namespace DRT
{
  class ImmersedFieldExchangeManager;

  namespace ELEMENTS
  {
    class FluidImmersedBase;

    /*!
      \brief Class for evaluating modifications terms in case of immersed problems.

      If a deformable structure is immersed in the porous medium, we have to make sure, that
      the compressibility of this structure is compatible with the compressibilities of the
      phases in the porous medium. the fluid phase on the porous medium is supposed to be
      incompressible, i.e., div(v^f)=0. In our immersed method, the velocities of the immersed
      structure are projected to the nodes of the porous medium. If the projected velocity field
      is not divergence free, this would result in spurious oscillations, unreasonable pressure
      solutions, and deterioration of fluid mass conservation. In this class we modify the
     continuity equation of the porous medium for the artificial porous medium domain, i.e., the
     porous medium elements covered by the immersed structure. In addition, we modifiy the equation
     for each integration point in the poro elements intersecting the immersed boundary.

     \author rauch 02/17
     */
    template <DRT::Element::DiscretizationType distype>
    class FluidEleCalcPoroP1Immersed : public FluidEleCalcPoroP1<distype>
    {
      typedef FluidEleCalc<distype> my;
      typedef FluidEleCalcPoro<distype> my_p;
      typedef FluidEleCalcPoroP1<distype> my_p1;

     protected:
      /// private Constructor since we are a Singleton.
      FluidEleCalcPoroP1Immersed();

      /// pointer to exchange manager
      ImmersedFieldExchangeManager* exchange_manager_;

     public:
      virtual ~FluidEleCalcPoroP1Immersed() {}

      /// Singleton access method
      static FluidEleCalcPoroP1Immersed<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();


     protected:
      /*!
        \brief Calls evaluate of base class with different number of integration points on elements
        which are IsBoundaryImmered() (ngp defined in input file)

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
        \param offdiag          (i) flag indicating whether diagonal or off diagonal blocks are to
        be calculated

       */
      int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          bool offdiag = false);

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
          const int& gp, const LINALG::Matrix<my::nen_, 1>& shapfct,
          const LINALG::Matrix<my::nen_, 1>* myporosity, double& porosity, double* dphi_dp,
          double* dphi_dJ, double* dphi_dJdp, double* dphi_dJJ, double* dphi_dpp, bool save);

      /*!
       \brief evaluate pressure equation (i.e. continuity equation for standard poro elements)
              This function is overwritten by the poro_p1 elements

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
          const LINALG::Matrix<my::nen_, 1>* eporositydot,
          const LINALG::Matrix<my::nen_, 1>* eporositydotn,
          const LINALG::Matrix<my::nen_, 1>& echist,
          const LINALG::Matrix<my::nsd_, my::nen_>& dgradphi_dp,
          LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>& estif_q_u,
          LINALG::Matrix<my::nen_, my::nen_>& ppmat, LINALG::Matrix<my::nen_, 1>& preforce);

      /*!
       \brief evaluate nontransient part of pressure equation (i.e. no time derivative of porosity)

        \param params         (i) element parameter list
        \param timefacfacpre  (i) fluid pressure at gauss point
        \param rhsfac         (i) jacobian determinant at gauss point
        \param dphi_dp        (i) derivative of porosity gradient w.r.t. fluid pressure
        \param dphi_dJ        (i) derivative of porosity gradient w.r.t. jacobian determinant
        \param dphi_dJdp      (i) mixed derivative of porosity gradient w.r.t. fluid pressure and
       jacobian determinant \param dphi_dJJ       (i) second derivative of porosity gradient w.r.t.
       jacobian determinant \param dphi_dpp       (i) second derivative of porosity gradient w.r.t.
       fluid pressure \param echist         (i) nodal values of history values of continuity
       equation \param dgradphi_dp    (i) derivative of porosity gradient w.r.t. fluid pressure
        \param estif_q_u      (o) element matrix (pressure - fluid velocity weighting) to be filled
        \param ppmat          (o) element matrix (pressure - pressure weighting) to be filled
        \param preforce       (o) element rhs vector to be filled
       * */
      void EvaluatePressureEquationNonTransient(Teuchos::ParameterList& params,
          const double& timefacfacpre, const double& rhsfac, const double& dphi_dp,
          const double& dphi_dJ, const double& dphi_dJdp, const double& dphi_dpp,
          const LINALG::Matrix<my::nsd_, my::nen_>& dgradphi_dp,
          LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>& estif_q_u,
          LINALG::Matrix<my::nen_, my::nen_>& ppmat, LINALG::Matrix<my::nen_, 1>& preforce);

      /*!
        \brief Gauss point loop for evaluation of off-diagonal terms

        \param params             (i) element parameter list
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
          const LINALG::Matrix<my::nsd_, my::nen_>& ebofoaf,
          const LINALG::Matrix<my::nsd_, my::nen_>& evelaf,
          const LINALG::Matrix<my::nsd_, my::nen_>& evelnp,
          const LINALG::Matrix<my::nsd_, my::nen_>& eveln,
          const LINALG::Matrix<my::nen_, 1>& epreaf, const LINALG::Matrix<my::nen_, 1>& eprenp,
          const LINALG::Matrix<my::nen_, 1>& epren,
          const LINALG::Matrix<my::nsd_, my::nen_>& emhist,
          const LINALG::Matrix<my::nen_, 1>& echist,
          const LINALG::Matrix<my::nen_, 1>& epressnp_timederiv,
          const LINALG::Matrix<my::nen_, 1>& epressam_timederiv,
          const LINALG::Matrix<my::nen_, 1>& epressn_timederiv,
          const LINALG::Matrix<my::nsd_, my::nen_>& eaccam,
          const LINALG::Matrix<my::nsd_, my::nen_>& edispnp,
          const LINALG::Matrix<my::nsd_, my::nen_>& edispn,
          const LINALG::Matrix<my::nsd_, my::nen_>& egridv,
          const LINALG::Matrix<my::nsd_, my::nen_>& egridvn,
          const LINALG::Matrix<my::nen_, 1>& escaaf, const LINALG::Matrix<my::nen_, 1>* eporositynp,
          LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>& eforce,
          LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>& ecoupl_u,
          LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>& ecoupl_p,
          LINALG::Matrix<my::nen_ * my::nsd_, my::nen_>& ecouplp1_u,
          LINALG::Matrix<my::nen_, my::nen_>& ecouplp1_p,
          Teuchos::RCP<const MAT::Material> material,
          const DRT::UTILS::GaussIntegration& intpoints);

      /// current immersed element
      DRT::ELEMENTS::FluidImmersedBase* immersedele_;
      /// id of current gp
      int gp_iquad_;
      /// modification type of flow equation in artifical domain
      int artifical_flow_mod_type;
      /*!
        \brief multiplicator is 0 for simplified formulation and 1 for full formulation

          use full formulation in case you expect porosity gradients and temporal
          changes in porosity which are not negligible.
        */
      double mod_type_multiplicator_;

    };  // class FluidEleCalcPoroP1Immersed
  }     // namespace ELEMENTS
}  // namespace DRT


#endif /* FLUID_ELE_CALC_PORO_P1_IMMERSED_H_ */
