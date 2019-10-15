/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of poro Fluid element (p1, mixed poro fluid)

\maintainer Johannes Kremheller

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_PORO_P1_H_
#define FLUID_ELE_CALC_PORO_P1_H_

/*---------------------------------------------------------------------*
 | headers                                                              |
 *---------------------------------------------------------------------*/
#include "fluid_ele_calc_poro.H"

/*---------------------------------------------------------------------*
 | forward declarations                                                 |
 *---------------------------------------------------------------------*/
namespace MAT
{
  class StructPoro;
}


namespace DRT
{
  namespace ELEMENTS
  {
    /// Class for Evaluating boundary integrals for porous media problems
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
    template <DRT::Element::DiscretizationType distype>
    class FluidEleCalcPoroP1 : public FluidEleCalcPoro<distype>
    {
      typedef FluidEleCalcPoro<distype> my;

     protected:
      /// private Constructor since we are a Singleton.
      FluidEleCalcPoroP1();

     public:
      virtual ~FluidEleCalcPoroP1() {}

      /// Singleton access method
      static FluidEleCalcPoroP1<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

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
      virtual int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          const DRT::UTILS::GaussIntegration& intpoints);

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
      virtual int EvaluateOD(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          const DRT::UTILS::GaussIntegration& intpoints);

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
          const LINALG::Matrix<my::nsd_, my::nen_>& ebofoaf,
          LINALG::Matrix<(my::nsd_ + 1) * my::nen_, (my::nsd_ + 1) * my::nen_>& elemat1,
          LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>& elevec1,
          const LINALG::Matrix<my::nsd_, my::nen_>& evelaf,
          const LINALG::Matrix<my::nen_, 1>& epreaf,
          const LINALG::Matrix<my::nsd_, my::nen_>& evelnp,
          const LINALG::Matrix<my::nsd_, my::nen_>& eveln,
          const LINALG::Matrix<my::nen_, 1>& eprenp, const LINALG::Matrix<my::nen_, 1>& epren,
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
          Teuchos::RCP<MAT::Material> mat, bool isale,
          const DRT::UTILS::GaussIntegration& intpoints);

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
      void SysmatOD(Teuchos::ParameterList& params,
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
          LINALG::Matrix<(my::nsd_ + 1) * my::nen_, (my::nsd_ + 1) * my::nen_>& ecoupl,
          LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>& eforce,
          Teuchos::RCP<const MAT::Material> material, bool isale,
          const DRT::UTILS::GaussIntegration& intpoints);

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
      virtual void EvaluatePressureEquation(Teuchos::ParameterList& params,
          const double& timefacfacpre, const double& rhsfac, const double& dphi_dp,
          const double& dphi_dJ, const double& dphi_dJdp, const double& dphi_dpp,
          const LINALG::Matrix<my::nen_, 1>* eporositydot,
          const LINALG::Matrix<my::nen_, 1>* eporositydotn,
          const LINALG::Matrix<my::nen_, 1>& echist,
          const LINALG::Matrix<my::nsd_, my::nen_>& dgradphi_dp,
          LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>& estif_q_u,
          LINALG::Matrix<my::nen_, my::nen_>& ppmat, LINALG::Matrix<my::nen_, 1>& preforce);

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
      virtual void ComputePorosity(Teuchos::ParameterList& params, const double& press,
          const double& J, const int& gp, const LINALG::Matrix<my::nen_, 1>& shapfct,
          const LINALG::Matrix<my::nen_, 1>* myporosity, double& porosity, double* dphi_dp,
          double* dphi_dJ, double* dphi_dJdp, double* dphi_dJJ, double* dphi_dpp, bool save);

      /*!
        \brief Compute spatial gradient of porosity

        \param dphidp         (i) derivative of porosity w.r.t. fluid pressure
        \param dphidJ         (i) derivative of porosity w.r.t. jacobian determinant
        \param eporositynp    (i) nodal porosities at n+1
        \param gradJ          (i) spatial gradient of jacobian determinant
        \param grad_porosity  (o) spatial gradient of porosity
      */
      virtual void ComputePorosityGradient(const double& dphidp, const double& dphidJ,
          const LINALG::Matrix<my::nsd_, 1>& gradJ, const LINALG::Matrix<my::nsd_, 1>& gradp,
          const LINALG::Matrix<my::nen_, 1>* eporositynp,
          LINALG::Matrix<my::nsd_, 1>& grad_porosity,
          LINALG::Matrix<my::nsd_, 1>& refgrad_porosity);

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
          const double& dphi_dJp, const LINALG::Matrix<my::nsd_, 1>& gradJ,
          LINALG::Matrix<my::nsd_, my::nen_>& dgradphi_dp);

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
          const double& dphi_dJp, const LINALG::Matrix<my::nsd_, my::nsd_>& defgrd_inv,
          const LINALG::Matrix<my::nsd_ * my::nsd_, 1>& defgrd_IT_vec,
          const LINALG::Matrix<my::nsd_ * my::nsd_, my::nsd_>& F_x,
          const LINALG::Matrix<my::nsd_ * my::nsd_, my::nsd_>& F_X,
          const LINALG::Matrix<my::nsd_, 1>& gradJ, LINALG::Matrix<1, my::nsd_ * my::nen_>& dJ_dus,
          LINALG::Matrix<1, my::nsd_ * my::nen_>& dphi_dus,
          LINALG::Matrix<my::nsd_, my::nen_ * my::nsd_>& dgradphi_dus);

      //! Compute element matrix entries: PSPG
      virtual void PSPG(LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>&
                            estif_q_u,                ///< block (weighting function q x u)
          LINALG::Matrix<my::nen_, my::nen_>& ppmat,  ///< block (weighting function q x p)
          LINALG::Matrix<my::nen_, 1>& preforce,      ///< rhs forces pressure
          const LINALG::Matrix<my::nsd_ * my::nsd_, my::nen_>&
              lin_resM_Du,  ///< linearisation of the stabilization residual
          const LINALG::Matrix<my::nsd_ * my::nsd_, my::nen_>& lin_resMRea_Du,
          const LINALG::Matrix<my::nsd_, my::nen_>&
              lin_resM_Dp,        ///< linearisation of the stabilization residual w.r.t. pressure
          const double& dphi_dp,  ///< linearisation of porosity w.r.t. pressure
          const double& fac3,     ///< factor for residual in current subgrid velocities
          const double& timefacfac,     ///< = timefac x fac
          const double& timefacfacpre,  ///< = timefacpre x fac
          const double& rhsfac          ///< right-hand-side factor for residuals
      );

      //! Compute element matrix entries: reactive stabilization
      virtual void ReacStab(LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>&
                                estif_u,  ///< block (weighting function v x u)
          LINALG::Matrix<my::nen_ * my::nsd_, my::nen_>&
              estif_p_v,                                 ///< block (weighting function v x p)
          LINALG::Matrix<my::nsd_, my::nen_>& velforce,  ///< rhs forces velocity
          LINALG::Matrix<my::nsd_ * my::nsd_, my::nen_>&
              lin_resM_Du,  ///< linearisation of the stabilization residual
          const LINALG::Matrix<my::nsd_, my::nen_>&
              lin_resM_Dp,        ///< linearisation of the stabilization residual w.r.t. pressure
          const double& dphi_dp,  ///< linearisation of porosity w.r.t. pressure
          const double& timefacfac,     ///< = timefac x fac
          const double& timefacfacpre,  ///< = timefacpre x fac
          const double& rhsfac,         ///< right-hand-side factor for residuals
          const double& fac3            ///< factor for residual in current subgrid velocities
      );

      virtual int ComputeVolume(Teuchos::ParameterList& params,
          DRT::ELEMENTS::Fluid* ele,            //< current fluid element
          DRT::Discretization& discretization,  //< fluid discretization
          std::vector<int>& lm,                 //< location vector for DOF management
          Epetra_SerialDenseVector& elevec1     //< reference to element vector to be filled
      );
    };

  }  // namespace ELEMENTS
}  // namespace DRT


#endif /* FLUID_ELE_CALC_PORO_P1_H_ */
