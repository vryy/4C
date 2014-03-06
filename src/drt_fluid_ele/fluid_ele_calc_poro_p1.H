/*----------------------------------------------------------------------*/
/*!
 \file fluid_ele_calc_poro_p1.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#ifndef FLUID_ELE_CALC_PORO_P1_H_
#define FLUID_ELE_CALC_PORO_P1_H_

#include "fluid_ele_calc_poro.H"

// forward declarations
namespace MAT
{
  class StructPoro;
}


namespace DRT
{
namespace ELEMENTS
{
  class So_Poro_Interface;

template<DRT::Element::DiscretizationType distype>
class FluidEleCalcPoroP1 : public FluidEleCalcPoro<distype>
{
  typedef FluidEleCalcPoro<distype> my;

private:

  /// private Constructor since we are a Singleton.
  FluidEleCalcPoroP1();

public:

    virtual ~FluidEleCalcPoroP1() {}

    /// Singleton access method
    static FluidEleCalcPoroP1<distype>* Instance( bool create=true );

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
    virtual int Evaluate(DRT::ELEMENTS::Fluid*                    ele,
                             DRT::Discretization &                discretization,
                             const std::vector<int> &             lm,
                             Teuchos::ParameterList&              params,
                             Teuchos::RCP<MAT::Material> &        mat,
                             Epetra_SerialDenseMatrix&            elemat1_epetra,
                             Epetra_SerialDenseMatrix&            elemat2_epetra,
                             Epetra_SerialDenseVector&            elevec1_epetra,
                             Epetra_SerialDenseVector&            elevec2_epetra,
                             Epetra_SerialDenseVector&            elevec3_epetra,
                             const DRT::UTILS::GaussIntegration & intpoints);

    /*!
    \brief Evaluate coupling terms (off diagonal terms) of the element at specified gauss points for porous flow

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
    virtual int EvaluateOD(DRT::ELEMENTS::Fluid*                  ele,
                             DRT::Discretization &                discretization,
                             const std::vector<int> &             lm,
                             Teuchos::ParameterList&              params,
                             Teuchos::RCP<MAT::Material> &        mat,
                             Epetra_SerialDenseMatrix&            elemat1_epetra,
                             Epetra_SerialDenseMatrix&            elemat2_epetra,
                             Epetra_SerialDenseVector&            elevec1_epetra,
                             Epetra_SerialDenseVector&            elevec2_epetra,
                             Epetra_SerialDenseVector&            elevec3_epetra,
                             const DRT::UTILS::GaussIntegration & intpoints);

protected:
    /*!
      \brief evaluate function for Fluid element for porous flow

      Specific evaluate function without any knowledge about DRT objects. This
      way the element evaluation is independent of the specific mesh storage.

      \param params           (i) element parameter list
      \param ebofoaf          (i) body force at n+alpha_F/n+1
      \param elemat1          (o) element matrix to be filled
      \param elevec1          (o) element rhs vector to be filled
      \param evelaf           (i) nodal velocities at n+alpha_F/n+1
      \param epreaf           (i) nodal pressure at n+alpha_F/n+1
      \param evelnp           (i) nodal velocities at n+1 (np_genalpha)
      \param eprenp           (i) nodal pressure at n+alpha_F/n+1
      \param epren            (i) nodal pressure at n
      \param emhist           (i) time rhs for momentum equation
      \param echist           (i) time rhs for continuity equation
      \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
      \param eaccam           (i) nodal accelerations at n+alpha_M
      \param edispnp          (i) nodal displacements at n+1 (on moving mesh)
      \param edispn           (i) nodal displacements at n (on moving mesh)
      \param egridv           (i) grid velocity at n+1
      \param escaaf           (i) nodal scalar at n+alpha_F/n+1
      \param material         (i) fluid material
      \param isale            (i) ALE flag
      \param intpoints        (i) Gaussian integration points
     */
    int Evaluate( Teuchos::ParameterList&                                         params,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       ebofoaf,
                  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> &   elemat1,
                  LINALG::Matrix<(my::nsd_+1)*my::nen_,1> &                       elevec1,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       evelaf,
                  const LINALG::Matrix<my::nen_,1>    &                           epreaf,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       evelnp,
                  const LINALG::Matrix<my::nen_,1>    &                           eprenp,
                  const LINALG::Matrix<my::nen_,1>    &                           epren,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       emhist,
                  const LINALG::Matrix<my::nen_,1>&                               echist,
                  const LINALG::Matrix<my::nen_,1>    &                           epressnp_timederiv,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       eaccam,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       edispnp,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       edispn,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       egridv,
                  const LINALG::Matrix<my::nen_,1>&                               escaaf,
                  const LINALG::Matrix<my::nen_,1>&                               eporositynp,
                  Teuchos::RCP<MAT::Material>                                     mat,
                  bool                                                            isale,
                  const DRT::UTILS::GaussIntegration &                            intpoints);

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
          \param edispnp          (i) nodal displacements at n+1 (on moving mesh)
          \param egridv           (i) grid velocity at n+1
          \param escaaf           (i) nodal scalar at n+alpha_F/n+1
          \param material         (i) fluid material
          \param isale            (i) ALE flag
          \param intpoints        (i) Gaussian integration points
     */
    int EvaluateOD(
        Teuchos::ParameterList&                                           params,
        const LINALG::Matrix<my::nsd_,my::nen_> &                         ebofoaf,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, (my::nsd_ + 1) * my::nen_> &  elemat1,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> &                    elevec1,
        const LINALG::Matrix<my::nsd_,my::nen_> &                         evelaf,
        const LINALG::Matrix<my::nen_, 1> &                               epreaf,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        evelnp,
        const LINALG::Matrix<my::nen_, 1> &                               eprenp,
        const LINALG::Matrix<my::nsd_,my::nen_> &                         emhist,
        const LINALG::Matrix<my::nen_,1>&                                 echist,
        const LINALG::Matrix<my::nen_, 1> &                               epressnp_timederiv,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        edispnp,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        egridv,
        const LINALG::Matrix<my::nen_,1>&                                 escaaf,
        const LINALG::Matrix<my::nen_,1>*                                 eporositynp,
        Teuchos::RCP<MAT::Material>                                       mat,
        bool                                                              isale,
        const DRT::UTILS::GaussIntegration &                              intpoints);

    /*!
      \brief calculate off diagonal element matrix and rhs for porous flow

      \param params           (i) element parameter list
      \param evelaf           (i) nodal velocities at n+alpha_F/n+1
      \param evelnp           (i) nodal velocities at n+1 (np_genalpha)
      \param epreaf           (i) nodal pressure at n+alpha_F/n+1
      \param eprenp           (i) nodal pressure at n+alpha_F/n+1
      \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
      \param edispnp          (i) nodal displacements at n+1 (on moving mesh)
      \param egridv           (i) grid velocity at n+1
      \param escaaf           (i) nodal scalar at n+alpha_F/n+1
      \param ecoupl           (o) element matrix to calculate
      \param eforce           (o) element rhs to calculate
      \param material         (i) fluid material
      \param isale            (i) ALE flag
      \param intpoints        (i) Gaussian integration points
    */
    void SysmatOD(
        Teuchos::ParameterList&                                         params,
        const LINALG::Matrix<my::nsd_,my::nen_>&                        ebofoaf,
        const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
        const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
        const LINALG::Matrix<my::nen_, 1>&                              epreaf,
        const LINALG::Matrix<my::nen_, 1>&                              eprenp,
        const LINALG::Matrix<my::nsd_,my::nen_> &                       emhist,
        const LINALG::Matrix<my::nen_,1>&                               echist,
        const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
        const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
        const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
        const LINALG::Matrix<my::nen_,1>&                               escaaf,
        const LINALG::Matrix<my::nen_,1>*                               eporositynp,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_,(my::nsd_ + 1) * my::nen_>&  ecoupl,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
        Teuchos::RCP<const MAT::Material>                               material,
        bool                                                            isale,
        const DRT::UTILS::GaussIntegration &                            intpoints);

    /*!
      \brief Gauss point loop for evaluation of off-diagonal terms

      \param params           ( i) element parameter list
      \param evelaf             (i) nodal velocities at n+alpha_F/n+1
      \param evelnp             (i) nodal velocities at n+1 (np_genalpha)
      \param epreaf             (i) nodal pressure at n+alpha_F/n+1
      \param eprenp             (i) nodal pressure at n+alpha_F/n+1
      \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
      \param edispnp            (i) nodal displacements at n+1 (on moving mesh)
      \param egridv             (i) grid velocity at n+1
      \param escaaf             (i) nodal scalar at n+alpha_F/n+1
      \param eporositynp        (i) nodal porosity at n+alpha_F/n+1
      \param eforce             (o) coupling rhs force vector
      \param ecoupl_u           (o) coupling element matrix for fluid velocity
      \param ecoupl_p           (o) coupling element matrix for fluid pressure
      \param ecouplp1_u         (o) coupling element matrix for fluid velocity (porosity dependent terms)
      \param ecouplp1_p         (o) coupling element matrix for fluid pressure (porosity dependent terms)
      \param material           (i) fluid material
      \param intpoints          (i) Gaussian integration points
    */
    void GaussPointLoopP1OD(
                            Teuchos::ParameterList&                                         params,
                            const LINALG::Matrix<my::nsd_,my::nen_>&                        ebofoaf,
                            const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
                            const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
                            const LINALG::Matrix<my::nen_, 1>&                              epreaf,
                            const LINALG::Matrix<my::nen_, 1>&                              eprenp,
                            const LINALG::Matrix<my::nsd_,my::nen_> &                       emhist,
                            const LINALG::Matrix<my::nen_,1>&                               echist,
                            const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
                            const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
                            const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
                            const LINALG::Matrix<my::nen_,1>&                               escaaf,
                            const LINALG::Matrix<my::nen_,1>*                               eporositynp,
                            LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
                            LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>&       ecoupl_u,
                            LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>&                  ecoupl_p,
                            LINALG::Matrix<my::nen_ * my::nsd_, my::nen_>&                  ecouplp1_u,
                            LINALG::Matrix<my::nen_, my::nen_>&                             ecouplp1_p,
                            Teuchos::RCP<const MAT::Material>                               material,
                            const DRT::UTILS::GaussIntegration &                            intpoints);

    /*!
     \brief evaluate pressure equation (i.e. continuity equation for standard poro elements and constitutive equation for poro_p2)
            This function is overwritten by the poro_p1 and poro_p2 element

      \param params         (i) element parameter list
      \param timefacfacpre  (i) fluid pressure at gauss point
      \param rhsfac         (i) jacobian determinant at gauss point
      \param dphi_dp        (i) derivative of porosity gradient w.r.t. fluid pressure
      \param dphi_dJ        (i) derivative of porosity gradient w.r.t. jacobian determinant
      \param dphi_dJdp      (i) mixed derivative of porosity gradient w.r.t. fluid pressure and jacobian determinant
      \param dphi_dJJ       (i) second derivative of porosity gradient w.r.t. jacobian determinant
      \param dphi_dpp       (i) second derivative of porosity gradient w.r.t. fluid pressure
      \param eporositydot   (i) nodal values of porosity time derivative at time step n+1
      \param eporositydotn  (i) nodal values of porosity time derivative at time step n
      \param echist         (i) nodal values of history values of continuity equation
      \param dgradphi_dp    (i) derivative of porosity gradient w.r.t. fluid pressure
      \param estif_q_u      (o) element matrix (pressure - fluid velocity weighting) to be filled
      \param ppmat          (o) element matrix (pressure - pressure weighting) to be filled
      \param preforce       (o) element rhs vector to be filled
     * */
    virtual void EvaluatePressureEquation(
        Teuchos::ParameterList&                       params,
        const double&                                 timefacfacpre,
        const double&                                 rhsfac,
        const double&                                 dphi_dp,
        const double&                                 dphi_dJ,
        const double&                                 dphi_dJdp,
        const double&                                 dphi_dpp,
        const LINALG::Matrix<my::nen_,1>*             eporositydot,
        const LINALG::Matrix<my::nen_,1>*             eporositydotn,
        const LINALG::Matrix<my::nen_,1>&             echist,
        const LINALG::Matrix<my::nsd_,my::nen_>&      dgradphi_dp,
        LINALG::Matrix<my::nen_, my::nen_*my::nsd_>&  estif_q_u,
        LINALG::Matrix<my::nen_,my::nen_>&            ppmat,
        LINALG::Matrix<my::nen_,1>&                   preforce);

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
      \param dphi_dJdp      (o) mixed derivative of porosity gradient w.r.t. fluid pressure and jacobian determinant
      \param dphi_dJJ       (o) second derivative of porosity gradient w.r.t. jacobian determinant
      \param dphi_dpp       (o) second derivative of porosity gradient w.r.t. fluid pressure
      \param save           (i) flag for saving porosity within structure material
    */
    virtual void ComputePorosity(
                               Teuchos::ParameterList&              params,
                               const double&                        press,
                               const double&                        J,
                               const int&                           gp,
                               const LINALG::Matrix<my::nen_,1>&    shapfct,
                               const LINALG::Matrix<my::nen_,1>*    myporosity,
                               double&                              porosity,
                               double*                              dphi_dp,
                               double*                              dphi_dJ,
                               double*                              dphi_dJdp,
                               double*                              dphi_dJJ,
                               double*                              dphi_dpp,
                               bool                                 save
                               );

    /*!
      \brief Compute spatial gradient of porosity

      \param dphidp         (i) derivative of porosity w.r.t. fluid pressure
      \param dphidJ         (i) derivative of porosity w.r.t. jacobian determinant
      \param eporositynp    (i) nodal porosities at n+1
      \param gradJ          (i) spatial gradient of jacobian determinant
      \param grad_porosity  (o) spatial gradient of porosity
    */
    virtual void ComputePorosityGradient(
        const double&                                      dphidp,
        const double&                                      dphidJ,
        const LINALG::Matrix<my::nsd_,1>&                  gradJ,
        const LINALG::Matrix<my::nen_,1>*                  eporositynp,
        LINALG::Matrix<my::nsd_,1>&                        grad_porosity);

    //! Compute element matrix entries: PSPG
    virtual void PSPG(
        LINALG::Matrix<my::nen_, my::nen_*my::nsd_> &        estif_q_u,      ///< block (weighting function q x u)
        LINALG::Matrix<my::nen_,my::nen_> &                  ppmat,          ///< block (weighting function q x p)
        LINALG::Matrix<my::nen_,1> &                         preforce,       ///< rhs forces pressure
        const LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> &   lin_resM_Du,    ///< linearisation of the stabilization residual
        const LINALG::Matrix<my::nsd_,my::nen_> &            lin_resM_Dp,    ///< linearisation of the stabilization residual w.r.t. pressure
        const double &                                       dphi_dp,        ///< linearisation of porosity w.r.t. pressure
        const double &                                       fac3,           ///< factor for residual in current subgrid velocities
        const double &                                       timefacfac,     ///< = timefac x fac
        const double &                                       timefacfacpre,  ///< = timefacpre x fac
        const double &                                       rhsfac          ///< right-hand-side factor for residuals
        );

    //! Compute element matrix entries: reactive stabilization
    virtual void ReacStab(
        LINALG::Matrix<my::nen_*my::nsd_,my::nen_*my::nsd_> &  estif_u,        ///< block (weighting function v x u)
        LINALG::Matrix<my::nen_*my::nsd_,my::nen_> &           estif_p_v,      ///< block (weighting function v x p)
        LINALG::Matrix<my::nsd_,my::nen_> &                    velforce,       ///< rhs forces velocity
        LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> &           lin_resM_Du,    ///< linearisation of the stabilization residual
        const LINALG::Matrix<my::nsd_,my::nen_> &              lin_resM_Dp,    ///< linearisation of the stabilization residual w.r.t. pressure
        const double &                                         dphi_dp,        ///< linearisation of porosity w.r.t. pressure
        const double &                                         timefacfac,     ///< = timefac x fac
        const double &                                         timefacfacpre,  ///< = timefacpre x fac
        const double &                                         rhsfac,         ///< right-hand-side factor for residuals
        const double &                                         fac3            ///< factor for residual in current subgrid velocities
        );

};

}//ELEMENTS
}//DRT


#endif /* FLUID_ELE_CALC_PORO_P1_H_ */
