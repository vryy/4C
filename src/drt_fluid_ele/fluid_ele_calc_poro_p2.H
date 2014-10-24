/*----------------------------------------------------------------------*/
/*!
 \file fluid_ele_calc_poro_p2.H

 \brief Internal implementation of poro Fluid element (p2 poro fluid)

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#ifndef FLUID_ELE_CALC_PORO_P2_H_
#define FLUID_ELE_CALC_PORO_P2_H_

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

 This is the poro P2 implementation, i.e. meant to be coupled with a scatra problem
 which solves for the porosity. The main difference here is, that the last equation
 evaluted by the fluid is not(!) the balance of mass/continuity equation, but
 the constitutive relation giving the fluid pressure.

 \author vuong 10/14
 */

template<DRT::Element::DiscretizationType distype>
class FluidEleCalcPoroP2 : public FluidEleCalcPoro<distype>
{
  typedef FluidEleCalcPoro<distype> my;

private:

  /// private Constructor since we are a Singleton.
  FluidEleCalcPoroP2();

public:

    virtual ~FluidEleCalcPoroP2() {}

    /// Singleton access method
    static FluidEleCalcPoroP2<distype>* Instance( bool create=true );

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
          \param eporositynp      (i) nodal porosity at n+1
          \param mat              (i) fluid material
          \param isale            (i) ALE flag
          \param intpoints        (i) Gaussian integration points
     */
    int EvaluateODPoroScatra(
        Teuchos::ParameterList&                                           params,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nen_> &             elemat1,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> &                    elevec1,
        const LINALG::Matrix<my::nsd_,my::nen_> &                         evelaf,
        const LINALG::Matrix<my::nen_, 1> &                               epreaf,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        evelnp,
        const LINALG::Matrix<my::nen_, 1> &                               eprenp,
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
    void SysmatPoroScatraOD(
        Teuchos::ParameterList&                                         params,
        const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
        const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
        const LINALG::Matrix<my::nen_, 1>&                              epreaf,
        const LINALG::Matrix<my::nen_, 1>&                              eprenp,
        const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
        const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
        const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
        const LINALG::Matrix<my::nen_,1>&                               escaaf,
        const LINALG::Matrix<my::nen_,1>*                               eporositynp,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_,  my::nen_>&           ecoupl,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
        Teuchos::RCP<const MAT::Material>                               material,
        bool                                                            isale,
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
      \brief Evaluate off diagonal terms in continuity equation

      \param timefacfacpre      (i) time factor (pressure) * integration factor
      \param porosity           (i) porosity at gauss point
      \param dphi_dJ            (i) derivative of porosity gradient w.r.t. jacobian determinant at gauss point
      \param dphi_dJJ           (i) second derivative of porosity gradient w.r.t. jacobian determinant at gauss point
      \param dphi_dJdp          (i) mixed derivative of porosity gradient w.r.t. jacobian and pressure at gauss point
      \param gridvelint         (i) grid (structure) velocity at gauss point
      \param dgradphi_dus       (i) derivative of porosity gradient w.r.t. structural displacements at gauss point
      \param dphi_dus           (i) derivative of porosity w.r.t. structural displacements at gauss point
      \param dJ_dus             (i) derivative of jacobian w.r.t. structural displacements at gauss point
      \param ecoupl_u           (o) coupling element matrix of continuity equation
    */
    virtual void FillMatrixContiOD(
        const double&                                               timefacfacpre,
        const double &                                              dphi_dp,
        const double &                                              dphi_dJ,
        const double&                                               dphi_dJJ,
        const double&                                               dphi_dJdp,
        const double &                                              refporositydot,
        const LINALG::Matrix<my::nsd_,my::nen_*my::nsd_>&           dgradphi_dus,
        const LINALG::Matrix<1,my::nsd_*my::nen_>&                  dphi_dus,
        const LINALG::Matrix<1,my::nsd_*my::nen_>&                  dJ_dus,
        const LINALG::Matrix<my::nsd_, my::nen_>&                   egridv,
        const LINALG::Matrix<my::nsd_, my::nen_ * my::nsd_>&        lin_resM_Dus,
        LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>&              ecoupl_p);

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
        const LINALG::Matrix<my::nen_,1>*                  eporositynp);

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

    //! internal energy due to porosity change at gausspoint
    double W_;
    //! derivative of internal energy due to porosity change w.r.t. Jacobian at gausspoint
    double dW_dJ_;

};

}//ELEMENTS
}//DRT


#endif /* FLUID_ELE_CALC_PORO_P2_H_ */
