/*!----------------------------------------------------------------------
\file fluid_ele_calc_poro_p1_immersed.H

 \brief internal implementation of poro immersed fluid element (p1 poro fluid)

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

\level 2

*----------------------------------------------------------------------*/


#ifndef FLUID_ELE_CALC_PORO_P1_IMMERSED_H_
#define FLUID_ELE_CALC_PORO_P1_IMMERSED_H_

#include "fluid_ele_calc_poro_p1.H"


namespace DRT
{

class ImmersedFieldExchangeManager;

namespace ELEMENTS
{

class FluidImmersedBase;

template<DRT::Element::DiscretizationType distype>
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
  static FluidEleCalcPoroP1Immersed<distype>* Instance( bool create=true );

  /// called upon destruction
  virtual void Done();


protected:

  /*!
    \brief Calls evaluate of base class with different number of integration points on elements which are IsBoundaryImmered() (ngp defined in input file)

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
    \param offdiag          (i) flag indicating whether diagonal or off diagonal blocks are to be calculated

   */
  int Evaluate(DRT::ELEMENTS::Fluid*         ele,
                       DRT::Discretization &         discretization,
                       const std::vector<int> &      lm,
                       Teuchos::ParameterList&       params,
                       Teuchos::RCP<MAT::Material> & mat,
                       Epetra_SerialDenseMatrix&     elemat1_epetra,
                       Epetra_SerialDenseMatrix&     elemat2_epetra,
                       Epetra_SerialDenseVector&     elevec1_epetra,
                       Epetra_SerialDenseVector&     elevec2_epetra,
                       Epetra_SerialDenseVector&     elevec3_epetra,
                       bool                          offdiag = false);

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
  void ComputePorosity(
      Teuchos::ParameterList&                params,
      const double&                          press,
      const double&                          J,
      const int&                             gp,
      const LINALG::Matrix<my::nen_,1>&      shapfct,
      const LINALG::Matrix<my::nen_,1>*      myporosity,
      double&                                porosity,
      double*                                dphi_dp,
      double*                                dphi_dJ,
      double*                                dphi_dJdp,
      double*                                dphi_dJJ,
      double*                                dphi_dpp,
      bool                                   save);

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
  void EvaluatePressureEquation(
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
   \brief evaluate nontransient part of pressure equation (i.e. no time derivative of porosity)

    \param params         (i) element parameter list
    \param timefacfacpre  (i) fluid pressure at gauss point
    \param rhsfac         (i) jacobian determinant at gauss point
    \param dphi_dp        (i) derivative of porosity gradient w.r.t. fluid pressure
    \param dphi_dJ        (i) derivative of porosity gradient w.r.t. jacobian determinant
    \param dphi_dJdp      (i) mixed derivative of porosity gradient w.r.t. fluid pressure and jacobian determinant
    \param dphi_dJJ       (i) second derivative of porosity gradient w.r.t. jacobian determinant
    \param dphi_dpp       (i) second derivative of porosity gradient w.r.t. fluid pressure
    \param echist         (i) nodal values of history values of continuity equation
    \param dgradphi_dp    (i) derivative of porosity gradient w.r.t. fluid pressure
    \param estif_q_u      (o) element matrix (pressure - fluid velocity weighting) to be filled
    \param ppmat          (o) element matrix (pressure - pressure weighting) to be filled
    \param preforce       (o) element rhs vector to be filled
   * */
  void EvaluatePressureEquationNonTransient(
      Teuchos::ParameterList&                       params,
      const double&                                 timefacfacpre,
      const double&                                 rhsfac,
      const double&                                 dphi_dp,
      const double&                                 dphi_dJ,
      const double&                                 dphi_dJdp,
      const double&                                 dphi_dpp,
      const LINALG::Matrix<my::nsd_,my::nen_>&      dgradphi_dp,
      LINALG::Matrix<my::nen_, my::nen_*my::nsd_>&  estif_q_u,
      LINALG::Matrix<my::nen_,my::nen_>&            ppmat,
      LINALG::Matrix<my::nen_,1>&                   preforce);


  /// current immersed element
  DRT::ELEMENTS::FluidImmersedBase* immersedele_;

  /// id of current gp
  int gp_iquad_;

}; // class FluidEleCalcPoroP1Immersed
}  // namespace ELEMENTS
}  // namespace DRT


#endif /* FLUID_ELE_CALC_PORO_P1_IMMERSED_H_ */
