/*----------------------------------------------------------------------*/
/*!

\brief calc class for immersed problems

\maintainer Martin Kronbichler

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_IMMERSED_H
#define FLUID_ELE_CALC_IMMERSED_H

#include "fluid_ele_calc.H"


namespace DRT
{
  namespace ELEMENTS
  {
    class FluidImmersedBase;

    template <DRT::Element::DiscretizationType distype>
    class FluidEleCalcImmersed : public FluidEleCalc<distype>
    {
      typedef DRT::ELEMENTS::FluidEleCalc<distype> my;

     protected:
      /// private Constructor since we are a Singleton.
      FluidEleCalcImmersed();

     public:
      virtual ~FluidEleCalcImmersed() {}

      /// Singleton access method
      static FluidEleCalcImmersed<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();


     protected:
      /*!
        Evaluate

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
      virtual int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          bool offdiag = false);

      //! compute residual of momentum equation and subgrid-scale velocity
      void ComputeSubgridScaleVelocity(
          const LINALG::Matrix<my::nsd_, my::nen_>& eaccam,  ///< acceleration at time n+alpha_M
          double& fac1,                                      ///< factor for old s.-s. velocities
          double& fac2,                                      ///< factor for old s.-s. accelerations
          double& fac3,     ///< factor for residual in current s.-s. velocities
          double& facMtau,  ///< facMtau = modified tau_M (see code)
          int iquad,        ///< integration point
          double* saccn,    ///< s.-s. acceleration at time n+alpha_a / n
          double* sveln,    ///< s.-s. velocity at time n+alpha_a / n
          double* svelnp    ///< s.-s. velocity at time n+alpha_f / n+1
      );

      //! Provide linearization of Garlerkin momentum residual with respect to the velocities
      void LinGalMomResU(LINALG::Matrix<my::nsd_ * my::nsd_, my::nen_>&
                             lin_resM_Du,  ///< linearisation of the Garlerkin momentum residual
          const double& timefacfac         ///< = timefac x fac
      );

      //! Compute element matrix and rhs entries: inertia, convective andyn
      //! reactive terms of the Galerkin part
      void InertiaConvectionReactionGalPart(
          LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>&
              estif_u,                                   ///< block (weighting function v x u)
          LINALG::Matrix<my::nsd_, my::nen_>& velforce,  ///< rhs forces velocity
          LINALG::Matrix<my::nsd_ * my::nsd_, my::nen_>&
              lin_resM_Du,  ///< linearisation of the Garlerkin momentum residual
          LINALG::Matrix<my::nsd_, 1>&
              resM_Du,          ///< linearisation of the Garlerkin momentum residual
          const double& rhsfac  ///< right-hand-side factor
      );

      //! Compute element matrix entries: continuity terms of the Garlerkin part and rhs
      void ContinuityGalPart(LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>&
                                 estif_q_u,       ///< block (weighting function q x u)
          LINALG::Matrix<my::nen_, 1>& preforce,  ///< rhs forces pressure
          const double& timefacfac,               ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac  ///< right-hand-side factor
      );

      //! Compute element matrix entries: conservative formulation
      void ConservativeFormulation(LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>&
                                       estif_u,          ///< block (weighting function v x u)
          LINALG::Matrix<my::nsd_, my::nen_>& velforce,  ///< rhs forces velocity
          const double& timefacfac,                      ///< = timefac x fac
          const double& rhsfac                           ///< right-hand-side factor
      );

      // current element
      DRT::ELEMENTS::FluidImmersedBase* immersedele_;
      // number of current gp
      int gp_iquad_;

    };  // class FluidEleCalcImmersed
  }     // namespace ELEMENTS
}  // namespace DRT

#endif
