/*----------------------------------------------------------------------*/
/*! \file

\brief calc class for immersed problems


\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_FLUID_ELE_CALC_IMMERSED_HPP
#define BACI_FLUID_ELE_CALC_IMMERSED_HPP

#include "baci_config.hpp"

#include "baci_fluid_ele_calc.hpp"
#include "baci_utils_singleton_owner.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidImmersedBase;

    template <CORE::FE::CellType distype>
    class FluidEleCalcImmersed : public FluidEleCalc<distype>
    {
      typedef DRT::ELEMENTS::FluidEleCalc<distype> my;
      using my::nen_;
      using my::nsd_;

     protected:
      /// private Constructor since we are a Singleton.
      FluidEleCalcImmersed();

     public:
      /// Singleton access method
      static FluidEleCalcImmersed<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);


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
      int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra, bool offdiag = false) override;

      //! compute residual of momentum equation and subgrid-scale velocity
      void ComputeSubgridScaleVelocity(
          const CORE::LINALG::Matrix<nsd_, nen_>& eaccam,  ///< acceleration at time n+alpha_M
          double& fac1,                                    ///< factor for old s.-s. velocities
          double& fac2,                                    ///< factor for old s.-s. accelerations
          double& fac3,     ///< factor for residual in current s.-s. velocities
          double& facMtau,  ///< facMtau = modified tau_M (see code)
          int iquad,        ///< integration point
          double* saccn,    ///< s.-s. acceleration at time n+alpha_a / n
          double* sveln,    ///< s.-s. velocity at time n+alpha_a / n
          double* svelnp    ///< s.-s. velocity at time n+alpha_f / n+1
          ) override;

      //! Provide linearization of Garlerkin momentum residual with respect to the velocities
      void LinGalMomResU(CORE::LINALG::Matrix<nsd_ * nsd_, nen_>&
                             lin_resM_Du,  ///< linearisation of the Garlerkin momentum residual
          const double& timefacfac         ///< = timefac x fac
          ) override;

      //! Compute element matrix and rhs entries: inertia, convective andyn
      //! reactive terms of the Galerkin part
      void InertiaConvectionReactionGalPart(CORE::LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>&
                                                estif_u,  ///< block (weighting function v x u)
          CORE::LINALG::Matrix<nsd_, nen_>& velforce,     ///< rhs forces velocity
          CORE::LINALG::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,  ///< linearisation of the Garlerkin momentum residual
          CORE::LINALG::Matrix<nsd_, 1>&
              resM_Du,          ///< linearisation of the Garlerkin momentum residual
          const double& rhsfac  ///< right-hand-side factor
          ) override;

      //! Compute element matrix entries: continuity terms of the Garlerkin part and rhs
      void ContinuityGalPart(
          CORE::LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          CORE::LINALG::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          const double& timefacfac,                            ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac  ///< right-hand-side factor
          ) override;

      //! Compute element matrix entries: conservative formulation
      void ConservativeFormulation(CORE::LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>&
                                       estif_u,        ///< block (weighting function v x u)
          CORE::LINALG::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& timefacfac,                    ///< = timefac x fac
          const double& rhsfac                         ///< right-hand-side factor
          ) override;

      // current element
      DRT::ELEMENTS::FluidImmersedBase* immersedele_;
      // number of current gp
      int gp_iquad_;

    };  // class FluidEleCalcImmersed
  }     // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
