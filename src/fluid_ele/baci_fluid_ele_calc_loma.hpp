/*----------------------------------------------------------------------*/
/*! \file

\brief Low-Mach-number flow routines for calculation of fluid element


\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_CALC_LOMA_HPP
#define FOUR_C_FLUID_ELE_CALC_LOMA_HPP

#include "baci_config.hpp"

#include "baci_fluid_ele_calc.hpp"
#include "baci_fluid_ele_interface.hpp"
#include "baci_lib_utils.hpp"
#include "baci_utils_singleton_owner.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype>
    class FluidEleCalcLoma : public FluidEleCalc<distype>
    {
      typedef DRT::ELEMENTS::FluidEleCalc<distype> my;

      using my::nen_;
      using my::nsd_;

     public:
      /// Singleton access method
      static FluidEleCalcLoma<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra, bool offdiag = false) override;

      /// Evaluate the element at specified gauss points for porous flow
      virtual int EvaluateOD(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra,
          const CORE::FE::GaussIntegration& intpoints);

     private:
      /// private Constructor since we are a Singleton.
      FluidEleCalcLoma();

      /*!
          \brief evaluation of off-diagonal matrix block for monolithic loma solver for fluid3
         element

          Specific evaluate function without any knowledge about DRT objects. This
          way the element evaluation is independent of the specific mesh storage.
       */
      int EvaluateOD(Teuchos::ParameterList& params,
          const CORE::LINALG::Matrix<nsd_, nen_>& ebofoaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& eprescpgaf,
          CORE::LINALG::Matrix<(nsd_ + 1) * nen_, nen_>& elemat1,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelaf,
          const CORE::LINALG::Matrix<nen_, 1>& epreaf, const CORE::LINALG::Matrix<nen_, 1>& epream,
          const CORE::LINALG::Matrix<nen_, 1>& escaaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& emhist,
          const CORE::LINALG::Matrix<nsd_, nen_>& eaccam,
          const CORE::LINALG::Matrix<nen_, 1>& escadtam,
          const CORE::LINALG::Matrix<nen_, 1>& escabofoaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& eveln,
          const CORE::LINALG::Matrix<nen_, 1>& escaam,
          const CORE::LINALG::Matrix<nsd_, nen_>& edispnp,
          const CORE::LINALG::Matrix<nsd_, nen_>& egridv, Teuchos::RCP<MAT::Material> mat,
          bool isale, double CsDeltaSq, double CiDeltaSq,
          const CORE::FE::GaussIntegration& intpoints);

      /*!
          \brief calculate element matrix for off-diagonal matrix block for monolithic
         low-Mach-number solver

          \param ebofoaf          (i) body force at n+alpha_F/n+1
          \param eprescpgaf       (i) prescribed pressure gradient at n+alpha_F/n+1 (required for
         turbulent channel flow) \param evelaf           (i) nodal velocities at n+alpha_F/n+1
          \param epreaf           (i) nodal pressure at n+alpha_F/n+1
          \param eaccam           (i) nodal accelerations at n+alpha_M
          \param escaaf           (i) nodal scalar at n+alpha_F/n+1
          \param escaam           (i) nodal scalar at n+alpha_M/n
          \param emhist           (i) time rhs for momentum equation
          \param edispnp          (i) nodal displacements (on moving mesh)
          \param egridv           (i) grid velocity (on moving mesh)
          \param estif            (o) element matrix to calculate
          \param thermpressaf     (i) thermodynamic pressure at n+alpha_F/n+1
          \param thermpressam     (i) thermodynamic pressure at n+alpha_M/n
          \param thermpressdtaf   (i) thermodynamic pressure derivative at n+alpha_F/n+1
          \param thermpressdtam   (i) thermodynamic pressure derivative at n+alpha_M/n+1
          \param material         (i) fluid material
          \param Cs_delta_sq      (i) parameter for dynamic Smagorinsky model (Cs*h*h)
          \param isale            (i) ALE flag
          \param intpoints        (i) Gaussian integration points

       */
      void SysmatOD(const CORE::LINALG::Matrix<nsd_, nen_>& ebofoaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& eprescpgaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& eveln,
          const CORE::LINALG::Matrix<nen_, 1>& epreaf, const CORE::LINALG::Matrix<nen_, 1>& epream,
          const CORE::LINALG::Matrix<nsd_, nen_>& eaccam,
          const CORE::LINALG::Matrix<nen_, 1>& escaaf, const CORE::LINALG::Matrix<nen_, 1>& escaam,
          const CORE::LINALG::Matrix<nen_, 1>& escadtam,
          const CORE::LINALG::Matrix<nen_, 1>& escabofoaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& emhist,
          const CORE::LINALG::Matrix<nsd_, nen_>& edispnp,
          const CORE::LINALG::Matrix<nsd_, nen_>& egridv,
          CORE::LINALG::Matrix<(nsd_ + 1) * nen_, nen_>& estif, const double thermpressaf,
          const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
          Teuchos::RCP<const MAT::Material> material, double& Cs_delta_sq, double& Ci_delta_sq,
          bool isale, const CORE::FE::GaussIntegration& intpoints);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
