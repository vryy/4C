/*----------------------------------------------------------------------*/
/*! \file

\brief Low-Mach-number flow routines for calculation of fluid element


\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_CALC_LOMA_HPP
#define FOUR_C_FLUID_ELE_CALC_LOMA_HPP

#include "4C_config.hpp"

#include "4C_fluid_ele_calc.hpp"
#include "4C_fluid_ele_interface.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    template <Core::FE::CellType distype>
    class FluidEleCalcLoma : public FluidEleCalc<distype>
    {
      typedef Discret::ELEMENTS::FluidEleCalc<distype> my;

      using my::nen_;
      using my::nsd_;

     public:
      /// Singleton access method
      static FluidEleCalcLoma<distype>* Instance(
          Core::UTILS::SingletonAction action = Core::UTILS::SingletonAction::create);

      int Evaluate(Discret::ELEMENTS::Fluid* ele, Discret::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra, bool offdiag = false) override;

      /// Evaluate the element at specified gauss points for porous flow
      virtual int evaluate_od(Discret::ELEMENTS::Fluid* ele,
          Discret::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params, Teuchos::RCP<Core::Mat::Material>& mat,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          const Core::FE::GaussIntegration& intpoints);

     private:
      /// private Constructor since we are a Singleton.
      FluidEleCalcLoma();

      /*!
          \brief evaluation of off-diagonal matrix block for monolithic loma solver for fluid3
         element

          Specific evaluate function without any knowledge about DRT objects. This
          way the element evaluation is independent of the specific mesh storage.
       */
      int evaluate_od(Teuchos::ParameterList& params,
          const Core::LinAlg::Matrix<nsd_, nen_>& ebofoaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& eprescpgaf,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, nen_>& elemat1,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,
          const Core::LinAlg::Matrix<nen_, 1>& epreaf, const Core::LinAlg::Matrix<nen_, 1>& epream,
          const Core::LinAlg::Matrix<nen_, 1>& escaaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& emhist,
          const Core::LinAlg::Matrix<nsd_, nen_>& eaccam,
          const Core::LinAlg::Matrix<nen_, 1>& escadtam,
          const Core::LinAlg::Matrix<nen_, 1>& escabofoaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& eveln,
          const Core::LinAlg::Matrix<nen_, 1>& escaam,
          const Core::LinAlg::Matrix<nsd_, nen_>& edispnp,
          const Core::LinAlg::Matrix<nsd_, nen_>& egridv, Teuchos::RCP<Core::Mat::Material> mat,
          bool isale, double CsDeltaSq, double CiDeltaSq,
          const Core::FE::GaussIntegration& intpoints);

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
      void sysmat_od(const Core::LinAlg::Matrix<nsd_, nen_>& ebofoaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& eprescpgaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& eveln,
          const Core::LinAlg::Matrix<nen_, 1>& epreaf, const Core::LinAlg::Matrix<nen_, 1>& epream,
          const Core::LinAlg::Matrix<nsd_, nen_>& eaccam,
          const Core::LinAlg::Matrix<nen_, 1>& escaaf, const Core::LinAlg::Matrix<nen_, 1>& escaam,
          const Core::LinAlg::Matrix<nen_, 1>& escadtam,
          const Core::LinAlg::Matrix<nen_, 1>& escabofoaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& emhist,
          const Core::LinAlg::Matrix<nsd_, nen_>& edispnp,
          const Core::LinAlg::Matrix<nsd_, nen_>& egridv,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, nen_>& estif, const double thermpressaf,
          const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
          Teuchos::RCP<const Core::Mat::Material> material, double& Cs_delta_sq,
          double& Ci_delta_sq, bool isale, const Core::FE::GaussIntegration& intpoints);
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
