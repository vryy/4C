/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_loma.H

\brief Low-Mach-number flow routines for calculation of fluid element

\level 2

<pre>
\maintainer Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_LOMA_H
#define FLUID_ELE_CALC_LOMA_H

#include "fluid_ele_interface.H"
#include "fluid_ele_calc.H"

// includes gauss points
#include "../drt_lib/drt_utils.H"

namespace DRT
{
  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype>
    class FluidEleCalcLoma : public FluidEleCalc<distype>
    {
      typedef DRT::ELEMENTS::FluidEleCalc<distype> my;

     public:
      virtual ~FluidEleCalcLoma() {}

      /// Singleton access method
      static FluidEleCalcLoma<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

      virtual int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          bool offdiag = false);

      /// Evaluate the element at specified gauss points for porous flow
      virtual int EvaluateOD(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          const DRT::UTILS::GaussIntegration& intpoints);

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
          const LINALG::Matrix<my::nsd_, my::nen_>& ebofoaf,
          const LINALG::Matrix<my::nsd_, my::nen_>& eprescpgaf,
          LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nen_>& elemat1,
          const LINALG::Matrix<my::nsd_, my::nen_>& evelaf,
          const LINALG::Matrix<my::nen_, 1>& epreaf, const LINALG::Matrix<my::nen_, 1>& epream,
          const LINALG::Matrix<my::nen_, 1>& escaaf,
          const LINALG::Matrix<my::nsd_, my::nen_>& emhist,
          const LINALG::Matrix<my::nsd_, my::nen_>& eaccam,
          const LINALG::Matrix<my::nen_, 1>& escadtam,
          const LINALG::Matrix<my::nen_, 1>& escabofoaf,
          const LINALG::Matrix<my::nsd_, my::nen_>& eveln,
          const LINALG::Matrix<my::nen_, 1>& escaam,
          const LINALG::Matrix<my::nsd_, my::nen_>& edispnp,
          const LINALG::Matrix<my::nsd_, my::nen_>& egridv, Teuchos::RCP<MAT::Material> mat,
          bool isale, double CsDeltaSq, double CiDeltaSq,
          const DRT::UTILS::GaussIntegration& intpoints);

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
      void SysmatOD(const LINALG::Matrix<my::nsd_, my::nen_>& ebofoaf,
          const LINALG::Matrix<my::nsd_, my::nen_>& eprescpgaf,
          const LINALG::Matrix<my::nsd_, my::nen_>& evelaf,
          const LINALG::Matrix<my::nsd_, my::nen_>& eveln,
          const LINALG::Matrix<my::nen_, 1>& epreaf, const LINALG::Matrix<my::nen_, 1>& epream,
          const LINALG::Matrix<my::nsd_, my::nen_>& eaccam,
          const LINALG::Matrix<my::nen_, 1>& escaaf, const LINALG::Matrix<my::nen_, 1>& escaam,
          const LINALG::Matrix<my::nen_, 1>& escadtam,
          const LINALG::Matrix<my::nen_, 1>& escabofoaf,
          const LINALG::Matrix<my::nsd_, my::nen_>& emhist,
          const LINALG::Matrix<my::nsd_, my::nen_>& edispnp,
          const LINALG::Matrix<my::nsd_, my::nen_>& egridv,
          LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nen_>& estif, const double thermpressaf,
          const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
          Teuchos::RCP<const MAT::Material> material, double& Cs_delta_sq, double& Ci_delta_sq,
          bool isale, const DRT::UTILS::GaussIntegration& intpoints);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif
