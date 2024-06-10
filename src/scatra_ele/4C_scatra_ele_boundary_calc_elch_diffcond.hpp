/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for diffusion-conduction formulation


\level 2
 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_ELCH_DIFFCOND_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_ELCH_DIFFCOND_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_boundary_calc_elch_electrode.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleDiffManagerElchDiffCond;
    class ScaTraEleParameterElchDiffCond;

    // class implementation
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcElchDiffCond
        : public ScaTraEleBoundaryCalcElchElectrode<distype, probdim>
    {
      using my = Discret::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim>;
      using myelch = Discret::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>;
      using myelectrode = Discret::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype, probdim>;

     protected:
      using my::nen_;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);



     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcElchDiffCond(
          const int numdofpernode, const int numscal, const std::string& disname);

      int evaluate_action(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, ScaTra::BoundaryAction action,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra) override;

      int evaluate_neumann(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
          Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseVector& elevec1,
          const double scalar) override;

      void evaluate_elch_boundary_kinetics(const Core::Elements::Element* ele,
          Core::LinAlg::SerialDenseMatrix& emat, Core::LinAlg::SerialDenseVector& erhs,
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ephinp,
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ehist, double timefac,
          Teuchos::RCP<const Core::Mat::Material> material,
          Teuchos::RCP<Core::Conditions::Condition> cond, const int nume,
          const std::vector<int> stoich, const int kinetics, const double pot0, const double frt,
          const double scalar) override;

      void evaluate_s2_i_coupling(const Core::Elements::FaceElement* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseMatrix& eslavematrix,
          Core::LinAlg::SerialDenseMatrix& emastermatrix,
          Core::LinAlg::SerialDenseVector& eslaveresidual) override;

      void evaluate_s2_i_coupling_od(const Core::Elements::FaceElement* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& eslavematrix) override;

      double get_valence(
          const Teuchos::RCP<const Core::Mat::Material>& material, const int k) const override;

      //! diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> dmedc_;
    };  // class ScaTraEleBoundaryCalcElchDiffCond
  }     // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
