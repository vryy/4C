/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for low Mach number problems


\level 2
 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_LOMA_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_LOMA_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_boundary_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    // class implementation
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcLoma : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      typedef Discret::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim> my;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcLoma<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);


      //! evaluate action
      int evaluate_action(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, ScaTra::BoundaryAction action,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra) override;

     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcLoma(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate loma thermal press
      void calc_loma_therm_press(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la);

      //! calculate Neumann inflow boundary conditions
      void neumann_inflow(const Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& emat, Core::LinAlg::SerialDenseVector& erhs) override;

      //! integral of normal diffusive flux and velocity over boundary surface
      void norm_diff_flux_and_vel_integral(const Core::Elements::Element* ele,
          Teuchos::ParameterList& params, const std::vector<double>& enormdiffflux,
          const std::vector<double>& enormvel);

      //! thermodynamic pressure
      double thermpress_;
    };
  }  // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
