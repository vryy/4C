/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 2

 */
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_PORO_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_PORO_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_boundary_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcPoro : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      typedef Discret::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim> my;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      /// Singleton access method
      static ScaTraEleBoundaryCalcPoro<distype, probdim>* instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate action
      int evaluate_action(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, ScaTra::BoundaryAction action,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra) override;

     protected:
      //! compute integral of convective mass/heat flux over boundary surface
      std::vector<double> calc_convective_flux(const Core::Elements::FaceElement* ele,
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ephinp,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelnp,
          Core::LinAlg::SerialDenseVector& erhs) override;

      //! compute porosity based on solid, fluid and (potentially) scatra solution
      virtual double compute_porosity(
          const Core::Elements::FaceElement* ele  //!< the element we are dealing with
      );

     private:
      /// private constructor since we are singleton
      ScaTraEleBoundaryCalcPoro(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! nodal porosity values at t_(n+1)
      Core::LinAlg::Matrix<nen_, 1> eporosity_;

      //! nodal pressure values at t_(n+1) or t_(n+alpha_F)
      Core::LinAlg::Matrix<nen_, 1> eprenp_;

      //! flag indacting a node based porosity
      bool isnodalporosity_;

    };  // class ScaTraEleBoundaryCalcPoro
  }     // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
