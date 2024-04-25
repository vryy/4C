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

namespace DRT
{
  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcPoro : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim> my;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      /// Singleton access method
      static ScaTraEleBoundaryCalcPoro<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate action
      int EvaluateAction(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, SCATRA::BoundaryAction action,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

     protected:
      //! compute integral of convective mass/heat flux over boundary surface
      std::vector<double> CalcConvectiveFlux(const DRT::FaceElement* ele,
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephinp,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelnp,
          CORE::LINALG::SerialDenseVector& erhs) override;

      //! compute porosity based on solid, fluid and (potentially) scatra solution
      virtual double ComputePorosity(
          const DRT::FaceElement* ele  //!< the element we are dealing with
      );

     private:
      /// private constructor since we are singleton
      ScaTraEleBoundaryCalcPoro(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! nodal porosity values at t_(n+1)
      CORE::LINALG::Matrix<nen_, 1> eporosity_;

      //! nodal pressure values at t_(n+1) or t_(n+alpha_F)
      CORE::LINALG::Matrix<nen_, 1> eprenp_;

      //! flag indacting a node based porosity
      bool isnodalporosity_;

    };  // class ScaTraEleBoundaryCalcPoro
  }     // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
