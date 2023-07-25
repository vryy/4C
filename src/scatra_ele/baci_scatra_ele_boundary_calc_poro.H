/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 2

 */
/*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_ELE_BOUNDARY_CALC_PORO_H
#define BACI_SCATRA_ELE_BOUNDARY_CALC_PORO_H

#include "baci_scatra_ele_boundary_calc.H"

namespace DRT
{
  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype,
        int probdim = CORE::DRT::UTILS::DisTypeToDim<distype>::dim + 1>
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
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra);

     protected:
      //! compute integral of convective mass/heat flux over boundary surface
      virtual std::vector<double> CalcConvectiveFlux(const DRT::FaceElement* ele,
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephinp,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelnp, Epetra_SerialDenseVector& erhs);

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
#endif
