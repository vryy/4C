/*----------------------------------------------------------------------*/
/*! \file
 \brief implementation of evaluation routines of porofluidmultiphase boundary element

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_ELE_BOUNDARY_CALC_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_ELE_BOUNDARY_CALC_HPP


#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_lib_element.hpp"
#include "baci_porofluidmultiphase_ele_action.hpp"
#include "baci_porofluidmultiphase_ele_calc_utils.hpp"
#include "baci_porofluidmultiphase_ele_interface.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    class PoroFluidMultiPhaseEleParameter;

    /*!
    \brief implementation of evaluation routines of porous fluid multiphase boundary element

    This singleton class is responsible for evaluating boundary terms.
    It provides the method Evaluate(...) which performs the actual evaluation
    depending on the action provided by the global algorithm.


    \author vuong
    */
    template <CORE::FE::CellType distype>
    class PoroFluidMultiPhaseEleBoundaryCalc : public PoroFluidMultiPhaseEleInterface
    {
     public:
      /// Singleton access method
      static PoroFluidMultiPhaseEleBoundaryCalc<distype>* Instance(
          const int numdofpernode, const std::string& disname);


      //! number of element nodes (nomenclature: T. Hughes, The finite element method)
      static constexpr int nen_ = CORE::FE::num_nodes<distype>;

      //! number of boundary(!) space dimensions
      static constexpr int nsd_ = CORE::FE::dim<distype>;

      //! Evaluate the element (using location array)
      int Evaluate(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          std::vector<CORE::LINALG::SerialDenseMatrix*>& elemat,
          std::vector<CORE::LINALG::SerialDenseVector*>& elevec) override;

     protected:
      //! setup element evaluation
      int SetupCalc(
          DRT::Element* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization);

      //! extract element based or nodal values
      //  return extracted values of phinp
      virtual void ExtractElementAndNodeValues(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la);

      //! evaluate action
      virtual int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, POROFLUIDMULTIPHASE::BoundaryAction action,
          DRT::Element::LocationArray& la, std::vector<CORE::LINALG::SerialDenseMatrix*>& elemat,
          std::vector<CORE::LINALG::SerialDenseVector*>& elevec);

      //! evaluate Neumann boundary condition
      virtual int EvaluateNeumann(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseVector& elevec1);

      //! evaluate shape functions and derivatives at int. point
      double EvalShapeFuncAndIntFac(
          const CORE::FE::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
          const int iquad,                                       ///< id of current Gauss point
          CORE::LINALG::Matrix<1 + nsd_, 1>* normalvec =
              nullptr  ///< normal vector at Gauss point(optional)
      );

     private:
      /// private constructor since we are singleton
      PoroFluidMultiPhaseEleBoundaryCalc(const int numdofpernode, const std::string& disname);

      //! pointer to parameter list
      DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter* params_;

      //! number of dof per node
      const int numdofpernode_;

      //! node coordinates
      CORE::LINALG::Matrix<nsd_ + 1, nen_> xyze_;
      //! nodal displacement values for ALE
      CORE::LINALG::Matrix<nsd_ + 1, nen_> edispnp_;
      //! coordinates of current integration point in reference coordinates
      CORE::LINALG::Matrix<nsd_, 1> xsi_;
      //! array for shape functions
      CORE::LINALG::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives w.r.t r,s,t
      CORE::LINALG::Matrix<nsd_, nen_> deriv_;
      //! global derivatives of shape functions w.r.t x,y,z
      CORE::LINALG::Matrix<nsd_, nen_> derxy_;
      //! unit normal vector at integration point
      CORE::LINALG::Matrix<nsd_ + 1, 1> normal_;
      //! velocity vector in gausspoint
      CORE::LINALG::Matrix<nsd_ + 1, 1> velint_;
      //! metric tensor at integration point
      CORE::LINALG::Matrix<nsd_, nsd_> metrictensor_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT



FOUR_C_NAMESPACE_CLOSE

#endif
