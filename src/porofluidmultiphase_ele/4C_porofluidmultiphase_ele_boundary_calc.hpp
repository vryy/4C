/*----------------------------------------------------------------------*/
/*! \file
 \brief implementation of evaluation routines of porofluidmultiphase boundary element

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_ELE_BOUNDARY_CALC_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_ELE_BOUNDARY_CALC_HPP


#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_porofluidmultiphase_ele_action.hpp"
#include "4C_porofluidmultiphase_ele_calc_utils.hpp"
#include "4C_porofluidmultiphase_ele_interface.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    class PoroFluidMultiPhaseEleParameter;

    /*!
    \brief implementation of evaluation routines of porous fluid multiphase boundary element

    This singleton class is responsible for evaluating boundary terms.
    It provides the method evaluate(...) which performs the actual evaluation
    depending on the action provided by the global algorithm.


    \author vuong
    */
    template <Core::FE::CellType distype>
    class PoroFluidMultiPhaseEleBoundaryCalc : public PoroFluidMultiPhaseEleInterface
    {
     public:
      /// Singleton access method
      static PoroFluidMultiPhaseEleBoundaryCalc<distype>* Instance(
          const int numdofpernode, const std::string& disname);


      //! number of element nodes (nomenclature: T. Hughes, The finite element method)
      static constexpr int nen_ = Core::FE::num_nodes<distype>;

      //! number of boundary(!) space dimensions
      static constexpr int nsd_ = Core::FE::dim<distype>;

      //! Evaluate the element (using location array)
      int evaluate(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Elements::Element::LocationArray& la,
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,
          std::vector<Core::LinAlg::SerialDenseVector*>& elevec) override;

     protected:
      //! setup element evaluation
      int setup_calc(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization);

      //! extract element based or nodal values
      //  return extracted values of phinp
      virtual void extract_element_and_node_values(Core::Elements::Element* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Elements::Element::LocationArray& la);

      //! evaluate action
      virtual int evaluate_action(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, POROFLUIDMULTIPHASE::BoundaryAction action,
          Core::Elements::Element::LocationArray& la,
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,
          std::vector<Core::LinAlg::SerialDenseVector*>& elevec);

      //! evaluate Neumann boundary condition
      virtual int evaluate_neumann(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
          Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseVector& elevec1);

      //! evaluate shape functions and derivatives at int. point
      double eval_shape_func_and_int_fac(
          const Core::FE::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
          const int iquad,                                       ///< id of current Gauss point
          Core::LinAlg::Matrix<1 + nsd_, 1>* normalvec =
              nullptr  ///< normal vector at Gauss point(optional)
      );

     private:
      /// private constructor since we are singleton
      PoroFluidMultiPhaseEleBoundaryCalc(const int numdofpernode, const std::string& disname);

      //! pointer to parameter list
      Discret::ELEMENTS::PoroFluidMultiPhaseEleParameter* params_;

      //! number of dof per node
      const int numdofpernode_;

      //! node coordinates
      Core::LinAlg::Matrix<nsd_ + 1, nen_> xyze_;
      //! nodal displacement values for ALE
      Core::LinAlg::Matrix<nsd_ + 1, nen_> edispnp_;
      //! coordinates of current integration point in reference coordinates
      Core::LinAlg::Matrix<nsd_, 1> xsi_;
      //! array for shape functions
      Core::LinAlg::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives w.r.t r,s,t
      Core::LinAlg::Matrix<nsd_, nen_> deriv_;
      //! global derivatives of shape functions w.r.t x,y,z
      Core::LinAlg::Matrix<nsd_, nen_> derxy_;
      //! unit normal vector at integration point
      Core::LinAlg::Matrix<nsd_ + 1, 1> normal_;
      //! velocity vector in gausspoint
      Core::LinAlg::Matrix<nsd_ + 1, 1> velint_;
      //! metric tensor at integration point
      Core::LinAlg::Matrix<nsd_, nsd_> metrictensor_;
    };
  }  // namespace ELEMENTS
}  // namespace Discret



FOUR_C_NAMESPACE_CLOSE

#endif
