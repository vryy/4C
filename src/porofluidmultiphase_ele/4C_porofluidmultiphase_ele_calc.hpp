/*----------------------------------------------------------------------*/
/*! \file
 \brief implementation of the evaluation routines of the porofluidmultiphase element

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_ELE_CALC_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_ELE_CALC_HPP

#include "4C_config.hpp"

#include "4C_fem_general_utils_integration.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_porofluidmultiphase_ele_action.hpp"
#include "4C_porofluidmultiphase_ele_calc_utils.hpp"
#include "4C_porofluidmultiphase_ele_interface.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret
{
  namespace ELEMENTS
  {
    // forward declarations
    class PoroFluidMultiPhaseEleParameter;

    namespace PoroFluidEvaluator
    {
      template <int, int>
      class EvaluatorInterface;
    }

    namespace PoroFluidManager
    {
      class PhaseManagerInterface;
      template <int, int>
      class VariableManagerInterface;
    }  // namespace PoroFluidManager

    /*!
    \brief implementation of evaluation routines of porous fluid multiphase  element

    This singleton class is responsible for evaluating boundary terms.
    It provides the method Evaluate(...) which performs the actual evaluation
    depending on the action provided by the global algorithm.


    \author vuong
    */
    template <Core::FE::CellType distype>
    class PoroFluidMultiPhaseEleCalc : public PoroFluidMultiPhaseEleInterface
    {
     protected:
      /// protected constructor, since we are a Singleton.
      /// this constructor is called from a derived class
      /// -> therefore, it has to be protected instead of private
      PoroFluidMultiPhaseEleCalc(const int numdofpernode, const std::string& disname);

     public:
      //! Singleton access method
      static PoroFluidMultiPhaseEleCalc<distype>* Instance(
          const int numdofpernode, const std::string& disname);

      /*========================================================================*/
      //! @name static member variables
      /*========================================================================*/

      //! number of element nodes (nomenclature: T. Hughes, The finite element method)
      static constexpr int nen_ = Core::FE::num_nodes<distype>;

      //! number of space dimensions
      static constexpr int nsd_ = Core::FE::dim<distype>;

      //! number of components necessary to store second derivatives
      // 1 component  for nsd=1:  (N,xx)
      // 3 components for nsd=2:  (N,xx ; N,yy ; N,xy)
      // 6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
      static constexpr int numderiv2_ = Core::FE::DisTypeToNumDeriv2<distype>::numderiv2;

      //! element-type specific flag if second derivatives are needed
      static constexpr bool use2ndderiv_ =
          POROFLUIDMULTIPHASE::ElementUtils::Use2ndDerivs<distype>::use;

      /// Evaluate the element
      /*!
        Generic virtual interface function. Called via base pointer.
       */
      int Evaluate(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la,
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,
          std::vector<Core::LinAlg::SerialDenseVector*>& elevec) override;

     protected:
      /*========================================================================*/
      //! @name general framework
      /*========================================================================*/

      /// Setup element evaluation
      virtual int setup_calc(Core::Elements::Element* ele, Discret::Discretization& discretization,
          const POROFLUIDMULTIPHASE::Action& action);

      //! evaluate action
      virtual int evaluate_action(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, const POROFLUIDMULTIPHASE::Action& action,
          Core::Elements::Element::LocationArray& la,
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,
          std::vector<Core::LinAlg::SerialDenseVector*>& elevec);

      //! extract element based or nodal values
      //  return extracted values of phinp
      virtual void extract_element_and_node_values(Core::Elements::Element* ele,
          Teuchos::ParameterList& params, Discret::Discretization& discretization,
          Core::Elements::Element::LocationArray& la);

      /// Setup element evaluation
      virtual void prepare_gauss_point_loop(
          Core::Elements::Element* ele  ///< the element whose matrix is calculated
      );

      void gauss_point_loop(
          const Core::FE::IntPointsAndWeights<nsd_>& intpoints,   ///< integration points
          Core::Elements::Element* ele,                           //!< current element
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,  //!< element matrixes to calculate
          std::vector<Core::LinAlg::SerialDenseVector*>&
              elevec,                                 //!< element rhs vectors to calculate
          Discret::Discretization& discretization,    //!< discretization
          Core::Elements::Element::LocationArray& la  //!< location array
      );

      void gauss_point_loop_od_struct(
          const Core::FE::IntPointsAndWeights<nsd_>& intpoints,   ///< integration points
          Core::Elements::Element* ele,                           //!< current element
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,  //!< element matrixes to calculate
          std::vector<Core::LinAlg::SerialDenseVector*>&
              elevec,                                 //!< element rhs vectors to calculate
          Discret::Discretization& discretization,    //!< discretization
          Core::Elements::Element::LocationArray& la  //!< location array
      );

      void gauss_point_loop_od_scatra(
          const Core::FE::IntPointsAndWeights<nsd_>& intpoints,   ///< integration points
          Core::Elements::Element* ele,                           //!< current element
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,  //!< element matrixes to calculate
          std::vector<Core::LinAlg::SerialDenseVector*>&
              elevec,                                 //!< element rhs vectors to calculate
          Discret::Discretization& discretization,    //!< discretization
          Core::Elements::Element::LocationArray& la  //!< location array
      );

      //! calculate matrix and rhs. Here the whole thing is hidden.
      void gauss_point_loop(Core::Elements::Element* ele,         //!< current element
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,  //!< element matrixes to calculate
          std::vector<Core::LinAlg::SerialDenseVector*>&
              elevec,                                 //!< element rhs vectors to calculate
          Discret::Discretization& discretization,    //!< discretization
          Core::Elements::Element::LocationArray& la  //!< location array
      );

      //! evaluate at all Gauss points and average
      void gauss_point_loop_average(Core::Elements::Element* ele,  //!< current element
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,  //!< element matrices to calculate
          std::vector<Core::LinAlg::SerialDenseVector*>&
              elevec,                                 //!< element rhs vectors to calculate
          Discret::Discretization& discretization,    //!< discretization
          Core::Elements::Element::LocationArray& la  //!< location array
      );

      //! calculate off-diagonal fluid-struct-coupling matrix. Here the whole thing is hidden.
      void gauss_point_loop_od_struct(Core::Elements::Element* ele,  //!< current element
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,  //!< element matrixes to calculate
          std::vector<Core::LinAlg::SerialDenseVector*>&
              elevec,                                 //!< element rhs vectors to calculate
          Discret::Discretization& discretization,    //!< discretization
          Core::Elements::Element::LocationArray& la  //!< location array
      );

      //! calculate off-diagonal fluid-scatra-coupling matrix. Here the whole thing is hidden.
      void gauss_point_loop_od_scatra(Core::Elements::Element* ele,  //!< current element
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,  //!< element matrixes to calculate
          std::vector<Core::LinAlg::SerialDenseVector*>&
              elevec,                                 //!< element rhs vectors to calculate
          Discret::Discretization& discretization,    //!< discretization
          Core::Elements::Element::LocationArray& la  //!< location array
      );

      //! evaluate shape functions and their derivatives at current integration point
      double eval_shape_func_and_derivs_at_int_point(
          const Core::FE::IntPointsAndWeights<nsd_>& intpoints,  //!< integration points
          const int iquad                                        //!< id of current Gauss point
      );

      //! evaluate shape functions and their derivatives at current integration point
      double eval_shape_func_and_derivs_in_parameter_space();

      // Compute Jacobian (determinant of deformation gradient) at node 'inode'
      void compute_jacobian_at_node(const int inode  //!< local node number
      );

      /*========================================================================*/
      //! @name routines for additional element evaluations (called from evaluate_action)
      /*========================================================================*/

      //! loop over nodes and evaluate element
      void node_loop(Core::Elements::Element* ele,                //!< current element
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,  //!< element matrixes to calculate
          std::vector<Core::LinAlg::SerialDenseVector*>&
              elevec,                                  //!< element rhs vectors to calculate
          Discret::Discretization& discretization,     //!< discretization
          Core::Elements::Element::LocationArray& la,  //!< location array
          const bool jacobian_needed                   //!< necessary to compute Jacobian at node
      );

      //! evaluate just the element
      void evaluate_only_element(Core::Elements::Element* ele,    //!< current element
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,  //!< element matrixes to calculate
          std::vector<Core::LinAlg::SerialDenseVector*>&
              elevec,                                 //!< element rhs vectors to calculate
          Discret::Discretization& discretization,    //!< discretization
          Core::Elements::Element::LocationArray& la  //!< location array
      );

     private:
      /*========================================================================*/
      //! @name can be very useful
      /*========================================================================*/

      //! element
      Core::Elements::Element* ele_;

      /*========================================================================*/
      //! @name dofs and nodes
      /*========================================================================*/

      //! number of dof per node (= number of fluid phases + number of volume fractions)
      const int totalnumdofpernode_;

      //! number of fluid phases
      int numfluidphases_;

      /*========================================================================*/
      //! @name parameter lists
      /*========================================================================*/

      //! pointer to general scalar transport parameter class
      Discret::ELEMENTS::PoroFluidMultiPhaseEleParameter* para_;

      /*========================================================================*/
      //! @name Galerkin approximation and related
      /*========================================================================*/

      //! coordinates of current integration point in reference coordinates
      Core::LinAlg::Matrix<nsd_, 1> xsi_;
      //! initial node coordinates
      Core::LinAlg::Matrix<nsd_, nen_> xyze0_;
      //! current node coordinates
      Core::LinAlg::Matrix<nsd_, nen_> xyze_;
      //! array for shape functions
      Core::LinAlg::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives w.r.t r,s,t
      Core::LinAlg::Matrix<nsd_, nen_> deriv_;
      //! array for second derivatives of shape function w.r.t r,s,t
      Core::LinAlg::Matrix<numderiv2_, nen_> deriv2_;
      //! global derivatives of shape functions w.r.t x,y,z
      Core::LinAlg::Matrix<nsd_, nen_> derxy_;
      //! global second derivatives of shape functions w.r.t x,y,z
      Core::LinAlg::Matrix<numderiv2_, nen_> derxy2_;

      //! transposed jacobian "dx/ds"
      Core::LinAlg::Matrix<nsd_, nsd_> xjm_;
      //! inverse of transposed jacobian "ds/dx"
      Core::LinAlg::Matrix<nsd_, nsd_> xij_;
      //! determinant of jacobian "dx/ds"
      double det_;
      //! determinant of deformation gradient "dx/dX"
      double j_;

      /*========================================================================*/
      //! @name scalar degrees of freedom and related
      /*========================================================================*/
      //! manager class for variables
      Teuchos::RCP<PoroFluidManager::VariableManagerInterface<nsd_, nen_>> variablemanager_;

      //! manager class for handling phases and corresponding DOFs
      Teuchos::RCP<PoroFluidManager::PhaseManagerInterface> phasemanager_;

      //! manager class for evaluation
      Teuchos::RCP<PoroFluidEvaluator::EvaluatorInterface<nsd_, nen_>> evaluator_;
    };

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
