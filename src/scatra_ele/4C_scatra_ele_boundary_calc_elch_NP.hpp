/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for Nernst-Planck formulation


\level 2
 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_ELCH_NP_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_ELCH_NP_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_boundary_calc_elch.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    // class implementation
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcElchNP : public ScaTraEleBoundaryCalcElch<distype, probdim>
    {
      typedef Discret::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim> my;
      typedef Discret::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim> myelch;
      using my::nen_;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcElchNP<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);


     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcElchNP(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate action
      int evaluate_action(Core::Elements::FaceElement* ele,  //!< boundary element
          Teuchos::ParameterList& params,                    //!< parameter list
          Discret::Discretization& discretization,           //!< discretization
          ScaTra::BoundaryAction action,                     //!< action
          Core::Elements::Element::LocationArray& la,        //!< location array
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,   //!< element matrix 1
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,   //!< element matrix 2
          Core::LinAlg::SerialDenseVector& elevec1_epetra,   //!< element right-hand side vector 1
          Core::LinAlg::SerialDenseVector& elevec2_epetra,   //!< element right-hand side vector 2
          Core::LinAlg::SerialDenseVector& elevec3_epetra    //!< element right-hand side vector 3
          ) override;

      //! evaluate Neumann boundary condition
      int evaluate_neumann(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, Core::Conditions::Condition& condition,
          Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseVector& elevec1,
          const double scalar) override;

      //! evaluate an electrode kinetics boundary condition
      void evaluate_elch_boundary_kinetics(const Core::Elements::Element* ele,  ///< current element
          Core::LinAlg::SerialDenseMatrix& emat,                                ///< element matrix
          Core::LinAlg::SerialDenseVector& erhs,  ///< element right-hand side vector
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>&
              ephinp,  ///< nodal values of concentration and electric potential
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ehist,  ///< nodal history vector
          double timefac,                                           ///< time factor
          Teuchos::RCP<const Core::Mat::Material> material,         ///< material
          Teuchos::RCP<Core::Conditions::Condition>
              cond,                       ///< electrode kinetics boundary condition
          const int nume,                 ///< number of transferred electrons
          const std::vector<int> stoich,  ///< stoichiometry of the reaction
          const int kinetics,             ///< desired electrode kinetics model
          const double pot0,              ///< electrode potential on metal side
          const double frt,               ///< factor F/RT
          const double
              scalar  ///< scaling factor for element matrix and right-hand side contributions
          ) override;

      //! extract valence of species k from element material
      double get_valence(
          const Teuchos::RCP<const Core::Mat::Material>& material,  //! element material
          const int k                                               //! species number
      ) const override;
    };  // class ScaTraEleBoundaryCalcElchNP
  }     // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
