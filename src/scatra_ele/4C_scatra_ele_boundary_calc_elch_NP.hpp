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

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcElchNP : public ScaTraEleBoundaryCalcElch<distype, probdim>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim> my;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim> myelch;
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
      int EvaluateAction(DRT::FaceElement* ele,             //!< boundary element
          Teuchos::ParameterList& params,                   //!< parameter list
          DRT::Discretization& discretization,              //!< discretization
          SCATRA::BoundaryAction action,                    //!< action
          DRT::Element::LocationArray& la,                  //!< location array
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
          CORE::LINALG::SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
          CORE::LINALG::SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
          CORE::LINALG::SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
          ) override;

      //! evaluate Neumann boundary condition
      int evaluate_neumann(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseVector& elevec1,
          const double scalar) override;

      //! evaluate an electrode kinetics boundary condition
      void evaluate_elch_boundary_kinetics(const DRT::Element* ele,  ///< current element
          CORE::LINALG::SerialDenseMatrix& emat,                     ///< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  ///< element right-hand side vector
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephinp,  ///< nodal values of concentration and electric potential
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ehist,  ///< nodal history vector
          double timefac,                                           ///< time factor
          Teuchos::RCP<const CORE::MAT::Material> material,         ///< material
          Teuchos::RCP<CORE::Conditions::Condition>
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
      double GetValence(
          const Teuchos::RCP<const CORE::MAT::Material>& material,  //! element material
          const int k                                               //! species number
      ) const override;
    };  // class ScaTraEleBoundaryCalcElchNP
  }     // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
