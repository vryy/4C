/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluation of a scatra element that does not contain any physics

\level 2


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_SCATRA_ELE_CALC_NO_PHYSICS_HPP
#define FOUR_C_SCATRA_ELE_CALC_NO_PHYSICS_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    /**
     * \brief ScaTra ImplType containing no physics.
     *
     * Evaluation of a scatra element that does not contain any physics. Currently, this class only
     * implements the minimal set of actions needed for reading the scatra results from a restart
     * file and simulating a one-way coupling to the structure. This ImplType can currently not be
     * used in solving the scatra equations, as the needed actions are not implemented yet.
     *
     * @tparam distype Element shape type
     * @tparam probdim Number of space dimensions of the problem
     */
    template <CORE::FE::CellType distype, int probdim>
    class ScaTraEleCalcNoPhysics : public ScaTraEleCalc<distype, probdim>
    {
     public:
      //! abbreviation
      typedef ScaTraEleCalc<distype, probdim> my;

      //! singleton access method
      static ScaTraEleCalcNoPhysics<distype, probdim>* Instance(
          int numdofpernode, int numscal, const std::string& disname);

      //! evaluate the element
      int evaluate_action(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

     protected:
      //! protected constructor for singletons
      ScaTraEleCalcNoPhysics(int numdofpernode, int numscal, const std::string& disname);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
