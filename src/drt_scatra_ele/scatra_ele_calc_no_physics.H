/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluation of a scatra element that does not contain any physics

\level 2


*/
/*----------------------------------------------------------------------*/


#ifndef BACI_SCATRA_ELE_CALC_NO_PHYSICS_H
#define BACI_SCATRA_ELE_CALC_NO_PHYSICS_H

#include "scatra_ele_calc.H"

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
    template <DRT::Element::DiscretizationType distype, int probdim>
    class ScaTraEleCalcNoPhysics : public ScaTraEleCalc<distype, probdim>
    {
     public:
      //! abbreviation
      typedef ScaTraEleCalc<distype, probdim> my;

      //! destructor
      virtual ~ScaTraEleCalcNoPhysics() = default;

      //! singleton access method
      static ScaTraEleCalcNoPhysics<distype, probdim>* Instance(int numdofpernode, int numscal,
          const std::string& disname, const ScaTraEleCalcNoPhysics* delete_me = NULL);

      //! called upon destruction
      virtual void Done();

      //! evaluate the element
      int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra);

     protected:
      //! protected constructor for singletons
      ScaTraEleCalcNoPhysics(int numdofpernode, int numscal, const std::string& disname);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // BACI_SCATRA_ELE_CALC_NO_PHYSICS_H
