/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_ELCH_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_ELCH_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_boundary_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleParameterElch;
    template <CORE::FE::CellType distype>
    class ScaTraEleUtilsElch;

    // class implementation
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcElch : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim> my;

     protected:
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      //! singleton access method
      // not needed, since class is purely virtual



     protected:
      //! protected constructor for singletons
      ScaTraEleBoundaryCalcElch(
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

      //! evaluate an electrode kinetics boundary condition
      virtual void EvaluateElchBoundaryKinetics(const DRT::Element* ele,  ///< current element
          CORE::LINALG::SerialDenseMatrix& emat,                          ///< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  ///< element right-hand side vector
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephinp,  ///< nodal values of concentration and electric potential
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ehist,  ///< nodal history vector
          double timefac,                                           ///< time factor
          Teuchos::RCP<const CORE::MAT::Material> material,         ///< material
          Teuchos::RCP<DRT::Condition> cond,  ///< electrode kinetics boundary condition
          const int nume,                     ///< number of transferred electrons
          const std::vector<int> stoich,      ///< stoichiometry of the reaction
          const int kinetics,                 ///< desired electrode kinetics model
          const double pot0,                  ///< electrode potential on metal side
          const double frt,                   ///< factor F/RT
          const double
              scalar  ///< scaling factor for element matrix and right-hand side contributions
      );

      //! process an electrode kinetics boundary condition
      void CalcElchBoundaryKinetics(DRT::FaceElement* ele,  ///< current element
          Teuchos::ParameterList& params,                   ///< parameter list
          DRT::Discretization& discretization,              ///< discretization
          DRT::Element::LocationArray& la,                  ///< location array
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,  ///< element matrix
          CORE::LINALG::SerialDenseVector& elevec1_epetra,  ///< element right-hand side vector
          const double
              scalar  ///< scaling factor for element matrix and right-hand side contributions
      );

      //! evaluate electrode kinetics status information
      void EvaluateElectrodeStatus(const DRT::Element* ele,  ///< current element
          CORE::LINALG::SerialDenseVector& scalars,          ///< scalars to be integrated
          Teuchos::ParameterList& params,                    ///< parameter list
          Teuchos::RCP<DRT::Condition> cond,                 ///< condition
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephinp,  ///< nodal values of concentration and electric potential
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephidtnp,                   ///< nodal time derivative vector
          const int kinetics,             ///< desired electrode kinetics model
          const std::vector<int> stoich,  ///< stoichiometry of the reaction
          const int nume,                 ///< number of transferred electrons
          const double pot0,              ///< electrode potential on metal side
          const double frt,               ///< factor F/RT
          const double timefac,           ///< time factor
          const double scalar             ///< scaling factor for current related quantities
      );

      //! evaluate linearization of nernst equation
      void CalcNernstLinearization(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra);

      //! calculate cell voltage
      void CalcCellVoltage(const DRT::Element* ele,  //!< the element we are dealing with
          Teuchos::ParameterList& params,            //!< parameter list
          DRT::Discretization& discretization,       //!< discretization
          DRT::Element::LocationArray& la,           //!< location array
          CORE::LINALG::SerialDenseVector&
              scalars  //!< result vector for scalar integrals to be computed
      );

      //! extract valence of species k from element material
      virtual double GetValence(
          const Teuchos::RCP<const CORE::MAT::Material>& material,  //! element material
          const int k                                               //! species number
      ) const = 0;

      //! parameter class for electrochemistry problems
      const DRT::ELEMENTS::ScaTraEleParameterElch* elchparams_;

      //! utility class supporting element evaluation
      const DRT::ELEMENTS::ScaTraEleUtilsElch<distype>* utils_;
    };  // class ScaTraEleBoundaryCalcElch
  }     // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
