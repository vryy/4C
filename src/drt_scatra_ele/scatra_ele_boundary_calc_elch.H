/*----------------------------------------------------------------------*/
/*!

\brief evaluation of scatra boundary terms at integration points

\level 2

\maintainer Christoph Schmidt
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_H

#include "scatra_ele_boundary_calc.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleParameterElch;
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleUtilsElch;

    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleBoundaryCalcElch : public ScaTraEleBoundaryCalc<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;

     public:
      //! singleton access method
      // not needed, since class is purely virtual

      //! destructor
      virtual ~ScaTraEleBoundaryCalcElch() { return; };

      //! singleton destruction
      virtual void Done() = 0;

     protected:
      //! protected constructor for singletons
      ScaTraEleBoundaryCalcElch(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate action
      virtual int EvaluateAction(DRT::FaceElement* ele,  //!< boundary element
          Teuchos::ParameterList& params,                //!< parameter list
          DRT::Discretization& discretization,           //!< discretization
          SCATRA::BoundaryAction action,                 //!< action
          DRT::Element::LocationArray& la,               //!< location array
          Epetra_SerialDenseMatrix& elemat1_epetra,      //!< element matrix 1
          Epetra_SerialDenseMatrix& elemat2_epetra,      //!< element matrix 2
          Epetra_SerialDenseVector& elevec1_epetra,      //!< element right-hand side vector 1
          Epetra_SerialDenseVector& elevec2_epetra,      //!< element right-hand side vector 2
          Epetra_SerialDenseVector& elevec3_epetra       //!< element right-hand side vector 3
      );

      //! evaluate an electrode kinetics boundary condition
      virtual void EvaluateElchBoundaryKinetics(const DRT::Element* ele,  ///< current element
          Epetra_SerialDenseMatrix& emat,                                 ///< element matrix
          Epetra_SerialDenseVector& erhs,  ///< element right-hand side vector
          const std::vector<LINALG::Matrix<my::nen_, 1>>&
              ephinp,  ///< nodal values of concentration and electric potential
          const std::vector<LINALG::Matrix<my::nen_, 1>>& ehist,  ///< nodal history vector
          double timefac,                                         ///< time factor
          Teuchos::RCP<const MAT::Material> material,             ///< material
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
          Epetra_SerialDenseMatrix& elemat1_epetra,         ///< element matrix
          Epetra_SerialDenseVector& elevec1_epetra,         ///< element right-hand side vector
          const double
              scalar  ///< scaling factor for element matrix and right-hand side contributions
      );

      //! evaluate electrode kinetics status information
      void EvaluateElectrodeStatus(const DRT::Element* ele,  ///< current element
          Epetra_SerialDenseVector& scalars,                 ///< scalars to be integrated
          Teuchos::ParameterList& params,                    ///< parameter list
          Teuchos::RCP<DRT::Condition> cond,                 ///< condition
          const std::vector<LINALG::Matrix<my::nen_, 1>>&
              ephinp,  ///< nodal values of concentration and electric potential
          const std::vector<LINALG::Matrix<my::nen_, 1>>&
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
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseVector& elevec1_epetra);

      //! calculate cell voltage
      void CalcCellVoltage(const DRT::Element* ele,  //!< the element we are dealing with
          Teuchos::ParameterList& params,            //!< parameter list
          DRT::Discretization& discretization,       //!< discretization
          DRT::Element::LocationArray& la,           //!< location array
          Epetra_SerialDenseVector& scalars  //!< result vector for scalar integrals to be computed
      );

      //! extract valence of species k from element material
      virtual double GetValence(
          const Teuchos::RCP<const MAT::Material>& material,  //! element material
          const int k                                         //! species number
          ) const = 0;

      //! parameter class for electrochemistry problems
      const DRT::ELEMENTS::ScaTraEleParameterElch* elchparams_;

      //! utility class supporting element evaluation
      const DRT::ELEMENTS::ScaTraEleUtilsElch<distype>* utils_;
    };  // class ScaTraEleBoundaryCalcElch
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
