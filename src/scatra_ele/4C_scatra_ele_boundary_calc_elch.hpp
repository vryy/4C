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

namespace Discret
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleParameterElch;
    template <Core::FE::CellType distype>
    class ScaTraEleUtilsElch;

    // class implementation
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcElch : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      typedef Discret::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim> my;

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
      int evaluate_action(Core::Elements::FaceElement* ele,  //!< boundary element
          Teuchos::ParameterList& params,                    //!< parameter list
          Core::FE::Discretization& discretization,          //!< discretization
          ScaTra::BoundaryAction action,                     //!< action
          Core::Elements::Element::LocationArray& la,        //!< location array
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,   //!< element matrix 1
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,   //!< element matrix 2
          Core::LinAlg::SerialDenseVector& elevec1_epetra,   //!< element right-hand side vector 1
          Core::LinAlg::SerialDenseVector& elevec2_epetra,   //!< element right-hand side vector 2
          Core::LinAlg::SerialDenseVector& elevec3_epetra    //!< element right-hand side vector 3
          ) override;

      //! evaluate an electrode kinetics boundary condition
      virtual void evaluate_elch_boundary_kinetics(
          const Core::Elements::Element* ele,     ///< current element
          Core::LinAlg::SerialDenseMatrix& emat,  ///< element matrix
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
      );

      //! process an electrode kinetics boundary condition
      void calc_elch_boundary_kinetics(Core::Elements::FaceElement* ele,  ///< current element
          Teuchos::ParameterList& params,                                 ///< parameter list
          Core::FE::Discretization& discretization,                       ///< discretization
          Core::Elements::Element::LocationArray& la,                     ///< location array
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,                ///< element matrix
          Core::LinAlg::SerialDenseVector& elevec1_epetra,  ///< element right-hand side vector
          const double
              scalar  ///< scaling factor for element matrix and right-hand side contributions
      );

      //! evaluate electrode kinetics status information
      void evaluate_electrode_status(const Core::Elements::Element* ele,  ///< current element
          Core::LinAlg::SerialDenseVector& scalars,        ///< scalars to be integrated
          Teuchos::ParameterList& params,                  ///< parameter list
          Teuchos::RCP<Core::Conditions::Condition> cond,  ///< condition
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>&
              ephinp,  ///< nodal values of concentration and electric potential
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>&
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
      void calc_nernst_linearization(Core::Elements::FaceElement* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra);

      //! calculate cell voltage
      void calc_cell_voltage(
          const Core::Elements::Element* ele,          //!< the element we are dealing with
          Teuchos::ParameterList& params,              //!< parameter list
          Core::FE::Discretization& discretization,    //!< discretization
          Core::Elements::Element::LocationArray& la,  //!< location array
          Core::LinAlg::SerialDenseVector&
              scalars  //!< result vector for scalar integrals to be computed
      );

      //! extract valence of species k from element material
      virtual double get_valence(
          const Teuchos::RCP<const Core::Mat::Material>& material,  //! element material
          const int k                                               //! species number
      ) const = 0;

      //! parameter class for electrochemistry problems
      const Discret::ELEMENTS::ScaTraEleParameterElch* elchparams_;

      //! utility class supporting element evaluation
      const Discret::ELEMENTS::ScaTraEleUtilsElch<distype>* utils_;
    };  // class ScaTraEleBoundaryCalcElch
  }     // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
