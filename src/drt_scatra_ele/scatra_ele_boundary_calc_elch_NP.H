/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for Nernst-Planck formulation

\maintainer Christoph Schmidt

\level 2
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_NP_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_NP_H

#include "scatra_ele_boundary_calc_elch.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleBoundaryCalcElchNP : public ScaTraEleBoundaryCalcElch<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype> myelch;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcElchNP<distype>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname,
          const ScaTraEleBoundaryCalcElchNP* delete_me = NULL);

      //! destructor
      virtual ~ScaTraEleBoundaryCalcElchNP() { return; };

      //! singleton destruction
      virtual void Done();

     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcElchNP(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate action
      int EvaluateAction(DRT::FaceElement* ele,      //!< boundary element
          Teuchos::ParameterList& params,            //!< parameter list
          DRT::Discretization& discretization,       //!< discretization
          SCATRA::BoundaryAction action,             //!< action
          DRT::Element::LocationArray& la,           //!< location array
          Epetra_SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
          Epetra_SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
          Epetra_SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
          Epetra_SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
          Epetra_SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
      );

      //! evaluate Neumann boundary condition
      int EvaluateNeumann(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition,
          DRT::Element::LocationArray& la, Epetra_SerialDenseVector& elevec1, const double scalar);

      //! evaluate an electrode kinetics boundary condition
      void EvaluateElchBoundaryKinetics(const DRT::Element* ele,  ///< current element
          Epetra_SerialDenseMatrix& emat,                         ///< element matrix
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

      //! extract valence of species k from element material
      double GetValence(const Teuchos::RCP<const MAT::Material>& material,  //! element material
          const int k                                                       //! species number
          ) const;
    };  // class ScaTraEleBoundaryCalcElchNP
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
