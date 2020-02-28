/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for isothermal electrodes exhibiting surface layer
growth, e.g., lithium plating

\level 2

\maintainer Christoph Schmidt
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_GROWTH_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_GROWTH_H

#include "scatra_ele_boundary_calc_elch_electrode.H"

// forward declaration
namespace DRT
{
  namespace ELEMENTS
  {
    class ScaTraEleBoundaryCalcElchElectrodeGrowthUtils;
  }
}  // namespace DRT

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleBoundaryCalcElchElectrodeGrowth
        : public ScaTraEleBoundaryCalcElchElectrode<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype> myelch;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype> myelectrode;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowthUtils myelectrodegrowthutils;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname,
          const ScaTraEleBoundaryCalcElchElectrodeGrowth* delete_me = NULL);

      //! destructor
      virtual ~ScaTraEleBoundaryCalcElchElectrodeGrowth() { return; };

      //! singleton destruction
      virtual void Done();

     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcElchElectrodeGrowth(
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
          ) override;

      //! evaluate minimum and maximum interfacial overpotential associated with scatra-scatra
      //! interface layer growth
      void EvaluateMinMaxOverpotential(const DRT::FaceElement* ele,  //!< current boundary element
          Teuchos::ParameterList& params,                            //!< parameter list
          DRT::Discretization& discretization,                       //!< discretization
          DRT::Element::LocationArray& la                            //!< location array
      );

      /*!
       * \brief evaluate scatra-scatra interface coupling condition
       *
       * @param[in] ele              current boundary element
       * @param[in] params           parameter list
       * @param[in] discretization   discretization
       * @param[in] la               location array
       * @param[out] eslavematrix    element matrix for slave side
       * @param[out] emastermatrix   element matrix for master side
       * @param[out] eslaveresidual  element residual for slave side
       */
      void EvaluateS2ICoupling(const DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& eslavematrix, Epetra_SerialDenseMatrix& emastermatrix,
          Epetra_SerialDenseVector& eslaveresidual) override;

      //! evaluate global growth-growth matrix block for scatra-scatra interface coupling involving
      //! interface layer growth
      void EvaluateS2ICouplingGrowthGrowth(
          const DRT::FaceElement* ele,              ///< current boundary element
          Teuchos::ParameterList& params,           ///< parameter list
          DRT::Discretization& discretization,      ///< discretization
          DRT::Element::LocationArray& la,          ///< location array
          Epetra_SerialDenseMatrix& eslavematrix,   ///< element matrix for slave side
          Epetra_SerialDenseVector& eslaveresidual  ///< element residual for slave side
      );

      //! evaluate global growth-scatra matrix block for scatra-scatra interface coupling involving
      //! interface layer growth
      void EvaluateS2ICouplingGrowthScatra(
          const DRT::FaceElement* ele,             ///< current boundary element
          Teuchos::ParameterList& params,          ///< parameter list
          DRT::Discretization& discretization,     ///< discretization
          DRT::Element::LocationArray& la,         ///< location array
          Epetra_SerialDenseMatrix& eslavematrix,  ///< element matrix for slave side
          Epetra_SerialDenseMatrix& emastermatrix  ///< element matrix for master side
      );

      //! evaluate global scatra-growth matrix block for scatra-scatra interface coupling involving
      //! interface layer growth
      void EvaluateS2ICouplingScatraGrowth(
          const DRT::FaceElement* ele,            ///< current boundary element
          Teuchos::ParameterList& params,         ///< parameter list
          DRT::Discretization& discretization,    ///< discretization
          DRT::Element::LocationArray& la,        ///< location array
          Epetra_SerialDenseMatrix& eslavematrix  ///< element matrix for slave side
      );

      //! extract nodal state variables associated with boundary element
      void ExtractNodeValues(const DRT::Discretization& discretization,  //!< discretization
          DRT::Element::LocationArray& la                                //!< location array
          ) override;

      //! nodal growth variables
      LINALG::Matrix<my::nen_, 1> egrowth_;
    };  // class ScaTraEleBoundaryCalcElchElectrodeGrowth
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
