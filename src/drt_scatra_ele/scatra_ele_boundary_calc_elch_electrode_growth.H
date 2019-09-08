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
      );

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
      );

      /*!
       * \brief compute Butler-Volmer current density via Newton-Raphson iteration
       *
       * @param[in] i0         exchange current density
       * @param[in] alphaa     symmetry coefficient of anodic intercalation reaction
       * @param[in] alphac     symmetry coefficient of cathodic intercalation reaction
       * @param[in] frt        factor F/RT
       * @param[in] pot_ed     electrode-side electric potential
       * @param[in] pot_el     electrolyte-side electric potential
       * @param[in] epd        half-cell open-circuit potential
       * @param[in] resistance scatra-scatra interface layer resistance
       * @param[in] thickness  scatra-scatra interface layer thickness
       * @return Butler-Volmer current density
       */
      double GetButlerVolmerCurrentDensity(const double i0, const double alphaa,
          const double alphac, const double frt, const double pot_ed, const double pot_el,
          const double epd, const double resistance, const double thickness) const;

      /*!
       * \brief compute regularization factor for lithium stripping
       *
       * @param[in] thickness  scatra-scatra interface layer thickness
       * @param[in] eta        electrode-electrolyte overpotential at integration point
       * @return  return the regularization factor if regularization is applied
       */
      double GetRegularizationFactor(const double thickness, const double eta) const;

      /*!
       * \brief compute derivative of regularization factor for lithium stripping w.r.t. thickness
       * of plated lithium
       *
       * @param[in] thickness  scatra-scatra interface layer thickness
       * @param[in] eta        electrode-electrolyte overpotential at integration point
       * @return  return the derivative of the regularization factor w.r.t. thickness of deposited
       * material if regularization is applied
       */
      double GetRegularizationFactorDerivative(const double thickness, const double eta) const;

      //! nodal growth variables
      LINALG::Matrix<my::nen_, 1> egrowth_;
    };  // class ScaTraEleBoundaryCalcElchElectrodeGrowth
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
