/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_electrode_growth.H

\brief evaluation of ScaTra boundary elements for isothermal electrodes exhibiting surface layer
growth, e.g., lithium plating

\level 2

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
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

      //! evaluate scatra-scatra interface coupling condition
      void EvaluateS2ICoupling(const DRT::FaceElement* ele,  ///< current boundary element
          Teuchos::ParameterList& params,                    ///< parameter list
          DRT::Discretization& discretization,               ///< discretization
          DRT::Element::LocationArray& la,                   ///< location array
          Epetra_SerialDenseMatrix& eslavematrix,            ///< element matrix for slave side
          Epetra_SerialDenseMatrix& emastermatrix,           ///< element matrix for master side
          Epetra_SerialDenseVector& eslaveresidual           ///< element residual for slave side
      );

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

      //! compute Butler-Volmer current density via Newton-Raphson iteration
      double GetButlerVolmerCurrentDensity(const double& i0,  //!< exchange current density
          const double& alphaa,                               //!< anodic transfer coefficient
          const double& alphac,                               //!< cathodic transfer coefficient
          const double& frt,                                  //!< factor F/RT
          const double& pot_ed,                               //!< electrode-side electric potential
          const double& pot_el,      //!< electrolyte-side electric potential
          const double& epd,         //!< half-cell open-circuit potential
          const double& resistance,  //!< scatra-scatra interface layer resistance
          const double& thickness,   //!< scatra-scatra interface layer thickness
          const Teuchos::RCP<DRT::Condition>&
              condition  //!< scatra-scatra interface coupling condition
          ) const;

      //! compute regularization factor for lithium stripping
      double GetRegularizationFactor(
          const double& thickness,  //!< scatra-scatra interface layer thickness
          const Teuchos::RCP<DRT::Condition>&
              condition,    //!< scatra-scatra interface coupling condition
          const double eta  //!< electrode-electrolyte overpotential at integration point
          ) const;

      //! compute derivative of regularization factor for lithium stripping w.r.t. thickness of
      //! plated lithium
      double GetRegularizationFactorDerivative(
          const double& thickness,  //!< scatra-scatra interface layer thickness
          const Teuchos::RCP<DRT::Condition>&
              condition,    //!< scatra-scatra interface coupling condition
          const double eta  //!< electrode-electrolyte overpotential at integration point
          ) const;

      //! nodal growth variables
      LINALG::Matrix<my::nen_, 1> egrowth_;
    };  // class ScaTraEleBoundaryCalcElchElectrodeGrowth
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
