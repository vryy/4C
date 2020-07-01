/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluate heat transport within binary, concentrated electrolytes on element level

\level 2

*/
/*--------------------------------------------------------------------------*/
#ifndef SCATRA_ELE_CALC_STI_DIFFCOND_H
#define SCATRA_ELE_CALC_STI_DIFFCOND_H

#include "scatra_ele_calc.H"
#include "scatra_ele_sti_elch.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleCalcElch;
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleCalcElchDiffCond;
    class ScaTraEleDiffManagerElchDiffCond;
    class ScaTraEleDiffManagerSTIThermo;
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchDiffCond;
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleUtilsElchDiffCond;

    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleCalcSTIDiffCond : public ScaTraEleCalc<distype>, public ScaTraEleSTIElch<distype>
    {
     public:
      //! singleton access method
      static ScaTraEleCalcSTIDiffCond<distype>* Instance(const int numdofpernode, const int numscal,
          const std::string& disname, const ScaTraEleCalcSTIDiffCond* delete_me = NULL);

      //! singleton destruction
      virtual void Done();

     private:
      //! abbreviations
      typedef ScaTraEleCalc<distype> my;
      typedef ScaTraEleCalcElch<distype> myelch;
      typedef ScaTraEleCalcElchDiffCond<distype> mydiffcond;
      typedef ScaTraEleSTIElch<distype> mystielch;

      //! private constructor for singletons
      ScaTraEleCalcSTIDiffCond(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! destructor
      virtual ~ScaTraEleCalcSTIDiffCond() { return; };

      //! evaluate action for off-diagonal system matrix block
      int EvaluateActionOD(DRT::Element* ele,        //!< current element
          Teuchos::ParameterList& params,            //!< parameter list
          DRT::Discretization& discretization,       //!< discretization
          const SCATRA::Action& action,              //!< action parameter
          DRT::Element::LocationArray& la,           //!< location array
          Epetra_SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
          Epetra_SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
          Epetra_SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
          Epetra_SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
          Epetra_SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
      );

      //! calculate element matrix and element right-hand side vector
      void Sysmat(DRT::Element* ele,           ///< current element
          Epetra_SerialDenseMatrix& emat,      ///< element matrix
          Epetra_SerialDenseVector& erhs,      ///< element right-hand side vector
          Epetra_SerialDenseVector& subgrdiff  ///< subgrid diffusivity scaling vector
      );

      //! element matrix and right-hand side vector contributions arising from Joule's heat
      void CalcMatAndRhsJoule(Epetra_SerialDenseMatrix& emat,  //!< element matrix
          Epetra_SerialDenseVector& erhs,                      //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
      );

      //! element matrix and right-hand side vector contributions arising from heat of mixing
      void CalcMatAndRhsMixing(Epetra_SerialDenseMatrix& emat,  //!< element matrix
          Epetra_SerialDenseVector& erhs,                       //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
      );

      //! element matrix and right-hand side vector contributions arising from Soret effect
      void CalcMatAndRhsSoret(Epetra_SerialDenseMatrix& emat,  //!< element matrix
          Epetra_SerialDenseVector& erhs,                      //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
      );

      //! fill element matrix with linearizations of discrete thermo residuals w.r.t. scatra dofs
      void SysmatODThermoScatra(DRT::Element* ele,  //!< current element
          Epetra_SerialDenseMatrix& emat            //!< element matrix
      );

      //! provide element matrix with linearizations of Joule's heat term in discrete thermo
      //! residuals w.r.t. scatra dofs
      void CalcMatJouleOD(Epetra_SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
      );

      //! provide element matrix with linearizations of heat of mixing term in discrete thermo
      //! residuals w.r.t. scatra dofs
      void CalcMatMixingOD(Epetra_SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
      );

      //! provide element matrix with linearizations of Soret effect term in discrete thermo
      //! residuals w.r.t. scatra dofs
      void CalcMatSoretOD(Epetra_SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
      );

      //! extract quantities for element evaluation
      void ExtractElementAndNodeValues(DRT::Element* ele,  //!< current element
          Teuchos::ParameterList& params,                  //!< parameter list
          DRT::Discretization& discretization,             //!< discretization
          DRT::Element::LocationArray& la                  //!< location array
      );

      //! get material parameters
      void GetMaterialParams(const DRT::Element* ele,  //!< current element
          std::vector<double>& densn,                  //!< density at t_(n)
          std::vector<double>& densnp,                 //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,                 //!< density at t_(n+alpha_M)
          double& visc,                                //!< fluid viscosity
          const int iquad = -1                         //!< ID of current integration point
      );

      //! evaluate Soret material
      void MatSoret(const Teuchos::RCP<const MAT::Material> material,  //!< Soret material
          double& densn,                                               //!< density at time t_(n)
          double& densnp,  //!< density at time t_(n+1) or t_(n+alpha_F)
          double& densam   //!< density at time t_(n+alpha_M)
      );

      //! set internal variables for element evaluation
      void SetInternalVariablesForMatAndRHS();

      //! get thermo diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerSTIThermo> DiffManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerSTIThermo>(my::diffmanager_);
      };

      //! get internal variable manager for heat transfer within electrochemical substances
      Teuchos::RCP<ScaTraEleInternalVariableManagerSTIElch<my::nsd_, my::nen_>> VarManager()
      {
        return Teuchos::rcp_static_cast<
            ScaTraEleInternalVariableManagerSTIElch<my::nsd_, my::nen_>>(my::scatravarmanager_);
      };

      //! diffusion manager for diffusion-conduction formulation
      Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> diffmanagerdiffcond_;

      //! utility class supporting element evaluation for diffusion-conduction formulation
      DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>* utils_;
    };  // class ScaTraEleCalcSTIDiffCond
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
