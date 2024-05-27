/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluate heat transport within electrodes on element level

\level 2

*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_CALC_STI_ELECTRODE_HPP
#define FOUR_C_SCATRA_ELE_CALC_STI_ELECTRODE_HPP

#include "4C_config.hpp"

#include "4C_inpar_elch.hpp"
#include "4C_mat_electrode.hpp"
#include "4C_scatra_ele_calc.hpp"
#include "4C_scatra_ele_calc_elch_electrode.hpp"
#include "4C_scatra_ele_sti_elch.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleDiffManagerElchElectrode;
    class ScaTraEleDiffManagerSTIElchElectrode;
    class ScaTraEleDiffManagerSTIThermo;
    template <CORE::FE::CellType distype>
    class ScaTraEleUtilsElchElectrode;

    // class implementation
    template <CORE::FE::CellType distype>
    class ScaTraEleCalcSTIElectrode : public ScaTraEleCalc<distype>,
                                      public ScaTraEleSTIElch<distype>
    {
     public:
      //! singleton access method
      static ScaTraEleCalcSTIElectrode<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);



     private:
      //! abbreviations
      using my = ScaTraEleCalc<distype>;
      using mystielch = ScaTraEleSTIElch<distype>;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

      //! private constructor for singletons
      ScaTraEleCalcSTIElectrode(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate action for off-diagonal system matrix block
      int EvaluateActionOD(DRT::Element* ele,               //!< current element
          Teuchos::ParameterList& params,                   //!< parameter list
          DRT::Discretization& discretization,              //!< discretization
          const SCATRA::Action& action,                     //!< action parameter
          DRT::Element::LocationArray& la,                  //!< location array
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
          CORE::LINALG::SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
          CORE::LINALG::SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
          CORE::LINALG::SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
          ) override;

      //! calculate element matrix and element right-hand side vector
      void sysmat(DRT::Element* ele,                  ///< current element
          CORE::LINALG::SerialDenseMatrix& emat,      ///< element matrix
          CORE::LINALG::SerialDenseVector& erhs,      ///< element right-hand side vector
          CORE::LINALG::SerialDenseVector& subgrdiff  ///< subgrid diffusivity scaling vector
          ) override;

      //! fill element matrix with linearizations of discrete thermo residuals w.r.t. scatra dofs
      void sysmat_od_thermo_scatra(DRT::Element* ele,  //!< current element
          CORE::LINALG::SerialDenseMatrix& emat        //!< element matrix
      );

      //! element matrix and right-hand side vector contributions arising from Joule's heat
      void calc_mat_and_rhs_joule(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
          ) override;

      //! element matrix and right-hand side vector contributions arising from heat of mixing
      void calc_mat_and_rhs_mixing(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
          ) override;

      //! element matrix and right-hand side vector contributions arising from Soret effect
      void calc_mat_and_rhs_soret(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
          ) override;

      //! provide element matrix with linearizations of Joule's heat term in discrete thermo
      //! residuals w.r.t. scatra dofs
      void calc_mat_joule_od(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
          ) override;

      //! provide element matrix with linearizations of heat of mixing term in discrete thermo
      //! residuals w.r.t. scatra dofs
      void calc_mat_mixing_od(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
          ) override;

      //! provide element matrix with linearizations of Soret effect term in discrete thermo
      //! residuals w.r.t. scatra dofs
      void calc_mat_soret_od(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
          ) override;

      //! extract quantities for element evaluation
      void extract_element_and_node_values(DRT::Element* ele,  //!< current element
          Teuchos::ParameterList& params,                      //!< parameter list
          DRT::Discretization& discretization,                 //!< discretization
          DRT::Element::LocationArray& la                      //!< location array
          ) override;

      //! get material parameters
      void get_material_params(const DRT::Element* ele,  //!< current element
          std::vector<double>& densn,                    //!< density at t_(n)
          std::vector<double>& densnp,                   //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,                   //!< density at t_(n+alpha_M)
          double& visc,                                  //!< fluid viscosity
          const int iquad                                //!< ID of current integration point
          ) override;

      //! evaluate Soret material
      void mat_soret(const Teuchos::RCP<const CORE::MAT::Material> material,  //!< Soret material
          double& densn,   //!< density at time t_(n)
          double& densnp,  //!< density at time t_(n+1) or t_(n+alpha_F)
          double& densam   //!< density at time t_(n+alpha_M)
      );

      void mat_fourier(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< Fourier material
          double& densn,                                           //!< density at time t_(n)
          double& densnp,  //!< density at time t_(n+1) or t_(n+alpha_F)
          double& densam   //!< density at time t_(n+alpha_M)
      );

      //! set internal variables for element evaluation
      void set_internal_variables_for_mat_and_rhs() override;

      //! get thermo diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerSTIThermo> diff_manager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerSTIThermo>(my::diffmanager_);
      };

      //! get internal variable manager for heat transfer within electrochemical substances
      Teuchos::RCP<ScaTraEleInternalVariableManagerSTIElch<nsd_, nen_>> var_manager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerSTIElch<nsd_, nen_>>(
            my::scatravarmanager_);
      };

      //! diffusion manager for thermodynamic electrodes
      Teuchos::RCP<ScaTraEleDiffManagerSTIElchElectrode> diffmanagerstielectrode_;

      //! utility class supporting element evaluation for electrodes
      DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>* utils_;
    };  // class ScaTraEleCalcSTIElectrode


    //! implementation of ScaTraEleDiffManagerSTIElchElectrode
    class ScaTraEleDiffManagerSTIElchElectrode : public ScaTraEleDiffManagerElchElectrode
    {
     public:
      //! constructor
      ScaTraEleDiffManagerSTIElchElectrode(int numscal)
          :  // constructor of base class
            ScaTraEleDiffManagerElchElectrode(numscal),

            // initialize internal member variables
            ocp_(0.),
            ocpderiv_(0.),
            ocpderiv2_(0.){};


      //! compute and store half cell open circuit potential and its first and second derivatives
      //! w.r.t. concentration
      void SetOCPAndDerivs(
          const DRT::Element* ele, const double& concentration, const double& temperature)
      {
        const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
        const double gasconstant =
            DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
        // factor F/RT
        const double frt = faraday / (gasconstant * temperature);

        // access electrode material
        const Teuchos::RCP<const MAT::Electrode> matelectrode =
            Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->Material(1));
        if (matelectrode == Teuchos::null) FOUR_C_THROW("Invalid electrode material!");

        // no deformation available in this code part
        const double dummy_detF(1.0);
        // evaluate material
        ocp_ =
            matelectrode->compute_open_circuit_potential(concentration, faraday, frt, dummy_detF);
        ocpderiv_ = matelectrode->compute_d_open_circuit_potential_d_concentration(
            concentration, faraday, frt, dummy_detF);
        ocpderiv2_ =
            matelectrode->compute_d2_open_circuit_potential_d_concentration_d_concentration(
                concentration, faraday, frt, dummy_detF);
      };

      //! return half cell open circuit potential
      const double& GetOCP() { return ocp_; };

      //! return first derivative of half cell open circuit potential w.r.t. concentration
      const double& GetOCPDeriv() { return ocpderiv_; };

      //! return second derivative of half cell open circuit potential w.r.t. concentration
      const double& GetOCPDeriv2() { return ocpderiv2_; };

     protected:
      //! half cell open circuit potential
      double ocp_;

      //! first derivative of half cell open circuit potential w.r.t. concentration
      double ocpderiv_;

      //! second derivative of half cell open circuit potential w.r.t. concentration
      double ocpderiv2_;
    };  // class ScaTraEleDiffManagerSTIElchElectrode
  }     // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
