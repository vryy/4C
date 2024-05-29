/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for conservation of mass concentration and electronic charge
within thermodynamic electrodes

\level 2

*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_CALC_ELCH_ELECTRODE_STI_THERMO_HPP
#define FOUR_C_SCATRA_ELE_CALC_ELCH_ELECTRODE_STI_THERMO_HPP

#include "4C_config.hpp"

#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_scatra_ele_calc_elch_electrode.hpp"
#include "4C_scatra_ele_calc_sti_electrode.hpp"
#include "4C_scatra_ele_sti_thermo.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchElectrodeSTIThermo;

    // class implementation
    template <CORE::FE::CellType distype>
    class ScaTraEleCalcElchElectrodeSTIThermo : public ScaTraEleCalcElchElectrode<distype>,
                                                public ScaTraEleSTIThermo<distype>
    {
     public:
      //! singleton access method
      static ScaTraEleCalcElchElectrodeSTIThermo<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     private:
      //! abbreviations
      using my = ScaTraEleCalc<distype>;
      using myelch = ScaTraEleCalcElch<distype>;
      using myelectrode = ScaTraEleCalcElchElectrode<distype>;
      using mythermo = ScaTraEleSTIThermo<distype>;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

      //! private constructor for singletons
      ScaTraEleCalcElchElectrodeSTIThermo(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate action for off-diagonal system matrix block
      int EvaluateActionOD(CORE::Elements::Element* ele,    //!< current element
          Teuchos::ParameterList& params,                   //!< parameter list
          DRT::Discretization& discretization,              //!< discretization
          const SCATRA::Action& action,                     //!< action parameter
          CORE::Elements::Element::LocationArray& la,       //!< location array
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
          CORE::LINALG::SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
          CORE::LINALG::SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
          CORE::LINALG::SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
          ) override;

      //! extract quantities for element evaluation
      void extract_element_and_node_values(CORE::Elements::Element* ele,  //!< current element
          Teuchos::ParameterList& params,                                 //!< parameter list
          DRT::Discretization& discretization,                            //!< discretization
          CORE::Elements::Element::LocationArray& la                      //!< location array
          ) override;

      //! get material parameters
      void get_material_params(const CORE::Elements::Element* ele,  //!< current element
          std::vector<double>& densn,                               //!< density at t_(n)
          std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,  //!< density at t_(n+alpha_M)
          double& visc,                 //!< fluid viscosity
          const int iquad               //!< ID of current integration point
          ) override;

      //! calculate element matrix and element right-hand side vector
      void calc_mat_and_rhs(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to calculate
          CORE::LINALG::SerialDenseVector& erhs,                    //!< element rhs to calculate+
          const int k,                                              //!< index of current scalar
          const double fac,                                         //!< domain-integration factor
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double taufac,  //!< tau times domain-integration factor
          const double
              timetaufac,  //!< domain-integration factor times tau times time-integration factor
          const double rhstaufac,  //!< time-integration factor for rhs times tau times
                                   //!< domain-integration factor
          CORE::LINALG::Matrix<nen_, 1>&
              tauderpot,  //!< derivatives of stabilization parameter w.r.t. electric potential
          double& rhsint  //!< rhs at Gauss point
          ) override;

      //! fill element matrix with linearizations of discrete scatra residuals w.r.t. thermo dofs
      void sysmat_od_scatra_thermo(CORE::Elements::Element* ele,  //!< current element
          CORE::LINALG::SerialDenseMatrix& emat                   //!< element matrix
      );

      //! set internal variables for element evaluation
      void set_internal_variables_for_mat_and_rhs() override;

      //! get internal variable manager for thermodynamic electrodes
      Teuchos::RCP<ScaTraEleInternalVariableManagerElchElectrodeSTIThermo<nsd_, nen_>> var_manager()
      {
        return Teuchos::rcp_static_cast<
            ScaTraEleInternalVariableManagerElchElectrodeSTIThermo<nsd_, nen_>>(
            my::scatravarmanager_);
      };

     private:
      // material type for evaluation
      CORE::Materials::MaterialType materialtype_;

    };  // class ScaTraEleCalcElchElectrodeSTIThermo


    //! implementation of ScaTraEleInternalVariableManagerElchElectrodeSTIThermo
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchElectrodeSTIThermo
        : public ScaTraEleInternalVariableManagerElchElectrode<NSD, NEN>,
          public ScaTraEleInternalVariableManagerSTIThermo<NSD, NEN>
    {
     public:
      //! abbreviations
      using vmelch = ScaTraEleInternalVariableManagerElch<NSD, NEN>;
      using vmelchelectrode = ScaTraEleInternalVariableManagerElchElectrode<NSD, NEN>;
      using vmthermo = ScaTraEleInternalVariableManagerSTIThermo<NSD, NEN>;

      //! constructor
      ScaTraEleInternalVariableManagerElchElectrodeSTIThermo(
          int numscal, const DRT::ELEMENTS::ScaTraEleParameterElch* elchpara)
          :  // call base class constructors
            ScaTraEleInternalVariableManagerElchElectrode<NSD, NEN>(numscal, elchpara),
            ScaTraEleInternalVariableManagerSTIThermo<NSD, NEN>(){};


      //! set internal variables for element evaluation
      void set_internal_variables(const CORE::LINALG::Matrix<NEN, 1>& funct,  //!< shape functions
          const CORE::LINALG::Matrix<NSD, NEN>& derxy,  //!< spatial derivatives of shape functions
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              ephinp,  //!< nodal concentration and electric potential values at time t_(n+1) or
                       //!< t_(n+alpha_F)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              ephin,  //!< nodal concentration and electric potential values at time t_(n)
          const CORE::LINALG::Matrix<NEN, 1>&
              etempnp,  //!< nodal temperature values at time t_(n+1) or t_(n+alpha_F)
          const CORE::LINALG::Matrix<NSD, NEN>&
              econvelnp,  //!< nodal convective velocity values at time t_(n+1) or t_(n+alpha_F)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>& ehist  //!< nodal history values
      )
      {
        // set thermo variables
        vmthermo::set_internal_variables_sti_thermo(funct, derxy, etempnp);

        // set scatra variables
        // this requires the temperature to be already set
        vmelchelectrode::set_internal_variables_elch_electrode(
            funct, derxy, ephinp, ephin, econvelnp, ehist);
      }

      //! set factor F/RT
      void SetFRT() override
      {
        const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
        const double gasconstant =
            DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();

        vmelch::frt_ = faraday / (gasconstant * vmthermo::Temp());
      }

      //! get GP temperature
      double Temp() { return vmthermo::Temp(); }
    };  // class ScaTraEleInternalVariableManagerElchElectrodeSTIThermo
  }     // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
