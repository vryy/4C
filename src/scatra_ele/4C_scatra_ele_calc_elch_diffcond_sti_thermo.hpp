/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for thermodynamic diffusion-conduction ion-transport equations

\level 2

*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_CALC_ELCH_DIFFCOND_STI_THERMO_HPP
#define FOUR_C_SCATRA_ELE_CALC_ELCH_DIFFCOND_STI_THERMO_HPP

#include "4C_config.hpp"

#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_scatra_ele_calc_elch_diffcond.hpp"
#include "4C_scatra_ele_sti_thermo.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchDiffCondSTIThermo;

    // class implementation
    template <CORE::FE::CellType distype>
    class ScaTraEleCalcElchDiffCondSTIThermo : public ScaTraEleCalcElchDiffCond<distype>,
                                               public ScaTraEleSTIThermo<distype>
    {
     public:
      //! singleton access method
      static ScaTraEleCalcElchDiffCondSTIThermo<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);


     private:
      //! abbreviations
      using my = ScaTraEleCalc<distype>;
      using myelch = ScaTraEleCalcElch<distype>;
      using myelectrode = ScaTraEleCalcElchElectrode<distype>;
      using mydiffcond = ScaTraEleCalcElchDiffCond<distype>;
      using mythermo = ScaTraEleSTIThermo<distype>;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

      //! private constructor for singletons
      ScaTraEleCalcElchDiffCondSTIThermo(
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
          const int iquad = -1                           //!< ID of current integration point
          ) override;

      //! calculate element matrix and element right-hand side vector
      void calc_mat_and_rhs(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  //!< element right-hand side vector
          const int k,                            //!< index of current scalar
          const double fac,                       //!< domain integration factor
          const double timefacfac,  //!< domain integration factor times time integration factor
          const double rhsfac,      //!< domain integration factor times time integration factor for
                                    //!< right-hand side vector
          const double taufac,      //!< domain integration factor times stabilization parameter
          const double timetaufac,  //!< domain integration factor times stabilization parameter
                                    //!< times time integration factor
          const double rhstaufac,   //!< domain integration factor times stabilization parameter
                                    //!< times time integration factor for right-hand side vector
          CORE::LINALG::Matrix<nen_, 1>&
              tauderpot,  //!< derivatives of stabilization parameter w.r.t. electric potential
          double& rhsint  //!< body force value
          ) override;

      //! fill element matrix with linearizations of discrete scatra residuals w.r.t. thermo dofs
      void sysmat_od_scatra_thermo(DRT::Element* ele,  //!< current element
          CORE::LINALG::SerialDenseMatrix& emat        //!< element matrix
      );

      //! set internal variables for element evaluation
      void set_internal_variables_for_mat_and_rhs() override;

      //! get internal variable manager for thermodynamic diffusion-conduction formulation
      Teuchos::RCP<ScaTraEleInternalVariableManagerElchDiffCondSTIThermo<nsd_, nen_>> var_manager()
      {
        return Teuchos::rcp_static_cast<
            ScaTraEleInternalVariableManagerElchDiffCondSTIThermo<nsd_, nen_>>(
            my::scatravarmanager_);
      };

     private:
      // material type for evaluation
      CORE::Materials::MaterialType materialtype_;

    };  // class ScaTraEleCalcElchDiffCondSTIThermo


    //! implementation of ScaTraEleInternalVariableManagerElchDiffCondSTIThermo
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchDiffCondSTIThermo
        : public ScaTraEleInternalVariableManagerElchDiffCond<NSD, NEN>,
          public ScaTraEleInternalVariableManagerSTIThermo<NSD, NEN>
    {
     public:
      //! abbreviations
      using vmelch = ScaTraEleInternalVariableManagerElch<NSD, NEN>;
      using vmdiffcond = ScaTraEleInternalVariableManagerElchDiffCond<NSD, NEN>;
      using vmthermo = ScaTraEleInternalVariableManagerSTIThermo<NSD, NEN>;

      //! constructor
      ScaTraEleInternalVariableManagerElchDiffCondSTIThermo(int numscal,
          const DRT::ELEMENTS::ScaTraEleParameterElch* elchparams,
          const DRT::ELEMENTS::ScaTraEleParameterElchDiffCond* diffcondparams)
          :  // call base class constructors
            ScaTraEleInternalVariableManagerElchDiffCond<NSD, NEN>(
                numscal, elchparams, diffcondparams),
            ScaTraEleInternalVariableManagerSTIThermo<NSD, NEN>(){};


      //! set internal variables for element evaluation
      void set_internal_variables(const CORE::LINALG::Matrix<NEN, 1>& funct,  //!< shape functions
          const CORE::LINALG::Matrix<NSD, NEN>& derxy,  //!< spatial derivatives of shape functions
          const CORE::LINALG::Matrix<NEN, 1>&
              etempnp,  //!< nodal temperature values at time t_(n+1) or t_(n+alpha_F)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              ephinp,  //!< nodal concentration and electric potential values at time t_(n+1) or
                       //!< t_(n+alpha_F)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              ephin,  //!< nodal concentration and electric potential values at time t_(n)
          const CORE::LINALG::Matrix<NSD, NEN>&
              econvelnp,  //!< nodal convective velocity values at time t_(n+1) or t_(n+alpha_F)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>& ehist  //!< nodal history values
      )
      {
        // set thermo variables
        vmthermo::set_internal_variables_sti_thermo(funct, derxy, etempnp);

        // set scatra variables
        // this requires the temperature to be already set
        vmdiffcond::set_internal_variables_elch_diff_cond(
            funct, derxy, ephinp, ephin, econvelnp, ehist);
      }

      //! set factor F/RT
      void SetFRT() override
      {
        vmelch::frt_ = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday() /
                       (DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant() *
                           vmthermo::Temp());
      }
    };  // class ScaTraEleInternalVariableManagerElchDiffCondSTIThermo
  }     // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
