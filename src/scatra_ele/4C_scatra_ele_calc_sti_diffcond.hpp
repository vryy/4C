/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluate heat transport within binary, concentrated electrolytes on element level

\level 2

*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_CALC_STI_DIFFCOND_HPP
#define FOUR_C_SCATRA_ELE_CALC_STI_DIFFCOND_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"
#include "4C_scatra_ele_sti_elch.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    // forward declarations
    template <Core::FE::CellType distype, int probdim>
    class ScaTraEleCalcElch;
    template <Core::FE::CellType distype, int probdim>
    class ScaTraEleCalcElchDiffCond;
    class ScaTraEleDiffManagerElchDiffCond;
    class ScaTraEleDiffManagerSTIThermo;
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchDiffCond;
    template <Core::FE::CellType distype>
    class ScaTraEleUtilsElchDiffCond;

    // class implementation
    template <Core::FE::CellType distype>
    class ScaTraEleCalcSTIDiffCond : public ScaTraEleCalc<distype>, public ScaTraEleSTIElch<distype>
    {
     public:
      //! singleton access method
      static ScaTraEleCalcSTIDiffCond<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);


     private:
      //! abbreviations
      using my = ScaTraEleCalc<distype>;
      using myelch = ScaTraEleCalcElch<distype, Core::FE::dim<distype>>;
      using mystielch = ScaTraEleSTIElch<distype>;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

      //! private constructor for singletons
      ScaTraEleCalcSTIDiffCond(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate action for off-diagonal system matrix block
      int EvaluateActionOD(Core::Elements::Element* ele,    //!< current element
          Teuchos::ParameterList& params,                   //!< parameter list
          Core::FE::Discretization& discretization,         //!< discretization
          const ScaTra::Action& action,                     //!< action parameter
          Core::Elements::Element::LocationArray& la,       //!< location array
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
          Core::LinAlg::SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
          Core::LinAlg::SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
          Core::LinAlg::SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
          ) override;

      //! calculate element matrix and element right-hand side vector
      void sysmat(Core::Elements::Element* ele,       ///< current element
          Core::LinAlg::SerialDenseMatrix& emat,      ///< element matrix
          Core::LinAlg::SerialDenseVector& erhs,      ///< element right-hand side vector
          Core::LinAlg::SerialDenseVector& subgrdiff  ///< subgrid diffusivity scaling vector
          ) override;

      //! element matrix and right-hand side vector contributions arising from Joule's heat
      void calc_mat_and_rhs_joule(Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix
          Core::LinAlg::SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
          ) override;

      //! element matrix and right-hand side vector contributions arising from Joule's heat
      void calc_mat_and_rhs_joule_solid(Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix
          Core::LinAlg::SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
      );

      //! element matrix and right-hand side vector contributions arising from heat of mixing
      void calc_mat_and_rhs_mixing(Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix
          Core::LinAlg::SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
          ) override;

      //! element matrix and right-hand side vector contributions arising from Soret effect
      void calc_mat_and_rhs_soret(Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix
          Core::LinAlg::SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
          ) override;

      //! fill element matrix with linearizations of discrete thermo residuals w.r.t. scatra dofs
      void sysmat_od_thermo_scatra(Core::Elements::Element* ele,  //!< current element
          Core::LinAlg::SerialDenseMatrix& emat                   //!< element matrix
      );

      //! provide element matrix with linearizations of Joule's heat term in discrete thermo
      //! residuals w.r.t. scatra dofs
      void calc_mat_joule_od(Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
          ) override;

      //! provide element matrix with linearizations of Joule's heat term in discrete thermo
      //! residuals w.r.t. scatra dofs
      void calc_mat_joule_solid_od(Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
      );

      //! provide element matrix with linearizations of heat of mixing term in discrete thermo
      //! residuals w.r.t. scatra dofs
      void calc_mat_mixing_od(Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
          ) override;

      //! provide element matrix with linearizations of Soret effect term in discrete thermo
      //! residuals w.r.t. scatra dofs
      void calc_mat_soret_od(Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
          ) override;

      //! extract quantities for element evaluation
      void extract_element_and_node_values(Core::Elements::Element* ele,  //!< current element
          Teuchos::ParameterList& params,                                 //!< parameter list
          Core::FE::Discretization& discretization,                       //!< discretization
          Core::Elements::Element::LocationArray& la                      //!< location array
          ) override;

      //! get material parameters
      void get_material_params(const Core::Elements::Element* ele,  //!< current element
          std::vector<double>& densn,                               //!< density at t_(n)
          std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,  //!< density at t_(n+alpha_M)
          double& visc,                 //!< fluid viscosity
          const int iquad = -1          //!< ID of current integration point
          ) override;

      //! evaluate Soret material
      void mat_soret(const Teuchos::RCP<const Core::Mat::Material> material,  //!< Soret material
          double& densn,   //!< density at time t_(n)
          double& densnp,  //!< density at time t_(n+1) or t_(n+alpha_F)
          double& densam   //!< density at time t_(n+alpha_M)
      );

      void mat_fourier(
          const Teuchos::RCP<const Core::Mat::Material> material,  //!< Fourier material
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

      //! diffusion manager for diffusion-conduction formulation
      Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> diffmanagerdiffcond_;

      //! utility class supporting element evaluation for diffusion-conduction formulation
      Discret::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>* utils_;
    };  // class ScaTraEleCalcSTIDiffCond
  }     // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
