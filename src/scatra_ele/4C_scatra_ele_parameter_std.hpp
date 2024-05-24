/*----------------------------------------------------------------------*/
/*! \file

\brief singleton class holding all static parameters required for the evaluation of a standard
scalar transport element

This singleton class holds all static parameters required for the evaluation of a standard scalar
transport element, e.g., stabilization parameters and finite difference check parameters. All
parameters are usually set only once at the beginning of a simulation, namely during initialization
of the global time integrator, and then never touched again throughout the simulation. Enhanced
scalar transport problems, such as electrochemistry and levelset problems, instantiate additional,
problem specific singleton classes holding additional static parameters required for element
evaluation. These additional singleton classes are not meant to be derived from, but rather to
coexist with this general class.

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_PARAMETER_STD_HPP
#define FOUR_C_SCATRA_ELE_PARAMETER_STD_HPP

#include "4C_config.hpp"

#include "4C_inpar_cardiac_monodomain.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_scatra_ele_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleParameterTimInt;

    // class implementation
    class ScaTraEleParameterStd : public ScaTraEleParameterBase
    {
     public:
      //! singleton access method
      static ScaTraEleParameterStd* Instance(const std::string& disname);

      //! set parameters
      void SetParameters(Teuchos::ParameterList& parameters) override;

      //! set the nodeset parameters
      void set_nodeset_parameters(Teuchos::ParameterList& parameters);

      //! @name return general parameters
      //! @{
      bool IsAle() const { return is_ale_; };
      bool IsConservative() const { return is_conservative_; };
      bool SphericalCoords() const { return sphericalcoords_; };
      INPAR::SCATRA::FluxType CalcFluxDomain() const { return calcflux_domain_; };
      Teuchos::RCP<std::vector<int>> WriteFluxIds() const { return writefluxids_; };
      INPAR::SCATRA::FdCheck fd_check() const { return fdcheck_; };
      double FDCheckEps() const { return fdcheckeps_; };
      double FDCheckTol() const { return fdchecktol_; };
      int ProbNum() const { return probnum_; };
      bool SemiImplicit() const { return semiimplicit_; };
      double int_layer_growth_conv_tol() const { return intlayergrowth_convtol_; };
      unsigned int_layer_growth_ite_max() const { return intlayergrowth_itemax_; };
      bool partitioned_multi_scale() const { return partitioned_multiscale_; };
      bool IsEMD() const { return is_emd_; };
      int EMDSource() const { return emd_source_; };
      //! return true if external force is applied
      [[nodiscard]] bool HasExternalForce() const { return has_external_force_; };
      //! number of dofset associated with displacement dofs
      int NdsDisp() const;
      //! number of dofset associated with interface growth dofs
      int NdsGrowth() const;
      //! number of dofset to write micro scale values on
      int NdsMicro() const;
      //! number of dofset associated with pressure dofs
      int NdsPres() const;
      //! number of dofset associated with scalar transport dofs
      int NdsScaTra() const;
      //! number of dofset associated with temperature dofs
      int NdsThermo() const;
      //! number of dofset associated with two-tensor quantity dofs, e.g. stresses, strains
      int nds_two_tensor_quantity() const;
      //! number of dofset associated with velocity related dofs
      int NdsVel() const;
      //! number of dofset associated with wall shear stress dofs
      int NdsWss() const;

      //! @}

      //! @name return stabilization parameters
      //! @{
      INPAR::SCATRA::StabType StabType() const { return stabtype_; };
      INPAR::SCATRA::TauType TauDef() const { return whichtau_; };
      INPAR::SCATRA::CharEleLength CharEleLength() const { return charelelength_; };
      double USFEMGLSFac() const { return diffreastafac_; };
      bool RBSubGrVel() const { return sgvel_; };
      bool ASSGD() const { return assgd_; };
      INPAR::SCATRA::AssgdType ASSGDType() const { return whichassgd_; };
      bool TauGP() const { return tau_gp_; };
      bool MatGP() const { return mat_gp_; };
      double TauValue() const { return tau_value_; };
      //! @}

     private:
      //! private constructor for singletons
      ScaTraEleParameterStd(const std::string& disname  //!< name of discretization
      );

      //! @name general parameters
      //! @{
      //! flag for ALE
      bool is_ale_;

      //! flag for conservative form
      bool is_conservative_;

      //! flag for use of spherical coordinates
      bool sphericalcoords_;

      //! flag for writing the flux vector fields
      INPAR::SCATRA::FluxType calcflux_domain_;

      //! ids of scalars for which flux vectors are written (starting with 1)
      Teuchos::RCP<std::vector<int>> writefluxids_;

      //! flag for finite difference check
      INPAR::SCATRA::FdCheck fdcheck_;

      //! perturbation magnitude for finite difference check
      double fdcheckeps_;

      //! relative tolerance for finite difference check
      double fdchecktol_;

      //! number of dofset associated with displacement dofs
      int nds_disp_;

      //! number of dofset associated with interface growth dofs
      int nds_growth_;

      //! number of dofset to write micro scale values on
      int nds_micro_;

      //! number of dofset associated with pressure dofs
      int nds_pres_;

      //! number of dofset associated with scalar transport dofs
      int nds_scatra_;

      //! number of dofset associated with temperature dofs
      int nds_thermo_;

      //! number of dofset associated with two-tensor quantity dofs, e.g. stresses, strains
      int nds_two_tensor_quantitiy_;

      //! number of dofset associated with velocity related dofs
      int nds_vel_;

      //! number of dofset associated with wall shear stress dofs
      int nds_wss_;

      //! problem number
      int probnum_;

      //! evaluation type
      bool semiimplicit_;

      //! local Newton-Raphson convergence tolerance for scatra-scatra interface coupling involving
      //! interface layer growth
      double intlayergrowth_convtol_;

      //! maximum number of local Newton-Raphson iterations for scatra-scatra interface coupling
      //! involving interface layer growth
      unsigned intlayergrowth_itemax_;

      // flag for truly partitioned multi-scale simulation
      bool partitioned_multiscale_;

      /// flag for electromagnetic diffusion simulation
      bool is_emd_;

      /// electromagnetic diffusion source function
      int emd_source_;

      /// flag for external force
      bool has_external_force_;

      //! @}

      //! @name stabilization parameters
      //! @{
      //! type of stabilization
      INPAR::SCATRA::StabType stabtype_;

      //! definition of stabilization parameter
      INPAR::SCATRA::TauType whichtau_;

      //! definition of characteristic element length
      INPAR::SCATRA::CharEleLength charelelength_;

      //! parameter to switch between SUPG, GLS and USFEM
      double diffreastafac_;

      //! flag to include residual-based subgrid-scale velocity
      bool sgvel_;

      //! flag to active artificial diffusion
      bool assgd_;

      //! definition of artificial diffusion
      INPAR::SCATRA::AssgdType whichassgd_;

      //! flag for evaluation of tau at Gauss point
      bool tau_gp_;

      //! flag for evaluation of material at Gauss point
      bool mat_gp_;

      //! numerical value for stabilization parameter
      double tau_value_;
      //! @}

      //! parameter class for time integration
      DRT::ELEMENTS::ScaTraEleParameterTimInt* scatraparatimint_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
