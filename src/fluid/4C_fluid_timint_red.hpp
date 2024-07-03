/*-----------------------------------------------------------*/
/*! \file

\brief Basic time integration driver for reduced models


\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_FLUID_TIMINT_RED_HPP
#define FOUR_C_FLUID_TIMINT_RED_HPP


#include "4C_config.hpp"

#include "4C_fluid_implicit_integration.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class ArtNet;
}
namespace Airway
{
  class RedAirwayImplicitTimeInt;
}

namespace FLD
{
  namespace UTILS
  {
    class MapExtractor;
    class FluidImpedanceWrapper;
    class FluidWkOptimizationWrapper;
    class FluidVolumetricSurfaceFlowWrapper;
    class TotalTractionCorrector;
    class FluidCouplingWrapperBase;
    class VolumetricFlowMapExtractor;
    template <class RedDTimeInt>
    class FluidCouplingWrapper;
  }  // namespace UTILS

  class TimIntRedModels : public virtual FluidImplicitTimeInt
  {
   public:
    /// Standard Constructor
    TimIntRedModels(const Teuchos::RCP<Core::FE::Discretization>& actdis,
        const Teuchos::RCP<Core::LinAlg::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void init() override;

    /*!
    \brief evaluate problem-dependent bounadry consitions: update the 3D-to-reduced_D coupling data
    for dirichlet bc in this context

    */
    void do_problem_specific_boundary_conditions() override;

    /*!
    \brief update the 3D-to-reduced_D coupling data in assemble_mat_and_rhs

    */
    virtual void update3_d_to_reduced_mat_and_rhs();

    /*!
    \brief read restart data
    */
    void read_restart(int step) override;

    /*!
    \brief read restart (some more RedModels-specific data)

    */
    virtual void ReadRestartReducedD(int step);

    /*!
    \brief update configuration and output to file/screen

    */
    virtual void OutputReducedD();

    /*!
    \brief Setup meshtying

    */
    void setup_meshtying() override;

    /*!
    \brief update configuration and output to file/screen

    */
    void output() override;

    /*!
    \brief Insert Womersley condition

    */
    void insert_volumetric_surface_flow_cond_vector(
        Teuchos::RCP<Epetra_Vector> vel, Teuchos::RCP<Epetra_Vector> res) override;

    /// prepare AVM3-based scale separation
    void av_m3_preparation() override;

    /// prepare time step
    void prepare_time_step() override;

    /*!
    \brief Additional function for RedModels in linear_relaxation_solve

    */
    void CustomSolve(Teuchos::RCP<Epetra_Vector> relax) override;

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    void set_custom_ele_params_assemble_mat_and_rhs(Teuchos::ParameterList& eleparams) override;

    /*!
    \brief call elements to calculate system matrix/rhs and assemble

    */
    void assemble_mat_and_rhs() override;

    /*!
    \brief apply Dirichlet boundary conditions to system of equations

    */
    void apply_dirichlet_to_system() override;


   protected:
    /// bio related special (in/outflow) traction velocity component adder
    Teuchos::RCP<UTILS::TotalTractionCorrector> traction_vel_comp_adder_bc_;

    /// bio related, 3D to reduced-D coupling
    Teuchos::RCP<UTILS::FluidCouplingWrapperBase> coupled3D_redDbc_art_;

    /// 1D arterial network time integration
    Teuchos::RCP<Adapter::ArtNet> ART_timeInt_;

    /// bio related, 3D to reduced-D coupling
    Teuchos::RCP<UTILS::FluidCouplingWrapperBase> coupled3D_redDbc_airways_;

    /// 1D arterial network time integration
    Teuchos::RCP<Airway::RedAirwayImplicitTimeInt> airway_imp_timeInt_;

    /// bio related special (in/outflow) boundaries
    Teuchos::RCP<UTILS::FluidVolumetricSurfaceFlowWrapper> vol_surf_flow_bc_;

    /// maps for womersley flow profile which is applied as a Dirichlet condition
    Teuchos::RCP<Epetra_Map> vol_surf_flow_bc_maps_;

    /// maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<FLD::UTILS::VolumetricFlowMapExtractor> vol_flow_rates_bc_extractor_;

    /// flag for potential 3D Reduced_D coupling
    bool strong_redD_3d_coupling_;

   private:
  };  // class TimIntRedModels

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
