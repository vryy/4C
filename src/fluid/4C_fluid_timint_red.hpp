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

namespace ADAPTER
{
  class ArtNet;
}
namespace AIRWAY
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
    template <class red_D_time_int>
    class FluidCouplingWrapper;
  }  // namespace UTILS

  class TimIntRedModels : public virtual FluidImplicitTimeInt
  {
   public:
    /// Standard Constructor
    TimIntRedModels(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void Init() override;

    /*!
    \brief evaluate problem-dependent bounadry consitions: update the 3D-to-reduced_D coupling data
    for dirichlet bc in this context

    */
    void DoProblemSpecificBoundaryConditions() override;

    /*!
    \brief update the 3D-to-reduced_D coupling data in AssembleMatAndRHS

    */
    virtual void Update3DToReducedMatAndRHS();

    /*!
    \brief read restart data
    */
    void ReadRestart(int step) override;

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
    void SetupMeshtying() override;

    /*!
    \brief update configuration and output to file/screen

    */
    void Output() override;

    /*!
    \brief Insert Womersley condition

    */
    void InsertVolumetricSurfaceFlowCondVector(
        Teuchos::RCP<Epetra_Vector> vel, Teuchos::RCP<Epetra_Vector> res) override;

    /// prepare AVM3-based scale separation
    void AVM3Preparation() override;

    /// prepare time step
    void PrepareTimeStep() override;

    /*!
    \brief Additional function for RedModels in LinearRelaxationSolve

    */
    void CustomSolve(Teuchos::RCP<Epetra_Vector> relax) override;

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams) override;

    /*!
    \brief call elements to calculate system matrix/rhs and assemble

    */
    void AssembleMatAndRHS() override;

    /*!
    \brief apply Dirichlet boundary conditions to system of equations

    */
    void ApplyDirichletToSystem() override;


   protected:
    /// bio related special (in/outflow) traction velocity component adder
    Teuchos::RCP<UTILS::TotalTractionCorrector> traction_vel_comp_adder_bc_;

    /// bio related, 3D to reduced-D coupling
    Teuchos::RCP<UTILS::FluidCouplingWrapperBase> coupled3D_redDbc_art_;

    /// 1D arterial network time integration
    Teuchos::RCP<ADAPTER::ArtNet> ART_timeInt_;

    /// bio related, 3D to reduced-D coupling
    Teuchos::RCP<UTILS::FluidCouplingWrapperBase> coupled3D_redDbc_airways_;

    /// 1D arterial network time integration
    Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt> airway_imp_timeInt_;

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
