/*!----------------------------------------------------------------------
\file fluid_timint_red.H
\brief

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_RED_H
#define FLD_TIMINT_RED_H


#include "fluidimplicitintegration.H"

#include "../linalg/linalg_utils.H"

namespace ART
{
  class ArtNetExplicitTimeInt;
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
  class Fluid_couplingWrapperBase;
  template <class red_D_time_int>
  class Fluid_couplingWrapper;
}

class TimIntRedAirways : public virtual FluidImplicitTimeInt
{
public:
  /// Standard Constructor
  TimIntRedAirways(const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

  /// Destructor
  virtual ~TimIntRedAirways();

  /*!
  \brief update the 3D-to-reduced_D coupling data for dirichlet bc

  */
  void Update3DToReducedDirichlet();

  /*!
  \brief update the 3D-to-reduced_D coupling data in AssembleMatAndRHS

  */
  void Update3DToReducedMatAndRHS();

  /*!
  \brief read restart data
  */
  void ReadRestart(int step);

  /*!
  \brief read restart (some more RedAirways-specific data)

  */
  void ReadRestartReducedD(int step);

  /*!
  \brief update configuration and output to file/screen

  */
  void OutputReducedD();

  /*!
  \brief Update the solution after convergence of the nonlinear
         iteration. Current solution becomes old solution of next
         timestep.
  */
  void TimeUpdate();

  /*!
  \brief Setup meshtying

  */
  void SetupMeshtying();

  /*!
  \brief update iterative increment

  */
  void UpdateIterIncrementally(
    Teuchos::RCP<const Epetra_Vector> vel  //!< input residual velocities
    );

  /*!
  \brief update configuration and output to file/screen

  */
  void Output();

  /*!
  \brief Insert Womersley condition

  */
  void InsertVolumetricSurfaceFlowCondVector(
      Teuchos::RCP<Epetra_Vector> vel ,
      Teuchos::RCP<Epetra_Vector> res);

  /// prepare AVM3-based scale separation
  void AVM3Preparation();

  /*!
  \brief Additional function for RedAirways in LinearRelaxationSolve

  */
  void RedALinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax);

  /*!
  \brief Set custom parameters in the respective time integration class (Loma, RedAirways...)

  */
  void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams);

  /*!
  \brief call elements to calculate system matrix/rhs and assemble

  */
  void AssembleMatAndRHS();

  /*!
  \brief apply Dirichlet boundary conditions to system of equations

  */
  void ApplyDirichletToSystem();

  /*!
  \brief Evaluate vel for Womersley conditions

  */
  void EvaluateWomersleyVel();


protected:

  /// bio related special (outflow) boundaries
  Teuchos::RCP<UTILS::FluidImpedanceWrapper> impedancebc_;

  /// bio related special (outflow) optimization of boundaries
  Teuchos::RCP<UTILS::FluidWkOptimizationWrapper> Wk_optimization_;

  /// bio related special (in/outflow) boundaries
  Teuchos::RCP<UTILS::FluidVolumetricSurfaceFlowWrapper> vol_surf_flow_bc_;

  /// bio related special (in/outflow) traction velocity component adder
  Teuchos::RCP<UTILS::TotalTractionCorrector> traction_vel_comp_adder_bc_;

  /// bio related, 3D to reduced-D coupling
  Teuchos::RCP<UTILS::Fluid_couplingWrapperBase> coupled3D_redDbc_art_;

  /// 1D arterial network time integration
  Teuchos::RCP<ART::ArtNetExplicitTimeInt>   ART_exp_timeInt_;

  /// bio related, 3D to reduced-D coupling
  Teuchos::RCP<UTILS::Fluid_couplingWrapperBase> coupled3D_redDbc_airways_;

  /// 1D arterial network time integration
  Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt>   airway_imp_timeInt_;

  /// maps for womersley flow profile which is applied as a Dirichlet condition
  Teuchos::RCP<Epetra_Map>                    vol_surf_flow_bcmaps_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<FLD::UTILS::MapExtractor> vol_flow_rates_bc_extractor_;

  /// flag for potential 3D Reduced_D coupling
  bool strong_redD_3d_coupling_;
private:


}; // class TimIntRedAirways

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_RED_H
