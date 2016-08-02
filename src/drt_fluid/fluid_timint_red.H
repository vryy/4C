/*!----------------------------------------------------------------------
\file fluid_timint_red.H

\brief TimIntRedModels

\level 3
<pre>
\maintainer Benjamin Krank & Martin Kronbichler
             {krank,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15252/-235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_RED_H
#define FLD_TIMINT_RED_H


#include "fluidimplicitintegration.H"

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
  class VolumetricFlowMapExtractor;
  template <class red_D_time_int>
  class Fluid_couplingWrapper;
}

class TimIntRedModels : public virtual FluidImplicitTimeInt
{
public:
  /// Standard Constructor
  TimIntRedModels(const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

  /// Destructor
  virtual ~TimIntRedModels();

  /*!
  \brief initialization

  */
  virtual void Init();

  /*!
  \brief evaluate problem-dependent bounadry consitions: update the 3D-to-reduced_D coupling data for dirichlet bc in this context

  */
  virtual void DoProblemSpecificBoundaryConditions();

  /*!
  \brief update the 3D-to-reduced_D coupling data in AssembleMatAndRHS

  */
  virtual void Update3DToReducedMatAndRHS();

  /*!
  \brief read restart data
  */
  virtual void ReadRestart(int step);

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
  virtual void SetupMeshtying();

  /*!
  \brief update configuration and output to file/screen

  */
  virtual void Output();

  /*!
  \brief Insert Womersley condition

  */
  virtual void InsertVolumetricSurfaceFlowCondVector(
      Teuchos::RCP<Epetra_Vector> vel ,
      Teuchos::RCP<Epetra_Vector> res);

  /// prepare AVM3-based scale separation
  virtual void AVM3Preparation();

  /// prepare time step
  virtual void PrepareTimeStep();

  /*!
  \brief Additional function for RedModels in LinearRelaxationSolve

  */
  virtual void CustomSolve(Teuchos::RCP<Epetra_Vector> relax);

  /*!
  \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

  */
  virtual void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams);

  /*!
  \brief call elements to calculate system matrix/rhs and assemble

  */
  virtual void AssembleMatAndRHS();

  /*!
  \brief apply Dirichlet boundary conditions to system of equations

  */
  virtual void ApplyDirichletToSystem();


protected:

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

  /// bio related special (in/outflow) boundaries
  Teuchos::RCP<UTILS::FluidVolumetricSurfaceFlowWrapper> vol_surf_flow_bc_;

  /// maps for womersley flow profile which is applied as a Dirichlet condition
  Teuchos::RCP<Epetra_Map>                    vol_surf_flow_bc_maps_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<FLD::UTILS::VolumetricFlowMapExtractor> vol_flow_rates_bc_extractor_;

  /// flag for potential 3D Reduced_D coupling
  bool strong_redD_3d_coupling_;
private:


}; // class TimIntRedModels

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_RED_H
