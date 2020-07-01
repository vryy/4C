/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for level-set problems

\level 2


*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_LEVELSET_H
#define INPAR_LEVELSET_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace SCATRA
  {
    /// type of reinitialization for level set function
    enum ReInitialAction
    {
      reinitaction_none,
      reinitaction_signeddistancefunction,
      reinitaction_sussman,
      reinitaction_ellipticeq
    };

    /// type of signum function for reinitialization based on solving Sussman's equation
    enum SmoothedSignType
    {
      signtype_nonsmoothed,
      signtype_SussmanFatemi1999,
      signtype_SussmanSmerekaOsher1994,
      signtype_PengEtAl1999
    };

    /// form of linearization for reinitialization based on solving Sussman's equation
    enum LinReinit
    {
      newton,
      fixed_point
    };

    /// type of characteristic element length for reinitialization based on solving Sussman's
    /// equation
    enum CharEleLengthReinit
    {
      root_of_volume_reinit,
      streamlength_reinit
    };

    /// form of velocity for reinitialization based on solving Sussman's equation
    enum VelReinit
    {
      vel_reinit_integration_point_based,
      vel_reinit_node_based
    };

    /// form of artificial diffusion for reinitialization based on solving Sussman's equation
    enum ArtDiff
    {
      artdiff_none,
      artdiff_isotropic,
      artdiff_crosswind
    };

    /// type of reinitialization convergence check
    enum ReInitialStationaryCheck
    {
      reinit_stationarycheck_L1normintegrated,
      reinit_stationarycheck_numsteps
    };

    /// compute error compared to analytical solution
    enum CalcErrorLevelSet
    {
      calcerror_no_ls,
      calcerror_initial_field
    };

    /// problem dimension in case of quasi 2D
    enum LSDim
    {
      ls_3D,
      ls_2Dx,
      ls_2Dy,
      ls_2Dz
    };

    /// diffusivity for elliptic reinitialization
    enum DiffFunc
    {
      hyperbolic,
      hyperbolic_smoothed_positive,
      hyperbolic_clipped_05,
      hyperbolic_clipped_1
    };

  }  // namespace SCATRA

  namespace LEVELSET
  {
    /// set the levelset parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific level set conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);
  }  // namespace LEVELSET
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_LEVELSET_H
