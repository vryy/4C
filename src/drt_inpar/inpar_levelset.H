/*----------------------------------------------------------------------*/
/*!
\file inpar_levelset.H

\brief Input parameters for level-set problems

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_LEVELSET_H
#define INPAR_LEVELSET_H

#include "inpar_parameterlist_utils.H"

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
      reinitaction_sussman
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

    /// type of characteristic element length for reinitialization based on solving Sussman's equation
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

    //TODO: Ursula
    /// type of signum function for reinitialization based on solving Sussman's equation
    enum PenaltyMethod
    {
      penalty_method_none,
      penalty_method_intersection_points,
      penalty_method_akkerman
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

  }  // namespace SCATRA

  namespace PARTICLE
  {
    /// problem dimsion in case of quasi 2D
    enum ParticleDim
    {
      particle_3D,
      particle_2Dx,
      particle_2Dy,
      particle_2Dz
    };

  }  // namespace PARTICLE

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_LEVELSET_H

