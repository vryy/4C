/*----------------------------------------------------------------------*/
/*!
\file inpar_wear.H

\level 1

<pre>
\maintainer Alexander Popp
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_WEAR_H
#define INPAR_WEAR_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{
  /// Input parameters for structural wear
  namespace WEAR
  {
    /// Type of contact wear law
    /// (this enum represents the input file parameter WEAR)
    enum WearLaw
    {
      wear_none,    ///< no wear
      wear_archard  ///< Archard wear law
    };

    /// Definition of contact wear surface
    /// (this enum represents the input file parameter WEAR_SIDE)
    enum WearSide
    {
      wear_slave,  ///< wear on slave side
      wear_both    ///< slave and master wear
    };

    /// Definition of contact wear algorithm
    /// (this enum represents the input file parameter WEARTYPE)
    enum WearType
    {
      wear_intstate,  ///< internal state variable approach for wear
      wear_primvar    ///< primary variable approach for wear
    };

    /// Definition of wear time integration
    /// (this enum represents the input file parameter WEARTIMINT)
    enum WearTimInt
    {
      wear_expl,  ///< implicit time integration
      wear_impl   ///< explicit time integration
    };

    /// Definition of wear shape functions (necessary for prim. var. approach)
    /// (this enum represents the input file parameter WEAR_SHAPEFCN)
    enum WearShape
    {
      wear_shape_dual,     ///< dual shape functions allowing for condensation
      wear_shape_standard  ///< std. shape functions
    };

    /// Definition of wear-ALE coupling algorithm
    /// (this enum represents the input file parameter WEAR_COUPALGO)
    enum WearCoupAlgo
    {
      wear_stagg,      ///< partitioned (fractional step) approach
      wear_iterstagg,  ///< partitioned approach
      wear_monolithic  ///< monolithic approach not (yet?) implemented
    };

    /// Definition of wear-ALE time scale coupling algorithm
    /// (this enum represents the input file parameter WEAR_TIMESCALE)
    enum WearTimeScale
    {
      wear_time_equal,     ///< shape evolution step after each structural step
      wear_time_different  ///< shape evolution for accumulated wear after predefined structural
                           ///< steps
    };

    /// Definition of configuration for wear coefficient
    /// (this enum represents the input file parameter WEARCOEFF_CONF)
    enum WearCoeffConf
    {
      wear_coeff_mat,  ///< wear coefficient in material conf. constant
      wear_coeff_sp    ///< wear coefficient in spatial conf. constant
    };

    /// Definition of configuration for shape evolution step
    /// (this enum represents the input file parameter WEAR_SHAPE_EVO)
    enum WearShapeEvo
    {
      wear_se_mat,  ///< shape evolution in material conf.
      wear_se_sp    ///< shape evolution in spat. conf.
    };

    /// set the wear parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);
  }  // namespace WEAR
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_WEAR_H
