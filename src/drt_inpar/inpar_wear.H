/*----------------------------------------------------------------------*/
/*!
\file inpar_wear.H

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
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
      wear_none,     ///< no wear
      wear_archard   ///< Archard wear law
    };

    /// Definition of contact wear surface
    /// (this enum represents the input file parameter BOTH_SIDED_WEAR)
    enum WearSide
    {
      wear_slave,     ///< wear on slave side
      wear_both_map,  ///< map slave-side wear to master side per mortar matrices
      wear_both_discr ///< primary variable wear approach for master and slave surface
    };

    /// Definition of contact wear algorithm
    /// (this enum represents the input file parameter WEARTYPE)
    enum WearType
    {
      wear_intstate_impl,  ///< internal state variable approach for wear (impl)
      wear_intstate_expl,  ///< internal state variable approach for wear (expl)
      wear_primvar         ///< primary variable approach for wear (impl)
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
      wear_stagg,        ///< partitioned (fractional step) approach
      wear_iterstagg,    ///< partitioned approach
      wear_monolithic    ///< monolithic approach not (yet?) implemented
    };

    /// Definition of wear-ALE time scale coupling algorithm
    /// (this enum represents the input file parameter WEAR_TIMESCALE)
    enum WearTimeScale
    {
      wear_time_equal,    ///< shape evolution step after each structural step
      wear_time_different ///< shape evolution for accumulated wear after predefined structural steps
    };

    /// Definition of configuration for wear coefficient
    /// (this enum represents the input file parameter WEARCOEFF_CONF)
    enum WearConf
    {
      wear_conf_mat,   ///< wear coefficient in material conf. constant
      wear_conf_sp     ///< wear coefficient in spatial conf. constant
    };
  }  // namespace WEAR
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_WEAR_H
