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
      wear_both_map,  ///< map slave-side wear to master side
      wear_both_discr
    };

    /// Definition of contact wear algorithm
    /// (this enum represents the input file parameter WEARTYPE)
    enum WearType
    {
      wear_impl,
      wear_expl,
      wear_discr
    };

    /// Definition of wear shape functions
    /// (this enum represents the input file parameter WEAR_SHAPEFCN)
    enum WearShape
    {
      wear_shape_dual,
      wear_shape_standard
    };

    /// Definition of wear-ALE coupling algorithm
    /// (this enum represents the input file parameter WEAR_COUPALGO)
    enum WearCoupAlgo
    {
      wear_stagg,
      wear_iterstagg,
      wear_monolithic
    };

    /// Definition of wear-ALE time scale coupling algorithm
    /// (this enum represents the input file parameter WEAR_TIMESCALE)
    enum WearTimeScale
    {
      wear_time_equal,
      wear_time_different
    };

    /// Definition of configuration for wear coefficient
    /// (this enum represents the input file parameter WEARCOEFF_CONF)
    enum WearConf
    {
      wear_conf_mat,
      wear_conf_sp
    };
  }  // namespace WEAR
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_WEAR_H
