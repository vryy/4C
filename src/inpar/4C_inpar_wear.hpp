/*----------------------------------------------------------------------*/
/*! \file

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_WEAR_HPP
#define FOUR_C_INPAR_WEAR_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace Wear
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

    /// Definition of wear-ALE time scale coupling algorithm
    /// (this enum represents the input file parameter WEAR_TIMESCALE)
    enum WearTimeScale
    {
      wear_time_equal,     ///< shape evolution step after each structural step
      wear_time_different  ///< shape evolution for accumulated wear after predefined structural
                           ///< steps
    };

    /// set the wear parameters
    void set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list);
  }  // namespace Wear
}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
