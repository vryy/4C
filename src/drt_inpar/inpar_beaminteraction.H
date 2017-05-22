/*-----------------------------------------------------------*/
/*!
\file inpar_beaminteraction.H

\brief input parameter for beaminteraction

\maintainer Jonas Eichinger, Maximilian Grill

\level 2

*/
/*-----------------------------------------------------------*/
#ifndef INPAR_BEAMINTERACTION_H
#define INPAR_BEAMINTERACTION_H

#include "inpar_parameterlist_utils.H"
#include <boost/algorithm/string/predicate.hpp>

#include "../drt_lib/drt_dserror.H"

// forward declaration
namespace DRT {
namespace INPUT {
class ConditionDefinition;
}
}

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters used in statistical mechanics
  namespace BEAMINTERACTION
  {

  /// type of the used submodel for beaminteraction
  enum SubModelType
  {
    submodel_crosslinking,         ///< evaluate the structural model
    submodel_beamcontact,          ///< evaluate the contact model
    submodel_potential,            ///< evaluate the model for potential-based interactions
    submodel_contractilecells,     ///< evaluate model for cell filament interactions
    submodel_vague                 ///< undefined model type
  };

  /// type of employed solving strategy for contact
  /// (this enum represents the input file parameter STRATEGY)
  enum Strategy
  {
    bstr_none,        ///< no beam contact
    bstr_penalty      ///< penalty method
  };

  /// type of filament
  enum FilamentType
  {
    filtype_arbitrary,        ///< no special type
    filtype_actin,           ///< actin type
    filtype_collagen         ///< collagen type
  };

  /// type of crosslinker
  enum CrosslinkerType
  {
    linkertype_arbitrary,        ///< binds to all filament
    linkertype_actin,           ///< only binds to actin filaments
    linkertype_collagen         ///< only binds to collagen filaments
  };

  //! Map type std::string to enum
  inline enum FilamentType String2FilamentType(const std::string& name)
  {
    FilamentType type = filtype_arbitrary;
    if ( boost::iequals( name, "arbitrary" ) )
      type = filtype_arbitrary;
    else if ( boost::iequals( name, "actin" ) )
      type = filtype_actin;
    else if ( boost::iequals( name, "collagen" ) )
      type = filtype_collagen;
    else
      dserror("invalid filament type std::string ");

    return type;
  };

  //! Map type std::string to enum
  inline enum CrosslinkerType String2CrosslinkerType(const std::string& name)
  {
    CrosslinkerType type = linkertype_arbitrary;
    if ( boost::iequals( name, "arbitrary" ) )
      type = linkertype_arbitrary;
    else if ( boost::iequals( name, "actin" ) )
      type = linkertype_actin;
    else if ( boost::iequals( name, "collagen" ) )
      type = linkertype_collagen;
    else
      dserror("invalid crosslinker type std::string ");

    return type;
  };


  /// set the beam interaction parameters
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  /// set beam interaction specific conditions
  void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  }

}

/*----------------------------------------------------------------------*/
#endif
