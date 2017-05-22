/*----------------------------------------------------------------------*/
/*!
\file binning_strategy_utils.H

\brief utils class for use of binning strategy

\level 2

\maintainer Jonas Eichinger
*----------------------------------------------------------------------*/


#ifndef BINNING_STRATEGY_UTILS_H
#define BINNING_STRATEGY_UTILS_H

#include "Teuchos_RCP.hpp"

// forward declarations
class Epetra_Vector;
class Epetra_Map;
namespace DRT
{
  class Discretization;
  class Element;
  class Node;
}

namespace BINSTRATEGY
{
namespace UTILS
{
  /// Type of elements assigned to a bin
  enum BinContentType
  {
    Scatra,                 ///< scatra element
    Fluid,                  ///< fluid element
    BELE3,                  ///< bele3 element
    Beam,                   ///< beam element
    RigidSphere,            ///< rigid sphere element
    Solid,                  ///< solid element (all elements derived from So_base, if more distinction is needed, split this type)

    enumsize                ///< NOTE: This entry gives the size of the enum and is used somewhere else in baci therefore put it ALWAYS in last place
  };

  void ExtendDiscretizationGhosting(
      Teuchos::RCP<DRT::Discretization> discret,
      Teuchos::RCP<Epetra_Map> const&   extendedelecolmap,
      bool                              assigndegreesoffreedom,
      bool                              initelements,
      bool                              doboundaryconditions );

  /// convert element to bin content type
  BinContentType ConvertElementToBinContentType(
      DRT::Element const * const eleptr );

}
}

#endif
