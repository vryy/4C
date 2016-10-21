/*----------------------------------------------------------------------*/
/*!
\file inpar_binstrategy.H

\brief Enums for problems using binstrategy

\level 3

\maintainer Georg Hammerl
*-----------------------------------------------------------------------*/
#ifndef INPAR_BINSTRATEGY_H_
#define INPAR_BINSTRATEGY_H_

#include "inpar_parameterlist_utils.H"


namespace INPAR
{
  /// Enums used in binstrategy
  namespace BINSTRATEGY
  {

    /// Type of elements assigned to a bin
    enum BinContent
    {
      Surface,                 ///< surface element
      Volume,                  ///< volume element
      Beam,                    ///< beam element
      enumsize                 ///< NOTE: This entry gives the size of the enum
                                 // and has therefore ALWAYS to be its last entry
    };

  } // namespace BINSTRATEGY
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif /* INPAR_BINSTRATEGY_H_ */
