/*----------------------------------------------------------------------*/
/*! \file

\brief preprocessor for exodusII format


\level 1

Pre_exodus contains classes to open and preprocess exodusII files into the
drt of Baci. It uses the "valid-parameters"-list defined in Baci for preparing
a up-to-date Baci header and another file specifying element and boundary
specifications. As result either a preliminary input file set is suggestioned,
or the well-known .dat file is created.

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_PRE_EXODUS_HPP
#define FOUR_C_PRE_EXODUS_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace EXODUS
{
  // forward declaration
  class Mesh;

  //! create default bc file
  int CreateDefaultBCFile(EXODUS::Mesh& mesh);

}  // namespace EXODUS

FOUR_C_NAMESPACE_CLOSE

#endif
