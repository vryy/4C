/*----------------------------------------------------------------------*/
/*! \file

\brief preprocessor for exodusII format

\maintainer Martin Kronbichler

\level 1

Pre_exodus contains classes to open and preprocess exodusII files into the
drt of Baci. It uses the "valid-parameters"-list defined in Baci for preparing
a up-to-date Baci header and another file specifying element and boundary
specifications. As result either a preliminary input file set is suggestioned,
or the well-known .dat file is created.

*/
/*----------------------------------------------------------------------*/
#ifndef PRE_EXODUS_H_
#define PRE_EXODUS_H_

namespace EXODUS
{
  // forward declaration
  class Mesh;

  //! create default bc file
  int CreateDefaultBCFile(EXODUS::Mesh& mesh);

}  // namespace EXODUS

#endif /*PRE_EXODUS_H_*/
