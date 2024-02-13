/*----------------------------------------------------------------------*/
/*! \file

\brief contains base class for a generic output filter (ensight and vtk are derived from this
 class)


\level 0
*/


#include "baci_post_writer_base.hpp"

#include "baci_post_common.hpp"

BACI_NAMESPACE_OPEN

PostWriterBase::PostWriterBase(PostField* field, const std::string& filename)
    : field_(field),
      filename_(filename),
      myrank_(field->problem()->comm()->MyPID()),
      numproc_(field->problem()->comm()->NumProc())
{
}

BACI_NAMESPACE_CLOSE
