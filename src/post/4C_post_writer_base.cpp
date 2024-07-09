/*----------------------------------------------------------------------*/
/*! \file

\brief contains base class for a generic output filter (ensight and vtk are derived from this
 class)


\level 0
*/


#include "4C_post_writer_base.hpp"

#include "4C_post_common.hpp"

FOUR_C_NAMESPACE_OPEN

PostWriterBase::PostWriterBase(PostField* field, const std::string& filename)
    : field_(field),
      filename_(filename),
      myrank_(field->problem()->get_comm()->MyPID()),
      numproc_(field->problem()->get_comm()->NumProc())
{
}

FOUR_C_NAMESPACE_CLOSE
