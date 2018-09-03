/*
 \file post_drt_writer_base.cpp

 \brief contains base class for a generic output filter (ensight and vtk are derived from this
 class)

 <pre>
 Maintainer: Martin Kronbichler
 kronbichler@lnm.mw.tum.de
 http://www.lnm.mw.tum.de/
 089 - 289-15235
 </pre>
 */


#include "post_writer_base.H"
#include "post_drt_common.H"


PostWriterBase::PostWriterBase(PostField* field, const std::string& filename)
    : field_(field),
      filename_(filename),
      myrank_(field->problem()->comm()->MyPID()),
      numproc_(field->problem()->comm()->NumProc())
{
}
