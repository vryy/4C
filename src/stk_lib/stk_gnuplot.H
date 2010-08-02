#ifdef STKADAPTIVE

#ifndef STK_GNUPLOT_H
#define STK_GNUPLOT_H

#include <string>

#include <stk_mesh/base/Field.hpp>


namespace STK
{
  class Mesh;

  void GnuplotDump( STK::Mesh & mesh, std::string basename, int step, int version, stk::mesh::FieldBase & field );

  void GnuplotDump( STK::Mesh & mesh, std::string filename, stk::mesh::FieldBase & field );
}

#endif
#endif
