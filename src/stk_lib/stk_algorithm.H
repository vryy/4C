#ifndef STK_ALGORITHM_H
#define STK_ALGORITHM_H

#ifdef STKADAPTIVE

#include <stk_mesh/base/MetaData.hpp>

namespace STK
{

  /// algorithm base class
  class Algorithm
  {
  public:

    virtual ~Algorithm() {}

    virtual void declare_fields( stk::mesh::MetaData & meta ) = 0;

    virtual void collect_unknowns( std::vector<stk::mesh::FieldBase*> & fields ) = 0;

    virtual void notify_state_changed() = 0;
  };
}

#endif
#endif
