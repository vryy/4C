#ifndef DRT_ELEMENTIMPL_H
#define DRT_ELEMENTIMPL_H

namespace DRT
{
namespace ELEMENTS
{

  /// base class of element implementation classes
  ///
  /// Element implementations are singleton classes that are never
  /// destructed. Nevertheless we might need to clean them up. This is what
  /// the Done() method is good for.
  class ElementImpl
  {
  public:

    ElementImpl();

    /// Empty destructor
    virtual ~ElementImpl() {}

    /// called upon destruction
    virtual void Done() {}
  };

}
}

#endif
