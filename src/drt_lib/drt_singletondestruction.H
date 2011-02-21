#ifndef DRT_SINGLETONDESTRUCTION_H
#define DRT_SINGLETONDESTRUCTION_H

namespace DRT
{

  /// base class of singletons that need an explicit destruction
  ///
  /// There is a cleanup phase when DRT::Problem goes down. Some singleton
  /// object need to do cleanup work for themselves. Thus these classes can
  /// register with DRT::Problem via this callback.
  class SingletonDestruction
  {
  public:

    SingletonDestruction();

    virtual ~SingletonDestruction() {}

    /// called upon destruction
    virtual void Done() {}
  };
}

#endif
