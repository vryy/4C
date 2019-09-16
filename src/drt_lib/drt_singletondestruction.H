/*---------------------------------------------------------------------*/
/*! \file

\brief Base class for singletons that need explicit destruction

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef DRT_SINGLETONDESTRUCTION_H
#define DRT_SINGLETONDESTRUCTION_H

namespace DRT
{
  /// base class of singletons that need an explicit destruction
  ///
  /// There is a cleanup phase when DRT::Problem goes down. Some singleton
  /// object need to do cleanup work for themselves. Thus these classes can
  /// register with DRT::Problem via this callback.
  ///
  /// This class is only meant for lazy initialization singletons. This is
  /// singletons that have a proper Instance() method that creates the
  /// instance at its first call. Do not use this class for static objects,
  /// since the creation order of static objects is arbitrary and the
  /// SingletonDestruction constructor needs to access the DRT::Problem
  /// instance.
  class SingletonDestruction
  {
   public:
    SingletonDestruction();

    virtual ~SingletonDestruction() {}

    /// called upon destruction (pure virtual to enforce specific implementation!)
    virtual void Done() = 0;
  };
}  // namespace DRT

#endif
