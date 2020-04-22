/*----------------------------------------------------------------------*/
/*! \file

\brief manage access to and provide data globally in immersed problems

\level 3

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/
#include "immersed_field_exchange_manager.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

///----------------------------------------------------------------------*/
/// the instance
///----------------------------------------------------------------------*/
DRT::ImmersedFieldExchangeManager* DRT::ImmersedFieldExchangeManager::instance_;

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
DRT::ImmersedFieldExchangeManager* DRT::ImmersedFieldExchangeManager::Instance(bool create)
{
  if (create)
  {
    if (instance_ == NULL)
    {
      instance_ = new ImmersedFieldExchangeManager();
    }
  }
  else
  {
    if (instance_ != NULL) delete instance_;
    instance_ = NULL;
  }
  return instance_;
}

//----------------------------------------------------------------------*/
//    destruction method
//----------------------------------------------------------------------*/
void DRT::ImmersedFieldExchangeManager::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}
