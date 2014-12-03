/*!----------------------------------------------------------------------
\file immersed_field_exchange_manager.cpp

\brief manage access to and provide data in immersed problems

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "immersed_field_exchange_manager.H"

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
IMMERSED::ImmersedFieldExchangeManager*  IMMERSED::ImmersedFieldExchangeManager::Instance( bool create )
{
  static ImmersedFieldExchangeManager* instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ImmersedFieldExchangeManager();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}

//----------------------------------------------------------------------*/
//    destruction method
//----------------------------------------------------------------------*/
void IMMERSED::ImmersedFieldExchangeManager::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( false );
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
IMMERSED::ImmersedFieldExchangeManager::ImmersedFieldExchangeManager()
{
  // empty
}


