#ifndef INPAR_PARAMETERLIST_UTILS_H
#define INPAR_PARAMETERLIST_UTILS_H

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

namespace DRT
{
  namespace INPUT
  {
    
    /// cast internal integer to enum type
    template <class T>
    T IntegralValue( const Teuchos::ParameterList & params, const std::string & name )
    {
      int value = Teuchos::getIntegralValue<int>( params, name );
      return static_cast<T>( value );
    }
  }
}

#endif
