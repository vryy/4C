/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_std.cpp

\brief Setting of general scatra parameter for element evaluation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_parameter_std.H"

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterStd* DRT::ELEMENTS::ScaTraEleParameterStd::Instance( const std::string& disname, bool create )
{
  static std::map<std::string,ScaTraEleParameterStd* >  instances;

  if(create)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleParameterStd(disname);
  }

  else if(instances.find(disname) != instances.end())
  {
    for( std::map<std::string,ScaTraEleParameterStd* >::iterator i=instances.begin(); i!=instances.end(); ++i )
     {
      delete i->second;
      i->second = NULL;
     }

    instances.clear();
    return NULL;
  }

  return instances[disname];
}

//----------------------------------------------------------------------*/
//    destruction method
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterStd::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( "",false );
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterStd::ScaTraEleParameterStd(const std::string& disname)
  : DRT::ELEMENTS::ScaTraEleParameter::ScaTraEleParameter(disname)
{
}


