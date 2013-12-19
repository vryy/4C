/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_std.H

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
DRT::ELEMENTS::ScaTraEleParameterStd* DRT::ELEMENTS::ScaTraEleParameterStd::Instance( bool create )
{
  static ScaTraEleParameterStd* instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleParameterStd();
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
void DRT::ELEMENTS::ScaTraEleParameterStd::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( false );
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterStd::ScaTraEleParameterStd()
  : DRT::ELEMENTS::ScaTraEleParameter::ScaTraEleParameter()
{
}


