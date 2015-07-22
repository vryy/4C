/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_elch.cpp

\brief element parameter class for electrochemistry problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_dserror.H"

#include "scatra_ele_parameter_elch.H"


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElch* DRT::ELEMENTS::ScaTraEleParameterElch::Instance(const std::string& disname,bool create)
{
  static std::map<std::string,ScaTraEleParameterElch* >  instances;

  if(create)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleParameterElch(disname);
  }

  else if(instances.find(disname) != instances.end())
  {
    for( std::map<std::string,ScaTraEleParameterElch* >::iterator i=instances.begin(); i!=instances.end(); ++i )
     {
      delete i->second;
      i->second = NULL;
     }

    instances.clear();
    return NULL;
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 02/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterElch::Done()
{
  // delete singleton
  Instance("",false);

  return;
}

/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 02/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElch::ScaTraEleParameterElch(const std::string& disname)
  : DRT::ELEMENTS::ScaTraEleParameter::ScaTraEleParameter(disname),
  equpot_(INPAR::ELCH::equpot_undefined),
  faraday_(INPAR::ELCH::faraday_const),
  epsilon_(INPAR::ELCH::epsilon_const),
  frt_(0.)
{
  return;
}

/*----------------------------------------------------------------------*
 | set parameters for electrochemistry problems              fang 02/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterElch::SetElementGeneralParameters(Teuchos::ParameterList& params)
{
  // call base class routine
  my::SetElementGeneralParameters(params);

  // type of closing equation for electric potential
  equpot_ = DRT::INPUT::get<INPAR::ELCH::EquPot>(params, "equpot");
  if(equpot_ == INPAR::ELCH::equpot_undefined)
    dserror("Invalid type of closing equation for electric potential!");

  // get parameter F/RT
  frt_ = params.get<double>("frt");

  // safety check - only stabilization of SUPG-type available
  if ((stabtype_ !=INPAR::SCATRA::stabtype_no_stabilization) and (stabtype_ !=INPAR::SCATRA::stabtype_SUPG))
    dserror("Only SUPG-type stabilization available for ELCH.");

  return;
}
