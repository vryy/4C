/*--------------------------------------------------------------------------*/
/*!
\file reynolds_ele_parameter.cpp

\brief singleton class holding all static parameters required for Reynolds element evaluation

This singleton class holds all static parameters required for Reynolds element evaluation. All
parameters are usually set only once at the beginning of a simulation, namely during initialization of the global
time integrator, and then never touched again throughout the simulation. This parameter class needs to coexist with
the general parameter class holding all general static parameters required for Reynolds element evaluation.

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "reynolds_ele_parameter.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 | singleton access method                                  wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ReynoldsEleParameter* DRT::ELEMENTS::ReynoldsEleParameter::Instance(
    const std::string&           disname,   //!< name of discretization
    const ReynoldsEleParameter* delete_me  //!< creation/destruction indication
    )
{
  // each discretization is associated with exactly one instance of this class according to a static map
  static std::map<std::string,ReynoldsEleParameter*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if not
  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ReynoldsEleParameter(disname);
  }

  // destruct instance given to the destructor
  else
  {
    for(std::map<std::string,ReynoldsEleParameter*>::iterator i=instances.begin(); i!=instances.end(); ++i)
      if ( i->second == delete_me )
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  // return existing or newly created instance
  return instances[disname];
}

/*----------------------------------------------------------------------*
 | singleton destruction                                    wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ReynoldsEleParameter::Done()
{
  // delete singleton
  Instance("",this);

  return;
}

/*----------------------------------------------------------------------*
 | private constructor for singletons                       wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ReynoldsEleParameter::ReynoldsEleParameter(
    const std::string& disname   //!< name of discretization
    ) :
    time_(-1.0)
{
  return;
}

//----------------------------------------------------------------------*/
// set parameters which are equal for every fluid           wirtz 10/15 |
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::ReynoldsEleParameter::SetTimeParameters(
    Teuchos::ParameterList& parameters   //!< parameter list
    )
{
  // get current time and time-step length
  time_ = parameters.get<double>("total time");
}

//----------------------------------------------------------------------*/
// set parameters which are equal for every fluid           wirtz 10/15 |
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::ReynoldsEleParameter::SetGeneralParameters(
    Teuchos::ParameterList& parameters   //!< parameter list
    )
{

}
