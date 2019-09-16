/*--------------------------------------------------------------------------*/
/*! \file

\brief singleton class holding all static parameters required for Lubrication element evaluation

This singleton class holds all static parameters required for Lubrication element evaluation. All
parameters are usually set only once at the beginning of a simulation, namely during initialization
of the global time integrator, and then never touched again throughout the simulation. This
parameter class needs to coexist with the general parameter class holding all general static
parameters required for Lubrication element evaluation.

\level 3

\maintainer Mostafa Faraji

*/
/*--------------------------------------------------------------------------*/

#include "../drt_lubrication_ele/lubrication_ele_parameter.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 | singleton access method                                  wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::LubricationEleParameter* DRT::ELEMENTS::LubricationEleParameter::Instance(
    const std::string& disname,               //!< name of discretization
    const LubricationEleParameter* delete_me  //!< creation/destruction indication
)
{
  // each discretization is associated with exactly one instance of this class according to a static
  // map
  static std::map<std::string, LubricationEleParameter*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if
  // not
  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new LubricationEleParameter(disname);
  }

  // destruct instance given to the destructor
  else
  {
    for (std::map<std::string, LubricationEleParameter*>::iterator i = instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
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
void DRT::ELEMENTS::LubricationEleParameter::Done()
{
  // delete singleton
  Instance("", this);

  return;
}

/*----------------------------------------------------------------------*
 | private constructor for singletons                       wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::LubricationEleParameter::LubricationEleParameter(
    const std::string& disname  //!< name of discretization
    )
    : time_(-1.0), modified_reynolds_(true), addsqz_(true), roughness_deviation_(0.0)
{
  return;
}

//----------------------------------------------------------------------*/
// set parameters which are equal for every lubrication     wirtz 10/15 |
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::LubricationEleParameter::SetTimeParameters(
    Teuchos::ParameterList& parameters  //!< parameter list
)
{
  // get current time and time-step length
  time_ = parameters.get<double>("total time");
}

//----------------------------------------------------------------------*/
// set parameters which are equal for every lubrication     wirtz 10/15 |
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::LubricationEleParameter::SetGeneralParameters(
    Teuchos::ParameterList& parameters  //!< parameter list
)
{
  modified_reynolds_ = parameters.get<bool>("ismodifiedrey");
  addsqz_ = parameters.get<bool>("addsqz");
  roughness_deviation_ = parameters.get<double>("roughnessdeviation");
}
