/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_elch.cpp

\brief singleton class holding all static electrochemistry parameters required for element
evaluation

This singleton class holds all static electrochemistry parameters required for element evaluation.
All parameters are usually set only once at the beginning of a simulation, namely during
initialization of the global time integrator, and then never touched again throughout the
simulation. This parameter class needs to coexist with the general parameter class holding all
general static parameters required for scalar transport element evaluation.

\level 2

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_dserror.H"

#include "scatra_ele_parameter_elch.H"


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElch* DRT::ELEMENTS::ScaTraEleParameterElch::Instance(
    const std::string& disname,              //!< name of discretization
    const ScaTraEleParameterElch* delete_me  //!< creation/destruction indication
)
{
  // each discretization is associated with exactly one instance of this class according to a static
  // map
  static std::map<std::string, ScaTraEleParameterElch*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if
  // not
  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleParameterElch(disname);
  }

  // destruct instance
  else
  {
    for (std::map<std::string, ScaTraEleParameterElch*>::iterator i = instances.begin();
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
 | singleton destruction                                     fang 02/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterElch::Done()
{
  // delete singleton
  Instance("", this);

  return;
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 02/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElch::ScaTraEleParameterElch(
    const std::string& disname  //!< name of discretization
    )
    : boundaryfluxcoupling_(true),
      equpot_(INPAR::ELCH::equpot_undefined),
      faraday_(0.),
      gas_constant_(0.),
      epsilon_(INPAR::ELCH::epsilon_const),
      frt_(0.)
{
  return;
}


/*----------------------------------------------------------------------*
 | set parameters                                            fang 02/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterElch::SetParameters(
    Teuchos::ParameterList& parameters  //!< parameter list
)
{
  // coupling of lithium-ion flux density and electric current density at Dirichlet and Neumann
  // boundaries
  boundaryfluxcoupling_ = parameters.get<bool>("boundaryfluxcoupling");

  // type of closing equation for electric potential
  equpot_ = DRT::INPUT::get<INPAR::ELCH::EquPot>(parameters, "equpot");
  if (equpot_ == INPAR::ELCH::equpot_undefined)
    dserror("Invalid type of closing equation for electric potential!");

  // get parameters
  faraday_ = parameters.get<double>("faraday", -1.0);
  gas_constant_ = parameters.get<double>("gas_constant", -1.0);
  frt_ = parameters.get<double>("frt", -1.0);

  // safety checks
  if (frt_ <= 0.) dserror("Factor F/RT is non-positive!");
  if (faraday_ <= 0.) dserror("Faraday constant is non-positive!");
  if (gas_constant_ <= 0.) dserror("(universal) gas constant is non-positive!");

  return;
}
