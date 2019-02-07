/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_elch_diffcond.cpp

\brief singleton class holding all static diffusion-conduction parameters required for element
evaluation

This singleton class holds all static diffusion-conduction parameters required for element
evaluation. All parameters are usually set only once at the beginning of a simulation, namely during
initialization of the global time integrator, and then never touched again throughout the
simulation. This parameter class needs to coexist with more general parameter classes holding
additional static parameters required for scalar transport element evaluation.

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251

\level 2
</pre>
*/
/*----------------------------------------------------------------------*/
#include "scatra_ele_parameter_elch_diffcond.H"
#include "scatra_ele_parameter_std.H"

#include <map>
#include <Teuchos_ParameterList.hpp>

#include "../drt_inpar/inpar_elch.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElchDiffCond*
DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::Instance(
    const std::string& disname,                      //!< name of discretization
    const ScaTraEleParameterElchDiffCond* delete_me  //!< creation/destruction indication
)
{
  // each discretization is associated with exactly one instance of this class according to a static
  // map
  static std::map<std::string, ScaTraEleParameterElchDiffCond*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if
  // not
  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleParameterElchDiffCond(disname);
  }

  // destruct instance
  else
  {
    for (std::map<std::string, ScaTraEleParameterElchDiffCond*>::iterator i = instances.begin();
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
void DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::Done()
{
  // delete singleton
  Instance("", this);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::ScaTraEleParameterElchDiffCond(
    const std::string& disname  //!< name of discretization
    )
    : cursolvar_(false),
      diffusioncoefbased_(true),
      newmanconsta_(2.),
      newmanconstb_(-2.),
      newmanconstc_(-1.)
{
  return;
}


/*---------------------------------------------------------------------*
 | set parameters                                           fang 02/15 |
 *---------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::SetParameters(
    Teuchos::ParameterList& parameters  //!< parameter list
)
{
  // access parameter sublist for diffusion-conduction formulation
  Teuchos::ParameterList& diffcondparams = parameters.sublist("DIFFCOND");

  // flag if current is used as a solution variable
  cursolvar_ = DRT::INPUT::IntegralValue<int>(diffcondparams, "CURRENT_SOLUTION_VAR");

  // mat_diffcond: flag if diffusion potential is based on diffusion coefficients or transference
  // number
  diffusioncoefbased_ =
      DRT::INPUT::IntegralValue<INPAR::ELCH::EquPot>(diffcondparams, "MAT_DIFFCOND_DIFFBASED");

  // switch for dilute and concentrated solution theory (diffusion potential in current equation):
  //    A          B
  //   |--|  |----------|
  //   z_1 + (z_2 - z_1) t_1
  // ------------------------ (RT/F kappa 1/c_k grad c_k)
  //      z_1 z_2
  //     |________|
  //         C
  newmanconsta_ = diffcondparams.get<double>("MAT_NEWMAN_CONST_A");
  newmanconstb_ = diffcondparams.get<double>("MAT_NEWMAN_CONST_B");
  newmanconstc_ = diffcondparams.get<double>("MAT_NEWMAN_CONST_C");

  return;
}
