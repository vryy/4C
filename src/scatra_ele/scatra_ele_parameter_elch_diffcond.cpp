/*----------------------------------------------------------------------*/
/*! \file

\brief singleton class holding all static diffusion-conduction parameters required for element
evaluation

This singleton class holds all static diffusion-conduction parameters required for element
evaluation. All parameters are usually set only once at the beginning of a simulation, namely during
initialization of the global time integrator, and then never touched again throughout the
simulation. This parameter class needs to coexist with more general parameter classes holding
additional static parameters required for scalar transport element evaluation.


\level 2
*/
/*----------------------------------------------------------------------*/
#include "scatra_ele_parameter_elch_diffcond.H"
#include "scatra_ele_parameter_std.H"

#include <map>
#include <Teuchos_ParameterList.hpp>

#include "inpar_elch.H"
#include "inpar_parameterlist_utils.H"

#include "lib_dserror.H"
#include "utils_singleton_owner.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElchDiffCond*
DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::Instance(
    const std::string& disname  //!< name of discretization
)
{
  static auto singleton_map = ::UTILS::MakeSingletonMap<std::string>(
      [](const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleParameterElchDiffCond>(
            new ScaTraEleParameterElchDiffCond(disname));
      });

  return singleton_map[disname].Instance(::UTILS::SingletonAction::create, disname);
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
