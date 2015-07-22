/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_elch_diffcond.cpp

\brief element parameter class for electrochemistry problems governed by diffusion-conduction formulation

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_dserror.H"

#include "scatra_ele_parameter_elch_diffcond.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElchDiffCond* DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::Instance(const std::string& disname, bool create)
{
  static std::map<std::string,ScaTraEleParameterElchDiffCond* >  instances;

  if(create)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleParameterElchDiffCond(disname);
  }

  else if(instances.find(disname) != instances.end())
  {
    for( std::map<std::string,ScaTraEleParameterElchDiffCond* >::iterator i=instances.begin(); i!=instances.end(); ++i )
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
void DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::Done()
{
  // delete singleton
  Instance( "", false);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::ScaTraEleParameterElchDiffCond(const std::string& disname)
  : DRT::ELEMENTS::ScaTraEleParameterElch::ScaTraEleParameterElch(disname),
  cursolvar_(false),
  diffusioncoefbased_(true),
  newmanconsta_(2.),
  newmanconstb_(-2.),
  newmanconstc_(-1.)
{
  return;
}


/*-------------------------------------------------------------------------------------------------*
 | set parameters for electrochemistry problems with diffusion-conduction formulation   fang 02/15 |
 *-------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::SetElementGeneralParameters(Teuchos::ParameterList& params)
{
  // call base class routine
  myelch::SetElementGeneralParameters(params);

  // safety checks
  if(stabtype_ !=INPAR::SCATRA::stabtype_no_stabilization or whichtau_ != INPAR::SCATRA::tau_zero)
    dserror("No stabilization available for the diffusion-conduction formulation \n"
            "since we had no problems so far.");
  if (mat_gp_ == false or tau_gp_ == false)
    dserror("Since most of the materials of the Diffusion-conduction formulation depend on the concentration,\n"
            "an evaluation of the material and the stabilization parameter at the element center is disabled.");

  // access parameter sublist for diffusion-conduction formulation
  Teuchos::ParameterList& diffcondparams = params.sublist("DIFFCOND");

  // flag if current is used as a solution variable
  cursolvar_ = DRT::INPUT::IntegralValue<int>(diffcondparams,"CURRENT_SOLUTION_VAR");

  // mat_diffcond: flag if diffusion potential is based on diffusion coefficients or transference number
  diffusioncoefbased_ = DRT::INPUT::IntegralValue<INPAR::ELCH::EquPot>(diffcondparams,"MAT_DIFFCOND_DIFFBASED");

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
