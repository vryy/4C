/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_turbulence.cpp

\brief singleton class holding all static turbulence parameters required for scalar transport
element evaluation

This singleton class holds all static turbulence parameters required for scalar transport element
evaluation. All parameters are usually set only once at the beginning of a simulation, namely during
initialization of the global time integrator, and then never touched again throughout the
simulation. This parameter class needs to coexist with the general parameter class holding all
general static parameters required for scalar transport element evaluation.

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "scatra_ele_parameter_turbulence.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 08/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterTurbulence* DRT::ELEMENTS::ScaTraEleParameterTurbulence::Instance(
    const std::string& disname,                    //!< name of discretization
    const ScaTraEleParameterTurbulence* delete_me  //!< creation/destruction indication
)
{
  // each discretization is associated with exactly one instance of this class according to a static
  // map
  static std::map<std::string, ScaTraEleParameterTurbulence*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if
  // not
  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleParameterTurbulence(disname);
  }

  // destruct instance
  else
  {
    for (std::map<std::string, ScaTraEleParameterTurbulence*>::iterator i = instances.begin();
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
 | singleton destruction                                     fang 08/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterTurbulence::Done()
{
  // delete singleton
  Instance("", this);

  return;
}

/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 08/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterTurbulence::ScaTraEleParameterTurbulence(
    const std::string& disname  //!< name of discretization
    )
    : turbmodel_(INPAR::FLUID::no_model),
      scalarforcing_(INPAR::FLUID::scalarforcing_no),
      fssgd_(false),
      whichfssgd_(INPAR::SCATRA::fssugrdiff_no),
      Cs_(0.0),
      tpn_(0.0),
      Cs_av_(false),
      Csgs_sgvel_(0.0),
      alpha_(0.0),
      calc_N_(false),
      N_vel_(0.0),
      refvel_(INPAR::FLUID::strainrate),
      reflength_(INPAR::FLUID::cube_edge),
      c_nu_(0.0),
      nwl_(false),
      nwl_scatra_(false),
      beta_(false),
      BD_gp_(false),
      Csgs_sgphi_(0.0),
      c_diff_(0.0),
      mfs_conservative_(false),
      meanCai_(0.0),
      adapt_Csgs_phi_(false),
      turbinflow_(false),
      timintparams_(DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(disname))
{
  return;
}


/*----------------------------------------------------------------------*
 | set parameters                                       rasthofer 11/11 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterTurbulence::SetParameters(
    Teuchos::ParameterList& parameters  //!< parameter list
)
{
  // get list with model-specific parameters
  Teuchos::ParameterList& turbulencelist = parameters.sublist("TURBULENCE MODEL");
  Teuchos::ParameterList& sgvisclist = parameters.sublist("SUBGRID VISCOSITY");
  Teuchos::ParameterList& mfslist = parameters.sublist("MULTIFRACTAL SUBGRID SCALES");

  // set flag for turbulence model
  if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "no_model")
    turbmodel_ = INPAR::FLUID::no_model;
  else if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "Smagorinsky")
    turbmodel_ = INPAR::FLUID::smagorinsky;
  else if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "Dynamic_Smagorinsky")
    turbmodel_ = INPAR::FLUID::dynamic_smagorinsky;
  else if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "Multifractal_Subgrid_Scales")
    turbmodel_ = INPAR::FLUID::multifractal_subgrid_scales;
  else if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "Dynamic_Vreman")
    turbmodel_ = INPAR::FLUID::dynamic_vreman;
  else
    dserror("Unknown turbulence model for scatra!");

  // define forcing for scalar field
  if (turbulencelist.get<std::string>("SCALAR_FORCING", "no") == "isotropic")
    scalarforcing_ = INPAR::FLUID::scalarforcing_isotropic;
  else if (turbulencelist.get<std::string>("SCALAR_FORCING", "no") == "mean_scalar_gradient")
    scalarforcing_ = INPAR::FLUID::scalarforcing_mean_scalar_gradient;
  else
    scalarforcing_ = INPAR::FLUID::scalarforcing_no;

  // set flag for fine-scale subgrid diffusivity and perform some checks
  whichfssgd_ = DRT::INPUT::get<INPAR::SCATRA::FSSUGRDIFF>(parameters, "fs subgrid diffusivity");
  if (whichfssgd_ == INPAR::SCATRA::fssugrdiff_artificial)
  {
    fssgd_ = true;

    // check for solver type
    if (timintparams_->IsIncremental())
      dserror(
          "Artificial fine-scale subgrid-diffusivity approach only in combination with "
          "non-incremental solver so far!");
  }
  else if (whichfssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_all or
           whichfssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_small)
  {
    fssgd_ = true;

    // check for solver type
    if (not timintparams_->IsIncremental())
      dserror(
          "Fine-scale subgrid-diffusivity approach using all/small-scale Smagorinsky model only in "
          "combination with incremental solver so far!");
  }

  // in some cases we may want to switch off the turbulence model in the scalar field
  if (not DRT::INPUT::IntegralValue<int>(turbulencelist, "TURBMODEL_LS"))
  {
    fssgd_ = false;
    whichfssgd_ = INPAR::SCATRA::fssugrdiff_no;
    turbmodel_ = INPAR::FLUID::no_model;
  }

  if (turbmodel_ != INPAR::FLUID::no_model or (timintparams_->IsIncremental() and fssgd_))
  {
    // get Smagorinsky constant and turbulent Prandtl number
    Cs_ = sgvisclist.get<double>("C_SMAGORINSKY");
    tpn_ = sgvisclist.get<double>("C_TURBPRANDTL");
    if (tpn_ <= 1.0E-16) dserror("Turbulent Prandtl number should be larger than zero!");

    Cs_av_ = DRT::INPUT::IntegralValue<int>(sgvisclist, "C_SMAGORINSKY_AVERAGED");

    if (turbmodel_ == INPAR::FLUID::dynamic_vreman)
      Cs_ = turbulencelist.get<double>("Dt_vreman", 1.0);

    // get model parameters
    if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    {
      // necessary parameters for subgrid-scale velocity estimation
      Csgs_sgvel_ = mfslist.get<double>("CSGS");
      if (mfslist.get<std::string>("SCALE_SEPARATION") == "algebraic_multigrid_operator")
        alpha_ = 3.0;
      else
        dserror("Scale-Separtion method not supported!");
      calc_N_ = DRT::INPUT::IntegralValue<int>(mfslist, "CALC_N");
      N_vel_ = mfslist.get<double>("N");
      if (mfslist.get<std::string>("REF_VELOCITY") == "strainrate")
        refvel_ = INPAR::FLUID::strainrate;
      else if (mfslist.get<std::string>("REF_VELOCITY") == "resolved")
        refvel_ = INPAR::FLUID::resolved;
      else if (mfslist.get<std::string>("REF_VELOCITY") == "fine_scale")
        refvel_ = INPAR::FLUID::fine_scale;
      else
        dserror("Unknown velocity!");
      if (mfslist.get<std::string>("REF_LENGTH") == "cube_edge")
        reflength_ = INPAR::FLUID::cube_edge;
      else if (mfslist.get<std::string>("REF_LENGTH") == "sphere_diameter")
        reflength_ = INPAR::FLUID::sphere_diameter;
      else if (mfslist.get<std::string>("REF_LENGTH") == "streamlength")
        reflength_ = INPAR::FLUID::streamlength;
      else if (mfslist.get<std::string>("REF_LENGTH") == "gradient_based")
        reflength_ = INPAR::FLUID::gradient_based;
      else if (mfslist.get<std::string>("REF_LENGTH") == "metric_tensor")
        reflength_ = INPAR::FLUID::metric_tensor;
      else
        dserror("Unknown length!");
      c_nu_ = mfslist.get<double>("C_NU");
      nwl_ = DRT::INPUT::IntegralValue<int>(mfslist, "NEAR_WALL_LIMIT");
      // necessary parameters for subgrid-scale scalar estimation
      Csgs_sgphi_ = mfslist.get<double>("CSGS_PHI");
      c_diff_ = mfslist.get<double>("C_DIFF");
      nwl_scatra_ = DRT::INPUT::IntegralValue<int>(mfslist, "NEAR_WALL_LIMIT_CSGS_PHI");
      // general parameters
      beta_ = mfslist.get<double>("BETA");
      if (beta_ != 0.0) dserror("Lhs terms for mfs not included! Fixed-point iteration only!");
      if (mfslist.get<std::string>("EVALUATION_B") == "element_center")
        BD_gp_ = false;
      else if (mfslist.get<std::string>("EVALUATION_B") == "integration_point")
        BD_gp_ = true;
      else
        dserror("Unknown evaluation point!");
      if (mfslist.get<std::string>("CONVFORM") == "convective")
        mfs_conservative_ = false;
      else if (mfslist.get<std::string>("CONVFORM") == "conservative")
        mfs_conservative_ = true;
      else
        dserror("Unknown form of convective term!");

      adapt_Csgs_phi_ = DRT::INPUT::IntegralValue<bool>(mfslist, "ADAPT_CSGS_PHI");
    }
  }

  turbinflow_ = parameters.get<bool>("turbulent inflow");

  return;
}
