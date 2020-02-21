/*----------------------------------------------------------------------*/
/*! \file

\brief Compute dissipation rates for Taylor-Green vortex
       and write them to files.

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/
#include <fstream>

#include "turbulence_statistics_tgv.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fluid/fluid_utils.H"

#include "../drt_fluid_ele/fluid_ele_action.H"

/*----------------------------------------------------------------------
                  Standard Constructor (public)
  ---------------------------------------------------------------------*/
FLD::TurbulenceStatisticsTgv::TurbulenceStatisticsTgv(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::ParameterList& params, const std::string& statistics_outfilename)
    : discret_(actdis), params_(params), statistics_outfilename_(statistics_outfilename)
{
  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim != 3) dserror("Evaluation of turbulence statistics only for 3d flow!");

  // get number of elements
  numele_ = discret_->NumGlobalElements();
  // and time-step size
  dt_ = params_.get<double>("time step size");

  // ---------------------------------------------------------------------
  // compute all planes for sampling

  // available planes: everything is written into one variable
  nodeplanes_ = Teuchos::rcp(new std::vector<double>);
  // -10 and -10 chosen such that planes lie outside of 2Pi-box
  nodeplanes_->resize(2);
  (*nodeplanes_)[0] = -10.0;
  (*nodeplanes_)[1] = +10.0;

  //----------------------------------------------------------------------
  // arrays for averaging of residual, subscales etc.

  // prepare time averaging for subscales and residual

  //--------------------------------------------------
  // local_incrtauC            (in plane) averaged values of stabilization parameter tauC
  // local_incrtauM            (in plane) averaged values of stabilization parameter tauM
  // local_incrres(_sq,abs)    (in plane) averaged values of resM (^2) (||.||)
  // local_incrsacc(_sq,abs)   (in plane) averaged values of sacc (^2) (||.||)
  // local_incrsvelaf(_sq,abs) (in plane) averaged values of svelaf (^2) (||.||)
  // local_incrresC(_sq)       (in plane) averaged values of resC (^2)
  // local_incrspressnp(_sq)   (in plane) averaged values of spressnp (^2)
  //--------------------------------------------------
  Teuchos::RCP<std::vector<double>> local_incrvol =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrhk =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrhbazilevs =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrstrle =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrgradle =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

  Teuchos::RCP<std::vector<double>> local_incrtauC =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrtauM =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

  Teuchos::RCP<std::vector<double>> local_incrmk =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

  Teuchos::RCP<std::vector<double>> local_incrres =
      Teuchos::rcp(new std::vector<double>(3 * (nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrres_sq =
      Teuchos::rcp(new std::vector<double>(3 * (nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrabsres =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

  Teuchos::RCP<std::vector<double>> local_incrtauinvsvel =
      Teuchos::rcp(new std::vector<double>(3 * (nodeplanes_->size() - 1), 0.0));

  Teuchos::RCP<std::vector<double>> local_incrsvelaf =
      Teuchos::rcp(new std::vector<double>(3 * (nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrsvelaf_sq =
      Teuchos::rcp(new std::vector<double>(3 * (nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrabssvelaf =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

  Teuchos::RCP<std::vector<double>> local_incrresC =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrresC_sq =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrspressnp =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrspressnp_sq =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

  Teuchos::RCP<std::vector<double>> local_incr_eps_pspg =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incr_eps_supg =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incr_eps_cross =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incr_eps_rey =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incr_eps_graddiv =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incr_eps_eddyvisc =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incr_eps_visc =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incr_eps_conv =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incr_eps_mfs =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incr_eps_mfscross =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incr_eps_mfsrey =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incr_eps_avm3 =
      Teuchos::rcp(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

  Teuchos::RCP<std::vector<double>> local_incrcrossstress =
      Teuchos::rcp(new std::vector<double>(6 * (nodeplanes_->size() - 1), 0.0));
  Teuchos::RCP<std::vector<double>> local_incrreystress =
      Teuchos::rcp(new std::vector<double>(6 * (nodeplanes_->size() - 1), 0.0));

  // pass pointers to local sum vectors to the element
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrvol", local_incrvol);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrhk", local_incrhk);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrhbazilevs", local_incrhbazilevs);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrstrle", local_incrstrle);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrgradle", local_incrgradle);

  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrmk", local_incrmk);

  eleparams_.set<Teuchos::RCP<std::vector<double>>>("planecoords_", nodeplanes_);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauC", local_incrtauC);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauM", local_incrtauM);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrres", local_incrres);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrres_sq", local_incrres_sq);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrabsres", local_incrabsres);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauinvsvel", local_incrtauinvsvel);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrsvelaf", local_incrsvelaf);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrsvelaf_sq", local_incrsvelaf_sq);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrabssvelaf", local_incrabssvelaf);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrresC", local_incrresC);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrresC_sq", local_incrresC_sq);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrspressnp", local_incrspressnp);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrspressnp_sq", local_incrspressnp_sq);

  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_pspg", local_incr_eps_pspg);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_supg", local_incr_eps_supg);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_cross", local_incr_eps_cross);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_rey", local_incr_eps_rey);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_graddiv", local_incr_eps_graddiv);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_eddyvisc", local_incr_eps_eddyvisc);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_visc", local_incr_eps_visc);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_conv", local_incr_eps_conv);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_mfs", local_incr_eps_mfs);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_mfscross", local_incr_eps_mfscross);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_mfsrey", local_incr_eps_mfsrey);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_avm3", local_incr_eps_avm3);

  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrcrossstress", local_incrcrossstress);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrreystress", local_incrreystress);

  // means for comparison of of residual and subscale acceleration

  sumres_ = Teuchos::rcp(new std::vector<double>);
  sumres_->resize(3 * (nodeplanes_->size() - 1), 0.0);
  sumres_sq_ = Teuchos::rcp(new std::vector<double>);
  sumres_sq_->resize(3 * (nodeplanes_->size() - 1), 0.0);
  sumabsres_ = Teuchos::rcp(new std::vector<double>);
  sumabsres_->resize((nodeplanes_->size() - 1), 0.0);
  sumtauinvsvel_ = Teuchos::rcp(new std::vector<double>);
  sumtauinvsvel_->resize(3 * (nodeplanes_->size() - 1), 0.0);

  sumsvelaf_ = Teuchos::rcp(new std::vector<double>);
  sumsvelaf_->resize(3 * (nodeplanes_->size() - 1), 0.0);
  sumsvelaf_sq_ = Teuchos::rcp(new std::vector<double>);
  sumsvelaf_sq_->resize(3 * (nodeplanes_->size() - 1), 0.0);
  sumabssvelaf_ = Teuchos::rcp(new std::vector<double>);
  sumabssvelaf_->resize((nodeplanes_->size() - 1), 0.0);

  sumresC_ = Teuchos::rcp(new std::vector<double>);
  sumresC_->resize(nodeplanes_->size() - 1, 0.0);
  sumresC_sq_ = Teuchos::rcp(new std::vector<double>);
  sumresC_sq_->resize(nodeplanes_->size() - 1, 0.0);

  sumspressnp_ = Teuchos::rcp(new std::vector<double>);
  sumspressnp_->resize(nodeplanes_->size() - 1, 0.0);
  sumspressnp_sq_ = Teuchos::rcp(new std::vector<double>);
  sumspressnp_sq_->resize(nodeplanes_->size() - 1, 0.0);

  sumhk_ = Teuchos::rcp(new std::vector<double>);
  sumhk_->resize(nodeplanes_->size() - 1, 0.0);
  sumhbazilevs_ = Teuchos::rcp(new std::vector<double>);
  sumhbazilevs_->resize(nodeplanes_->size() - 1, 0.0);
  sumstrle_ = Teuchos::rcp(new std::vector<double>);
  sumstrle_->resize(nodeplanes_->size() - 1, 0.0);
  sumgradle_ = Teuchos::rcp(new std::vector<double>);
  sumgradle_->resize(nodeplanes_->size() - 1, 0.0);
  sumtauM_ = Teuchos::rcp(new std::vector<double>);
  sumtauM_->resize(nodeplanes_->size() - 1, 0.0);
  sumtauC_ = Teuchos::rcp(new std::vector<double>);
  sumtauC_->resize(nodeplanes_->size() - 1, 0.0);

  summk_ = Teuchos::rcp(new std::vector<double>);
  summk_->resize(nodeplanes_->size() - 1, 0.0);

  sum_eps_pspg_ = Teuchos::rcp(new std::vector<double>);
  sum_eps_pspg_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_supg_ = Teuchos::rcp(new std::vector<double>);
  sum_eps_supg_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_cross_ = Teuchos::rcp(new std::vector<double>);
  sum_eps_cross_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_rey_ = Teuchos::rcp(new std::vector<double>);
  sum_eps_rey_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_graddiv_ = Teuchos::rcp(new std::vector<double>);
  sum_eps_graddiv_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_eddyvisc_ = Teuchos::rcp(new std::vector<double>);
  sum_eps_eddyvisc_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_visc_ = Teuchos::rcp(new std::vector<double>);
  sum_eps_visc_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_conv_ = Teuchos::rcp(new std::vector<double>);
  sum_eps_conv_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_mfs_ = Teuchos::rcp(new std::vector<double>);
  sum_eps_mfs_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_mfscross_ = Teuchos::rcp(new std::vector<double>);
  sum_eps_mfscross_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_mfsrey_ = Teuchos::rcp(new std::vector<double>);
  sum_eps_mfsrey_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_avm3_ = Teuchos::rcp(new std::vector<double>);
  sum_eps_avm3_->resize(nodeplanes_->size() - 1, 0.0);

  sum_crossstress_ = Teuchos::rcp(new std::vector<double>);
  sum_crossstress_->resize(6 * (nodeplanes_->size() - 1), 0.0);
  sum_reystress_ = Teuchos::rcp(new std::vector<double>);
  sum_reystress_->resize(6 * (nodeplanes_->size() - 1), 0.0);

  // output of residuals and subscale quantities
  std::string s_res(statistics_outfilename_);
  s_res.append(".dissipation");
  Teuchos::RCP<std::ofstream> log_res;
  log_res = Teuchos::rcp(new std::ofstream(s_res.c_str(), std::ios::out));
  (*log_res) << "# Statistics for Taylor-Green vortex (residuals and subscale quantities)\n";
  (*log_res) << "# All values are first averaged over the integration points in an element \n";
  (*log_res) << "# and after that averaged over the whole box\n\n";

  // clear statistics
  this->ClearStatistics();

  return;
}  // TurbulenceStatisticsCha::TurbulenceStatisticsTgv


/*----------------------------------------------------------------------*
                           Destructor
 -----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsTgv::~TurbulenceStatisticsTgv()
{
  return;
}  // TurbulenceStatisticsCha::~TurbulenceStatisticsTgv()


/*----------------------------------------------------------------------*
----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsTgv::EvaluateResiduals(
    std::map<std::string, Teuchos::RCP<Epetra_Vector>> statevecs,
    std::map<std::string, Teuchos::RCP<Epetra_MultiVector>> statetenss, const double thermpressaf,
    const double thermpressam, const double thermpressdtaf, const double thermpressdtam)
{
  //--------------------------------------------------------------------
  // set parameter list (time integration)

  // action for elements
  eleparams_.set<int>("action", FLD::calc_dissipation);

  // add velafgrad
  Teuchos::ParameterList* stabparams = &(params_.sublist("RESIDUAL-BASED STABILIZATION"));
  if (DRT::INPUT::IntegralValue<int>(*stabparams, "Reconstruct_Sec_Der"))
  {
    // add velafgrad
    Teuchos::ParameterList* stabparams = &(params_.sublist("RESIDUAL-BASED STABILIZATION"));
    if (DRT::INPUT::IntegralValue<int>(*stabparams, "Reconstruct_Sec_Der"))
    {
      for (std::map<std::string, Teuchos::RCP<Epetra_Vector>>::iterator state = statevecs.begin();
           state != statevecs.end(); ++state)
      {
        if (state->first == "velaf")
        {
          FLD::UTILS::ProjectGradientAndSetParam(
              discret_, eleparams_, state->second, "velafgrad", false);
          break;
        }
      }
    }
  }

  eleparams_.set<double>("thermpress at n+alpha_F/n+1", thermpressaf);
  eleparams_.set<double>("thermpress at n+alpha_M/n", thermpressam);
  eleparams_.set<double>("thermpressderiv at n+alpha_F/n+1", thermpressdtaf);
  eleparams_.set<double>("thermpressderiv at n+alpha_M/n+1", thermpressdtam);

  // set state vectors for element call
  for (std::map<std::string, Teuchos::RCP<Epetra_Vector>>::iterator state = statevecs.begin();
       state != statevecs.end(); ++state)
  {
    discret_->SetState(state->first, state->second);
  }

  // call loop over elements to compute means
  discret_->Evaluate(
      eleparams_, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  discret_->ClearState();

  // ------------------------------------------------
  // get results from element call via parameter list
  Teuchos::RCP<std::vector<double>> local_vol =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrvol");

  Teuchos::RCP<std::vector<double>> local_incrhk =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrhk");
  Teuchos::RCP<std::vector<double>> local_incrhbazilevs =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrhbazilevs");
  Teuchos::RCP<std::vector<double>> local_incrstrle =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrstrle");
  Teuchos::RCP<std::vector<double>> local_incrgradle =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrgradle");

  Teuchos::RCP<std::vector<double>> local_incrtauC =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrtauC");
  Teuchos::RCP<std::vector<double>> local_incrtauM =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrtauM");

  Teuchos::RCP<std::vector<double>> local_incrmk =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrmk");

  Teuchos::RCP<std::vector<double>> local_incrres =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrres");
  Teuchos::RCP<std::vector<double>> local_incrres_sq =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrres_sq");
  Teuchos::RCP<std::vector<double>> local_incrabsres =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrabsres");
  Teuchos::RCP<std::vector<double>> local_incrtauinvsvel =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrtauinvsvel");
  Teuchos::RCP<std::vector<double>> local_incrsvelaf =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrsvelaf");
  Teuchos::RCP<std::vector<double>> local_incrsvelaf_sq =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrsvelaf_sq");
  Teuchos::RCP<std::vector<double>> local_incrabssvelaf =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrabssvelaf");
  Teuchos::RCP<std::vector<double>> local_incrresC =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrresC");
  Teuchos::RCP<std::vector<double>> local_incrresC_sq =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrresC_sq");
  Teuchos::RCP<std::vector<double>> local_incrspressnp =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrspressnp");
  Teuchos::RCP<std::vector<double>> local_incrspressnp_sq =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrspressnp_sq");

  Teuchos::RCP<std::vector<double>> local_incr_eps_visc =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_visc");
  Teuchos::RCP<std::vector<double>> local_incr_eps_conv =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_conv");
  Teuchos::RCP<std::vector<double>> local_incr_eps_eddyvisc =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_eddyvisc");
  Teuchos::RCP<std::vector<double>> local_incr_eps_avm3 =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_avm3");
  Teuchos::RCP<std::vector<double>> local_incr_eps_mfs =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_mfs");
  Teuchos::RCP<std::vector<double>> local_incr_eps_mfscross =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_mfscross");
  Teuchos::RCP<std::vector<double>> local_incr_eps_mfsrey =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_mfsrey");
  Teuchos::RCP<std::vector<double>> local_incr_eps_supg =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_supg");
  Teuchos::RCP<std::vector<double>> local_incr_eps_cross =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_cross");
  Teuchos::RCP<std::vector<double>> local_incr_eps_rey =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_rey");
  Teuchos::RCP<std::vector<double>> local_incr_eps_graddiv =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_graddiv");
  Teuchos::RCP<std::vector<double>> local_incr_eps_pspg =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_pspg");

  Teuchos::RCP<std::vector<double>> local_incrcrossstress =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrcrossstress");
  Teuchos::RCP<std::vector<double>> local_incrreystress =
      eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrreystress");

  int presize = local_incrresC->size();
  int velsize = local_incrres->size();
  int stresssize = local_incrcrossstress->size();

  //--------------------------------------------------
  // vectors to sum over all procs

  // volume of element layers
  Teuchos::RCP<std::vector<double>> global_vol;
  global_vol = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // element sizes of element layers
  Teuchos::RCP<std::vector<double>> global_incrhk;
  global_incrhk = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // element sizes in Bazilevs parameter, viscous regime in element layers
  Teuchos::RCP<std::vector<double>> global_incrhbazilevs;
  global_incrhbazilevs = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // element sizes of element stream length
  Teuchos::RCP<std::vector<double>> global_incrstrle;
  global_incrstrle = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // element sizes based on gradient length
  Teuchos::RCP<std::vector<double>> global_incrgradle;
  global_incrgradle = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of tauM/tauC

  Teuchos::RCP<std::vector<double>> global_incrtauM;
  global_incrtauM = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  Teuchos::RCP<std::vector<double>> global_incrtauC;
  global_incrtauC = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // mk of element layers
  Teuchos::RCP<std::vector<double>> global_incrmk;
  global_incrmk = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of resM (^2) (abs)

  Teuchos::RCP<std::vector<double>> global_incrres;
  global_incrres = Teuchos::rcp(new std::vector<double>(velsize, 0.0));

  Teuchos::RCP<std::vector<double>> global_incrres_sq;
  global_incrres_sq = Teuchos::rcp(new std::vector<double>(velsize, 0.0));

  Teuchos::RCP<std::vector<double>> global_incrtauinvsvel;
  global_incrtauinvsvel = Teuchos::rcp(new std::vector<double>(velsize, 0.0));

  Teuchos::RCP<std::vector<double>> global_incrabsres;
  global_incrabsres = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of svelaf (^2) (abs)

  Teuchos::RCP<std::vector<double>> global_incrsvelaf;
  global_incrsvelaf = Teuchos::rcp(new std::vector<double>(velsize, 0.0));

  Teuchos::RCP<std::vector<double>> global_incrsvelaf_sq;
  global_incrsvelaf_sq = Teuchos::rcp(new std::vector<double>(velsize, 0.0));

  Teuchos::RCP<std::vector<double>> global_incrabssvelaf;
  global_incrabssvelaf = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of resC (^2)

  Teuchos::RCP<std::vector<double>> global_incrresC;
  global_incrresC = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  Teuchos::RCP<std::vector<double>> global_incrresC_sq;
  global_incrresC_sq = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of spressnp (^2)

  Teuchos::RCP<std::vector<double>> global_incrspressnp;
  global_incrspressnp = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  Teuchos::RCP<std::vector<double>> global_incrspressnp_sq;
  global_incrspressnp_sq = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of dissipation by pspg term

  Teuchos::RCP<std::vector<double>> global_incr_eps_pspg;
  global_incr_eps_pspg = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of dissipation by supg term

  Teuchos::RCP<std::vector<double>> global_incr_eps_supg;
  global_incr_eps_supg = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of dissipation by cross term

  Teuchos::RCP<std::vector<double>> global_incr_eps_cross;
  global_incr_eps_cross = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of dissipation by reynolds term

  Teuchos::RCP<std::vector<double>> global_incr_eps_rey;
  global_incr_eps_rey = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of dissipation by continuity stabilisation

  Teuchos::RCP<std::vector<double>> global_incr_eps_graddiv;
  global_incr_eps_graddiv = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of dissipation by eddy viscosity

  Teuchos::RCP<std::vector<double>> global_incr_eps_eddyvisc;
  global_incr_eps_eddyvisc = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  Teuchos::RCP<std::vector<double>> global_incr_eps_avm3;
  global_incr_eps_avm3 = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  // (in plane) averaged values of dissipation by mfs model
  Teuchos::RCP<std::vector<double>> global_incr_eps_mfs;
  global_incr_eps_mfs = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  Teuchos::RCP<std::vector<double>> global_incr_eps_mfscross;
  global_incr_eps_mfscross = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  Teuchos::RCP<std::vector<double>> global_incr_eps_mfsrey;
  global_incr_eps_mfsrey = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of dissipation by viscous forces

  Teuchos::RCP<std::vector<double>> global_incr_eps_visc;
  global_incr_eps_visc = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of dissipation/production by convection

  Teuchos::RCP<std::vector<double>> global_incr_eps_conv;
  global_incr_eps_conv = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  // (in plane) averaged values of subgrid stresses resulting from supg and cross term

  Teuchos::RCP<std::vector<double>> global_incrcrossstress;
  global_incrcrossstress = Teuchos::rcp(new std::vector<double>(stresssize, 0.0));

  // (in plane) averaged values of subgrid stresses resulting from reynolds stresses
  Teuchos::RCP<std::vector<double>> global_incrreystress;
  global_incrreystress = Teuchos::rcp(new std::vector<double>(stresssize, 0.0));


  //--------------------------------------------------
  // global sums

  // compute global sum, volume
  discret_->Comm().SumAll(&((*local_vol)[0]), &((*global_vol)[0]), presize);

  // compute global sum, element sizes
  discret_->Comm().SumAll(&((*local_incrhk)[0]), &((*global_incrhk)[0]), presize);

  // compute global sum, element sizes in viscous regime, Bazilevs parameter
  discret_->Comm().SumAll(&((*local_incrhbazilevs)[0]), &((*global_incrhbazilevs)[0]), presize);

  // compute global sum, element sizes
  discret_->Comm().SumAll(&((*local_incrstrle)[0]), &((*global_incrstrle)[0]), presize);

  // compute global sum, gradient based element sizes
  discret_->Comm().SumAll(&((*local_incrgradle)[0]), &((*global_incrgradle)[0]), presize);

  // compute global sums, stabilisation parameters
  discret_->Comm().SumAll(&((*local_incrtauM)[0]), &((*global_incrtauM)[0]), presize);
  discret_->Comm().SumAll(&((*local_incrtauC)[0]), &((*global_incrtauC)[0]), presize);

  // compute global sum, mk
  discret_->Comm().SumAll(&((*local_incrmk)[0]), &((*global_incrmk)[0]), presize);

  // compute global sums, momentum equation residuals
  discret_->Comm().SumAll(&((*local_incrres)[0]), &((*global_incrres)[0]), velsize);
  discret_->Comm().SumAll(&((*local_incrres_sq)[0]), &((*global_incrres_sq)[0]), velsize);
  discret_->Comm().SumAll(&((*local_incrtauinvsvel)[0]), &((*global_incrtauinvsvel)[0]), velsize);
  discret_->Comm().SumAll(&((*local_incrabsres)[0]), &((*global_incrabsres)[0]), presize);

  discret_->Comm().SumAll(&((*local_incrsvelaf)[0]), &((*global_incrsvelaf)[0]), velsize);
  discret_->Comm().SumAll(&((*local_incrsvelaf_sq)[0]), &((*global_incrsvelaf_sq)[0]), velsize);
  discret_->Comm().SumAll(&((*local_incrabssvelaf)[0]), &((*global_incrabssvelaf)[0]), presize);

  // compute global sums, incompressibility residuals
  discret_->Comm().SumAll(&((*local_incrresC)[0]), &((*global_incrresC)[0]), presize);
  discret_->Comm().SumAll(&((*local_incrresC_sq)[0]), &((*global_incrresC_sq)[0]), presize);

  discret_->Comm().SumAll(&((*local_incrspressnp)[0]), &((*global_incrspressnp)[0]), presize);
  discret_->Comm().SumAll(&((*local_incrspressnp_sq)[0]), &((*global_incrspressnp_sq)[0]), presize);

  // compute global sums, dissipation rates

  discret_->Comm().SumAll(&((*local_incr_eps_pspg)[0]), &((*global_incr_eps_pspg)[0]), presize);
  discret_->Comm().SumAll(&((*local_incr_eps_supg)[0]), &((*global_incr_eps_supg)[0]), presize);
  discret_->Comm().SumAll(&((*local_incr_eps_cross)[0]), &((*global_incr_eps_cross)[0]), presize);
  discret_->Comm().SumAll(&((*local_incr_eps_rey)[0]), &((*global_incr_eps_rey)[0]), presize);
  discret_->Comm().SumAll(
      &((*local_incr_eps_graddiv)[0]), &((*global_incr_eps_graddiv)[0]), presize);
  discret_->Comm().SumAll(
      &((*local_incr_eps_eddyvisc)[0]), &((*global_incr_eps_eddyvisc)[0]), presize);
  discret_->Comm().SumAll(&((*local_incr_eps_visc)[0]), &((*global_incr_eps_visc)[0]), presize);
  discret_->Comm().SumAll(&((*local_incr_eps_conv)[0]), &((*global_incr_eps_conv)[0]), presize);
  discret_->Comm().SumAll(&((*local_incr_eps_avm3)[0]), &((*global_incr_eps_avm3)[0]), presize);
  discret_->Comm().SumAll(&((*local_incr_eps_mfs)[0]), &((*global_incr_eps_mfs)[0]), presize);
  discret_->Comm().SumAll(
      &((*local_incr_eps_mfscross)[0]), &((*global_incr_eps_mfscross)[0]), presize);
  discret_->Comm().SumAll(&((*local_incr_eps_mfsrey)[0]), &((*global_incr_eps_mfsrey)[0]), presize);

  // compute global sums, subgrid stresses
  discret_->Comm().SumAll(
      &((*local_incrcrossstress)[0]), &((*global_incrcrossstress)[0]), stresssize);
  discret_->Comm().SumAll(&((*local_incrreystress)[0]), &((*global_incrreystress)[0]), stresssize);


  for (int rr = 0; rr < velsize; ++rr)
  {
    (*sumres_)[rr] += (*global_incrres)[rr];
    (*sumres_sq_)[rr] += (*global_incrres_sq)[rr];
    (*sumsvelaf_)[rr] += (*global_incrsvelaf)[rr];
    (*sumsvelaf_sq_)[rr] += (*global_incrsvelaf_sq)[rr];

    (*sumtauinvsvel_)[rr] += (*global_incrtauinvsvel)[rr];
  }
  for (int rr = 0; rr < presize; ++rr)
  {
    (*sumabsres_)[rr] += (*global_incrabsres)[rr];
    (*sumabssvelaf_)[rr] += (*global_incrabssvelaf)[rr];

    (*sumhk_)[rr] += (*global_incrhk)[rr];
    (*sumhbazilevs_)[rr] += (*global_incrhbazilevs)[rr];
    (*sumstrle_)[rr] += (*global_incrstrle)[rr];
    (*sumgradle_)[rr] += (*global_incrgradle)[rr];

    (*sumtauM_)[rr] += (*global_incrtauM)[rr];
    (*sumtauC_)[rr] += (*global_incrtauC)[rr];

    (*summk_)[rr] += (*global_incrmk)[rr];

    (*sumresC_)[rr] += (*global_incrresC)[rr];
    (*sumresC_sq_)[rr] += (*global_incrresC_sq)[rr];
    (*sumspressnp_)[rr] += (*global_incrspressnp)[rr];
    (*sumspressnp_sq_)[rr] += (*global_incrspressnp_sq)[rr];

    (*sum_eps_pspg_)[rr] += (*global_incr_eps_pspg)[rr];
    (*sum_eps_supg_)[rr] += (*global_incr_eps_supg)[rr];
    (*sum_eps_cross_)[rr] += (*global_incr_eps_cross)[rr];
    (*sum_eps_rey_)[rr] += (*global_incr_eps_rey)[rr];
    (*sum_eps_graddiv_)[rr] += (*global_incr_eps_graddiv)[rr];
    (*sum_eps_eddyvisc_)[rr] += (*global_incr_eps_eddyvisc)[rr];
    (*sum_eps_visc_)[rr] += (*global_incr_eps_visc)[rr];
    (*sum_eps_conv_)[rr] += (*global_incr_eps_conv)[rr];
    (*sum_eps_avm3_)[rr] += (*global_incr_eps_avm3)[rr];
    (*sum_eps_mfs_)[rr] += (*global_incr_eps_mfs)[rr];
    (*sum_eps_mfscross_)[rr] += (*global_incr_eps_mfscross)[rr];
    (*sum_eps_mfsrey_)[rr] += (*global_incr_eps_mfsrey)[rr];
  }

  for (int rr = 0; rr < stresssize; ++rr)
  {
    (*sum_crossstress_)[rr] += (*global_incrcrossstress)[rr];
    (*sum_reystress_)[rr] += (*global_incrreystress)[rr];
  }

  // reset working arrays
  local_vol = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  local_incrhk = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incrhbazilevs = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incrstrle = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incrgradle = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  local_incrtauC = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incrtauM = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  local_incrmk = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  local_incrres = Teuchos::rcp(new std::vector<double>(velsize, 0.0));
  local_incrres_sq = Teuchos::rcp(new std::vector<double>(velsize, 0.0));
  local_incrsvelaf = Teuchos::rcp(new std::vector<double>(velsize, 0.0));
  local_incrsvelaf_sq = Teuchos::rcp(new std::vector<double>(velsize, 0.0));
  local_incrtauinvsvel = Teuchos::rcp(new std::vector<double>(velsize, 0.0));

  local_incrabsres = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incrabssvelaf = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  local_incrresC = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incrresC_sq = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incrspressnp = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incrspressnp_sq = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  local_incr_eps_pspg = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incr_eps_supg = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incr_eps_cross = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incr_eps_rey = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incr_eps_graddiv = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incr_eps_eddyvisc = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incr_eps_visc = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incr_eps_conv = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incr_eps_avm3 = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incr_eps_mfs = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incr_eps_mfscross = Teuchos::rcp(new std::vector<double>(presize, 0.0));
  local_incr_eps_mfsrey = Teuchos::rcp(new std::vector<double>(presize, 0.0));

  local_incrcrossstress = Teuchos::rcp(new std::vector<double>(stresssize, 0.0));
  local_incrreystress = Teuchos::rcp(new std::vector<double>(stresssize, 0.0));

  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrvol", local_vol);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrhk", local_incrhk);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrhbazilevs", local_incrhbazilevs);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrstrle", local_incrstrle);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrgradle", local_incrgradle);

  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauC", local_incrtauC);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauM", local_incrtauM);

  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrmk", local_incrmk);

  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrres", local_incrres);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrres_sq", local_incrres_sq);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrabsres", local_incrabsres);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauinvsvel", local_incrtauinvsvel);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrsvelaf", local_incrsvelaf);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrsvelaf_sq", local_incrsvelaf_sq);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrabssvelaf", local_incrabssvelaf);

  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrresC", local_incrresC);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrresC_sq", local_incrresC_sq);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrspressnp", local_incrspressnp);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrspressnp_sq", local_incrspressnp_sq);

  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_pspg", local_incr_eps_pspg);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_supg", local_incr_eps_supg);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_cross", local_incr_eps_cross);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_rey", local_incr_eps_rey);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_graddiv", local_incr_eps_graddiv);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_eddyvisc", local_incr_eps_eddyvisc);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_visc", local_incr_eps_visc);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_conv", local_incr_eps_conv);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_avm3", local_incr_eps_avm3);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_mfs", local_incr_eps_mfs);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_mfscross", local_incr_eps_mfscross);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_mfsrey", local_incr_eps_mfsrey);

  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrcrossstress", local_incrcrossstress);
  eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrreystress", local_incrreystress);

  return;
}  // FLD::TurbulenceStatisticsTgv::EvaluateResiduals



/*----------------------------------------------------------------------*
       Dump the result to file.
  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsTgv::DumpStatistics(const int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->Comm().MyPID() == 0)
  {
    Teuchos::RCP<std::ofstream> log_res;

    // output of residuals and subscale quantities
    std::string s_res(statistics_outfilename_);
    s_res.append(".dissipation");

    log_res = Teuchos::rcp(new std::ofstream(s_res.c_str(), std::ios::app));

    if (step == 1)
    {
      (*log_res) << "\n\n\n";

      (*log_res) << "#     time    ";
      (*log_res) << "    step  ";
      (*log_res) << "    res_x   ";
      (*log_res) << "      res_y  ";
      (*log_res) << "      res_z  ";
      (*log_res) << "     svel_x  ";
      (*log_res) << "     svel_y  ";
      (*log_res) << "     svel_z  ";

      (*log_res) << "   res_sq_x  ";
      (*log_res) << "   res_sq_y  ";
      (*log_res) << "   res_sq_z  ";
      (*log_res) << "   svel_sq_x ";
      (*log_res) << "   svel_sq_y ";
      (*log_res) << "   svel_sq_z ";

      (*log_res) << " tauinvsvel_x";
      (*log_res) << " tauinvsvel_y";
      (*log_res) << " tauinvsvel_z";

      (*log_res) << "    ||res||  ";
      (*log_res) << "   ||svel||  ";

      (*log_res) << "      resC   ";
      (*log_res) << "    spresnp  ";

      (*log_res) << "    resC_sq  ";
      (*log_res) << "  spresnp_sq ";

      (*log_res) << "    tauM     ";
      (*log_res) << "    tauC     ";

      (*log_res) << "  eps_pspg   ";
      (*log_res) << "  eps_supg   ";
      (*log_res) << "  eps_cross  ";
      (*log_res) << "   eps_rey   ";
      (*log_res) << "  eps_graddiv  ";
      (*log_res) << " eps_eddyvisc";
      (*log_res) << "   eps_visc  ";
      (*log_res) << "   eps_conv  ";
      (*log_res) << "   eps_avm3  ";
      (*log_res) << "   eps_mfs   ";
      (*log_res) << " eps_mfscross";
      (*log_res) << " eps_mfsrey  ";

      (*log_res) << "     hk      ";
      (*log_res) << "   strle     ";
      (*log_res) << "   gradle    ";
      (*log_res) << " h_bazilevs  ";

      (*log_res) << " tau_cross_11";
      (*log_res) << " tau_cross_22";
      (*log_res) << " tau_cross_33";
      (*log_res) << " tau_cross_12";
      (*log_res) << " tau_cross_23";
      (*log_res) << " tau_cross_31";
      (*log_res) << " tau_rey_11  ";
      (*log_res) << " tau_rey_22  ";
      (*log_res) << " tau_rey_33  ";
      (*log_res) << " tau_rey_12  ";
      (*log_res) << " tau_rey_23  ";
      (*log_res) << " tau_rey_31  ";
      (*log_res) << " mk          ";
      (*log_res) << "\n";
    }

    (*log_res) << std::scientific;
    for (unsigned rr = 0; rr < nodeplanes_->size() - 1; ++rr)
    {
      (*log_res) << std::setw(11) << std::setprecision(4) << step * dt_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << step << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumres_)[3 * rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumres_)[3 * rr + 1] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumres_)[3 * rr + 2] / numele_
                 << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumsvelaf_)[3 * rr] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumsvelaf_)[3 * rr + 1] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumsvelaf_)[3 * rr + 2] / numele_
                 << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumres_sq_)[3 * rr] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumres_sq_)[3 * rr + 1] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumres_sq_)[3 * rr + 2] / numele_
                 << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumsvelaf_sq_)[3 * rr] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumsvelaf_sq_)[3 * rr + 1] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumsvelaf_sq_)[3 * rr + 2] / numele_
                 << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumtauinvsvel_)[3 * rr] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumtauinvsvel_)[3 * rr + 1] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumtauinvsvel_)[3 * rr + 2] / numele_
                 << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumabsres_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumabssvelaf_)[rr] / numele_ << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumresC_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumspressnp_)[rr] / numele_ << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumresC_sq_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumspressnp_sq_)[rr] / numele_
                 << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumtauM_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumtauC_)[rr] / numele_ << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_eps_pspg_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_eps_supg_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_eps_cross_)[rr] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_eps_rey_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_eps_graddiv_)[rr] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_eps_eddyvisc_)[rr] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_eps_visc_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_eps_conv_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_eps_avm3_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_eps_mfs_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_eps_mfscross_)[rr] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_eps_mfsrey_)[rr] / numele_
                 << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumhk_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumstrle_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumgradle_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumhbazilevs_)[rr] / numele_ << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_crossstress_)[6 * rr] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4)
                 << (*sum_crossstress_)[6 * rr + 1] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4)
                 << (*sum_crossstress_)[6 * rr + 2] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4)
                 << (*sum_crossstress_)[6 * rr + 3] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4)
                 << (*sum_crossstress_)[6 * rr + 4] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4)
                 << (*sum_crossstress_)[6 * rr + 5] / numele_ << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_reystress_)[6 * rr] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_reystress_)[6 * rr + 1] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_reystress_)[6 * rr + 2] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_reystress_)[6 * rr + 3] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_reystress_)[6 * rr + 4] / numele_
                 << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sum_reystress_)[6 * rr + 5] / numele_
                 << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*summk_)[rr] / numele_ << "  ";

      (*log_res) << &std::endl;
    }
    log_res->flush();

  }  // end myrank 0


  return;

}  // TurbulenceStatisticsTvg::TimeAverageMeansAndOutputOfStatistics


/*----------------------------------------------------------------------*
                  Reset sums and number of samples to 0
  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsTgv::ClearStatistics()
{
  // reset residuals and subscale averages
  for (unsigned rr = 0; rr < sumres_->size() / 3; ++rr)
  {
    (*sumres_)[3 * rr] = 0.0;
    (*sumres_)[3 * rr + 1] = 0.0;
    (*sumres_)[3 * rr + 2] = 0.0;

    (*sumsvelaf_)[3 * rr] = 0.0;
    (*sumsvelaf_)[3 * rr + 1] = 0.0;
    (*sumsvelaf_)[3 * rr + 2] = 0.0;

    (*sumres_sq_)[3 * rr] = 0.0;
    (*sumres_sq_)[3 * rr + 1] = 0.0;
    (*sumres_sq_)[3 * rr + 2] = 0.0;

    (*sumsvelaf_sq_)[3 * rr] = 0.0;
    (*sumsvelaf_sq_)[3 * rr + 1] = 0.0;
    (*sumsvelaf_sq_)[3 * rr + 2] = 0.0;

    (*sumtauinvsvel_)[3 * rr] = 0.0;
    (*sumtauinvsvel_)[3 * rr + 1] = 0.0;
    (*sumtauinvsvel_)[3 * rr + 2] = 0.0;

    for (int mm = 0; mm < 6; ++mm)
    {
      (*sum_crossstress_)[6 * rr + mm] = 0.0;
      (*sum_reystress_)[6 * rr + mm] = 0.0;
    }
  }
  for (unsigned rr = 0; rr < sumresC_->size(); ++rr)
  {
    (*sumabsres_)[rr] = 0.0;
    (*sumabssvelaf_)[rr] = 0.0;

    (*sumhk_)[rr] = 0.0;
    (*sumhbazilevs_)[rr] = 0.0;
    (*sumstrle_)[rr] = 0.0;
    (*sumgradle_)[rr] = 0.0;

    (*sumtauM_)[rr] = 0.0;
    (*sumtauC_)[rr] = 0.0;

    (*summk_)[rr] = 0.0;

    (*sum_eps_pspg_)[rr] = 0.0;
    (*sum_eps_supg_)[rr] = 0.0;
    (*sum_eps_cross_)[rr] = 0.0;
    (*sum_eps_rey_)[rr] = 0.0;
    (*sum_eps_graddiv_)[rr] = 0.0;
    (*sum_eps_eddyvisc_)[rr] = 0.0;
    (*sum_eps_visc_)[rr] = 0.0;
    (*sum_eps_conv_)[rr] = 0.0;
    (*sum_eps_mfs_)[rr] = 0.0;
    (*sum_eps_mfscross_)[rr] = 0.0;
    (*sum_eps_mfsrey_)[rr] = 0.0;
    (*sum_eps_avm3_)[rr] = 0.0;

    (*sumresC_)[rr] = 0.0;
    (*sumspressnp_)[rr] = 0.0;

    (*sumresC_sq_)[rr] = 0.0;
    (*sumspressnp_sq_)[rr] = 0.0;
  }

  return;
}  // TurbulenceStatisticsTvg::ClearStatistics
