// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_turbulence_statistics_tgv.hpp"

#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_utils.hpp"
#include "4C_global_data.hpp"

#include <fstream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------
                  Standard Constructor (public)
  ---------------------------------------------------------------------*/
FLD::TurbulenceStatisticsTgv::TurbulenceStatisticsTgv(
    std::shared_ptr<Core::FE::Discretization> actdis, Teuchos::ParameterList& params,
    const std::string& statistics_outfilename)
    : discret_(actdis), params_(params), statistics_outfilename_(statistics_outfilename)
{
  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim != 3) FOUR_C_THROW("Evaluation of turbulence statistics only for 3d flow!");

  // get number of elements
  numele_ = discret_->num_global_elements();
  // and time-step size
  dt_ = params_.get<double>("time step size");

  // ---------------------------------------------------------------------
  // compute all planes for sampling

  // available planes: everything is written into one variable
  nodeplanes_ = std::make_shared<std::vector<double>>();
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
  std::shared_ptr<std::vector<double>> local_incrvol =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrhk =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrhbazilevs =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrstrle =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrgradle =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);

  std::shared_ptr<std::vector<double>> local_incrtauC =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrtauM =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);

  std::shared_ptr<std::vector<double>> local_incrmk =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);

  std::shared_ptr<std::vector<double>> local_incrres =
      std::make_shared<std::vector<double>>(3 * (nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrres_sq =
      std::make_shared<std::vector<double>>(3 * (nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrabsres =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);

  std::shared_ptr<std::vector<double>> local_incrtauinvsvel =
      std::make_shared<std::vector<double>>(3 * (nodeplanes_->size() - 1), 0.0);

  std::shared_ptr<std::vector<double>> local_incrsvelaf =
      std::make_shared<std::vector<double>>(3 * (nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrsvelaf_sq =
      std::make_shared<std::vector<double>>(3 * (nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrabssvelaf =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);

  std::shared_ptr<std::vector<double>> local_incrresC =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrresC_sq =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrspressnp =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrspressnp_sq =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);

  std::shared_ptr<std::vector<double>> local_incr_eps_pspg =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incr_eps_supg =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incr_eps_cross =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incr_eps_rey =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incr_eps_graddiv =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incr_eps_eddyvisc =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incr_eps_visc =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incr_eps_conv =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incr_eps_mfs =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incr_eps_mfscross =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incr_eps_mfsrey =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incr_eps_avm3 =
      std::make_shared<std::vector<double>>((nodeplanes_->size() - 1), 0.0);

  std::shared_ptr<std::vector<double>> local_incrcrossstress =
      std::make_shared<std::vector<double>>(6 * (nodeplanes_->size() - 1), 0.0);
  std::shared_ptr<std::vector<double>> local_incrreystress =
      std::make_shared<std::vector<double>>(6 * (nodeplanes_->size() - 1), 0.0);

  // pass pointers to local sum vectors to the element
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrvol", local_incrvol);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrhk", local_incrhk);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrhbazilevs", local_incrhbazilevs);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrstrle", local_incrstrle);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrgradle", local_incrgradle);

  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrmk", local_incrmk);

  eleparams_.set<std::shared_ptr<std::vector<double>>>("planecoords_", nodeplanes_);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrtauC", local_incrtauC);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrtauM", local_incrtauM);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrres", local_incrres);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrres_sq", local_incrres_sq);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrabsres", local_incrabsres);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrtauinvsvel", local_incrtauinvsvel);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrsvelaf", local_incrsvelaf);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrsvelaf_sq", local_incrsvelaf_sq);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrabssvelaf", local_incrabssvelaf);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrresC", local_incrresC);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrresC_sq", local_incrresC_sq);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrspressnp", local_incrspressnp);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrspressnp_sq", local_incrspressnp_sq);

  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_pspg", local_incr_eps_pspg);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_supg", local_incr_eps_supg);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_cross", local_incr_eps_cross);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_rey", local_incr_eps_rey);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_graddiv", local_incr_eps_graddiv);
  eleparams_.set<std::shared_ptr<std::vector<double>>>(
      "incr_eps_eddyvisc", local_incr_eps_eddyvisc);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_visc", local_incr_eps_visc);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_conv", local_incr_eps_conv);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_mfs", local_incr_eps_mfs);
  eleparams_.set<std::shared_ptr<std::vector<double>>>(
      "incr_eps_mfscross", local_incr_eps_mfscross);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_mfsrey", local_incr_eps_mfsrey);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_avm3", local_incr_eps_avm3);

  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrcrossstress", local_incrcrossstress);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrreystress", local_incrreystress);

  // means for comparison of of residual and subscale acceleration

  sumres_ = std::make_shared<std::vector<double>>();
  sumres_->resize(3 * (nodeplanes_->size() - 1), 0.0);
  sumres_sq_ = std::make_shared<std::vector<double>>();
  sumres_sq_->resize(3 * (nodeplanes_->size() - 1), 0.0);
  sumabsres_ = std::make_shared<std::vector<double>>();
  sumabsres_->resize((nodeplanes_->size() - 1), 0.0);
  sumtauinvsvel_ = std::make_shared<std::vector<double>>();
  sumtauinvsvel_->resize(3 * (nodeplanes_->size() - 1), 0.0);

  sumsvelaf_ = std::make_shared<std::vector<double>>();
  sumsvelaf_->resize(3 * (nodeplanes_->size() - 1), 0.0);
  sumsvelaf_sq_ = std::make_shared<std::vector<double>>();
  sumsvelaf_sq_->resize(3 * (nodeplanes_->size() - 1), 0.0);
  sumabssvelaf_ = std::make_shared<std::vector<double>>();
  sumabssvelaf_->resize((nodeplanes_->size() - 1), 0.0);

  sumres_c_ = std::make_shared<std::vector<double>>();
  sumres_c_->resize(nodeplanes_->size() - 1, 0.0);
  sumres_c_sq_ = std::make_shared<std::vector<double>>();
  sumres_c_sq_->resize(nodeplanes_->size() - 1, 0.0);

  sumspressnp_ = std::make_shared<std::vector<double>>();
  sumspressnp_->resize(nodeplanes_->size() - 1, 0.0);
  sumspressnp_sq_ = std::make_shared<std::vector<double>>();
  sumspressnp_sq_->resize(nodeplanes_->size() - 1, 0.0);

  sumhk_ = std::make_shared<std::vector<double>>();
  sumhk_->resize(nodeplanes_->size() - 1, 0.0);
  sumhbazilevs_ = std::make_shared<std::vector<double>>();
  sumhbazilevs_->resize(nodeplanes_->size() - 1, 0.0);
  sumstrle_ = std::make_shared<std::vector<double>>();
  sumstrle_->resize(nodeplanes_->size() - 1, 0.0);
  sumgradle_ = std::make_shared<std::vector<double>>();
  sumgradle_->resize(nodeplanes_->size() - 1, 0.0);
  sumtau_m_ = std::make_shared<std::vector<double>>();
  sumtau_m_->resize(nodeplanes_->size() - 1, 0.0);
  sumtau_c_ = std::make_shared<std::vector<double>>();
  sumtau_c_->resize(nodeplanes_->size() - 1, 0.0);

  summk_ = std::make_shared<std::vector<double>>();
  summk_->resize(nodeplanes_->size() - 1, 0.0);

  sum_eps_pspg_ = std::make_shared<std::vector<double>>();
  sum_eps_pspg_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_supg_ = std::make_shared<std::vector<double>>();
  sum_eps_supg_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_cross_ = std::make_shared<std::vector<double>>();
  sum_eps_cross_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_rey_ = std::make_shared<std::vector<double>>();
  sum_eps_rey_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_graddiv_ = std::make_shared<std::vector<double>>();
  sum_eps_graddiv_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_eddyvisc_ = std::make_shared<std::vector<double>>();
  sum_eps_eddyvisc_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_visc_ = std::make_shared<std::vector<double>>();
  sum_eps_visc_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_conv_ = std::make_shared<std::vector<double>>();
  sum_eps_conv_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_mfs_ = std::make_shared<std::vector<double>>();
  sum_eps_mfs_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_mfscross_ = std::make_shared<std::vector<double>>();
  sum_eps_mfscross_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_mfsrey_ = std::make_shared<std::vector<double>>();
  sum_eps_mfsrey_->resize(nodeplanes_->size() - 1, 0.0);
  sum_eps_avm3_ = std::make_shared<std::vector<double>>();
  sum_eps_avm3_->resize(nodeplanes_->size() - 1, 0.0);

  sum_crossstress_ = std::make_shared<std::vector<double>>();
  sum_crossstress_->resize(6 * (nodeplanes_->size() - 1), 0.0);
  sum_reystress_ = std::make_shared<std::vector<double>>();
  sum_reystress_->resize(6 * (nodeplanes_->size() - 1), 0.0);

  // output of residuals and subscale quantities
  std::string s_res(statistics_outfilename_);
  s_res.append(".dissipation");
  std::shared_ptr<std::ofstream> log_res;
  log_res = std::make_shared<std::ofstream>(s_res.c_str(), std::ios::out);
  (*log_res) << "# Statistics for Taylor-Green vortex (residuals and subscale quantities)\n";
  (*log_res) << "# All values are first averaged over the integration points in an element \n";
  (*log_res) << "# and after that averaged over the whole box\n\n";

  // clear statistics
  this->clear_statistics();

  return;
}  // TurbulenceStatisticsCha::TurbulenceStatisticsTgv



/*----------------------------------------------------------------------*
----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsTgv::evaluate_residuals(
    std::map<std::string, std::shared_ptr<Core::LinAlg::Vector<double>>> statevecs,
    std::map<std::string, std::shared_ptr<Core::LinAlg::MultiVector<double>>> statetenss,
    const double thermpressaf, const double thermpressam, const double thermpressdtaf,
    const double thermpressdtam)
{
  //--------------------------------------------------------------------
  // set parameter list (time integration)

  // action for elements
  eleparams_.set<FLD::Action>("action", FLD::calc_dissipation);

  // add velafgrad
  Teuchos::ParameterList* stabparams = &(params_.sublist("RESIDUAL-BASED STABILIZATION"));
  if (stabparams->get<bool>("Reconstruct_Sec_Der"))
  {
    // add velafgrad
    Teuchos::ParameterList* stabparams = &(params_.sublist("RESIDUAL-BASED STABILIZATION"));
    if (stabparams->get<bool>("Reconstruct_Sec_Der"))
    {
      for (std::map<std::string, std::shared_ptr<Core::LinAlg::Vector<double>>>::iterator state =
               statevecs.begin();
           state != statevecs.end(); ++state)
      {
        if (state->first == "velaf")
        {
          FLD::Utils::project_gradient_and_set_param(
              *discret_, eleparams_, state->second, "velafgrad", false);
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
  for (std::map<std::string, std::shared_ptr<Core::LinAlg::Vector<double>>>::iterator state =
           statevecs.begin();
       state != statevecs.end(); ++state)
  {
    discret_->set_state(state->first, state->second);
  }

  // call loop over elements to compute means
  discret_->evaluate(eleparams_, nullptr, nullptr, nullptr, nullptr, nullptr);

  discret_->clear_state();

  // ------------------------------------------------
  // get results from element call via parameter list
  std::shared_ptr<std::vector<double>> local_vol =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrvol");

  std::shared_ptr<std::vector<double>> local_incrhk =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrhk");
  std::shared_ptr<std::vector<double>> local_incrhbazilevs =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrhbazilevs");
  std::shared_ptr<std::vector<double>> local_incrstrle =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrstrle");
  std::shared_ptr<std::vector<double>> local_incrgradle =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrgradle");

  std::shared_ptr<std::vector<double>> local_incrtauC =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrtauC");
  std::shared_ptr<std::vector<double>> local_incrtauM =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrtauM");

  std::shared_ptr<std::vector<double>> local_incrmk =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrmk");

  std::shared_ptr<std::vector<double>> local_incrres =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrres");
  std::shared_ptr<std::vector<double>> local_incrres_sq =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrres_sq");
  std::shared_ptr<std::vector<double>> local_incrabsres =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrabsres");
  std::shared_ptr<std::vector<double>> local_incrtauinvsvel =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrtauinvsvel");
  std::shared_ptr<std::vector<double>> local_incrsvelaf =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrsvelaf");
  std::shared_ptr<std::vector<double>> local_incrsvelaf_sq =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrsvelaf_sq");
  std::shared_ptr<std::vector<double>> local_incrabssvelaf =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrabssvelaf");
  std::shared_ptr<std::vector<double>> local_incrresC =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrresC");
  std::shared_ptr<std::vector<double>> local_incrresC_sq =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrresC_sq");
  std::shared_ptr<std::vector<double>> local_incrspressnp =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrspressnp");
  std::shared_ptr<std::vector<double>> local_incrspressnp_sq =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrspressnp_sq");

  std::shared_ptr<std::vector<double>> local_incr_eps_visc =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incr_eps_visc");
  std::shared_ptr<std::vector<double>> local_incr_eps_conv =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incr_eps_conv");
  std::shared_ptr<std::vector<double>> local_incr_eps_eddyvisc =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incr_eps_eddyvisc");
  std::shared_ptr<std::vector<double>> local_incr_eps_avm3 =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incr_eps_avm3");
  std::shared_ptr<std::vector<double>> local_incr_eps_mfs =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incr_eps_mfs");
  std::shared_ptr<std::vector<double>> local_incr_eps_mfscross =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incr_eps_mfscross");
  std::shared_ptr<std::vector<double>> local_incr_eps_mfsrey =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incr_eps_mfsrey");
  std::shared_ptr<std::vector<double>> local_incr_eps_supg =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incr_eps_supg");
  std::shared_ptr<std::vector<double>> local_incr_eps_cross =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incr_eps_cross");
  std::shared_ptr<std::vector<double>> local_incr_eps_rey =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incr_eps_rey");
  std::shared_ptr<std::vector<double>> local_incr_eps_graddiv =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incr_eps_graddiv");
  std::shared_ptr<std::vector<double>> local_incr_eps_pspg =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incr_eps_pspg");

  std::shared_ptr<std::vector<double>> local_incrcrossstress =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrcrossstress");
  std::shared_ptr<std::vector<double>> local_incrreystress =
      eleparams_.get<std::shared_ptr<std::vector<double>>>("incrreystress");

  int presize = local_incrresC->size();
  int velsize = local_incrres->size();
  int stresssize = local_incrcrossstress->size();

  //--------------------------------------------------
  // vectors to sum over all procs

  // volume of element layers
  std::shared_ptr<std::vector<double>> global_vol;
  global_vol = std::make_shared<std::vector<double>>(presize, 0.0);

  // element sizes of element layers
  std::shared_ptr<std::vector<double>> global_incrhk;
  global_incrhk = std::make_shared<std::vector<double>>(presize, 0.0);

  // element sizes in Bazilevs parameter, viscous regime in element layers
  std::shared_ptr<std::vector<double>> global_incrhbazilevs;
  global_incrhbazilevs = std::make_shared<std::vector<double>>(presize, 0.0);

  // element sizes of element stream length
  std::shared_ptr<std::vector<double>> global_incrstrle;
  global_incrstrle = std::make_shared<std::vector<double>>(presize, 0.0);

  // element sizes based on gradient length
  std::shared_ptr<std::vector<double>> global_incrgradle;
  global_incrgradle = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of tauM/tauC

  std::shared_ptr<std::vector<double>> global_incrtauM;
  global_incrtauM = std::make_shared<std::vector<double>>(presize, 0.0);

  std::shared_ptr<std::vector<double>> global_incrtauC;
  global_incrtauC = std::make_shared<std::vector<double>>(presize, 0.0);

  // mk of element layers
  std::shared_ptr<std::vector<double>> global_incrmk;
  global_incrmk = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of resM (^2) (abs)

  std::shared_ptr<std::vector<double>> global_incrres;
  global_incrres = std::make_shared<std::vector<double>>(velsize, 0.0);

  std::shared_ptr<std::vector<double>> global_incrres_sq;
  global_incrres_sq = std::make_shared<std::vector<double>>(velsize, 0.0);

  std::shared_ptr<std::vector<double>> global_incrtauinvsvel;
  global_incrtauinvsvel = std::make_shared<std::vector<double>>(velsize, 0.0);

  std::shared_ptr<std::vector<double>> global_incrabsres;
  global_incrabsres = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of svelaf (^2) (abs)

  std::shared_ptr<std::vector<double>> global_incrsvelaf;
  global_incrsvelaf = std::make_shared<std::vector<double>>(velsize, 0.0);

  std::shared_ptr<std::vector<double>> global_incrsvelaf_sq;
  global_incrsvelaf_sq = std::make_shared<std::vector<double>>(velsize, 0.0);

  std::shared_ptr<std::vector<double>> global_incrabssvelaf;
  global_incrabssvelaf = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of resC (^2)

  std::shared_ptr<std::vector<double>> global_incrresC;
  global_incrresC = std::make_shared<std::vector<double>>(presize, 0.0);

  std::shared_ptr<std::vector<double>> global_incrresC_sq;
  global_incrresC_sq = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of spressnp (^2)

  std::shared_ptr<std::vector<double>> global_incrspressnp;
  global_incrspressnp = std::make_shared<std::vector<double>>(presize, 0.0);

  std::shared_ptr<std::vector<double>> global_incrspressnp_sq;
  global_incrspressnp_sq = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of dissipation by pspg term

  std::shared_ptr<std::vector<double>> global_incr_eps_pspg;
  global_incr_eps_pspg = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of dissipation by supg term

  std::shared_ptr<std::vector<double>> global_incr_eps_supg;
  global_incr_eps_supg = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of dissipation by cross term

  std::shared_ptr<std::vector<double>> global_incr_eps_cross;
  global_incr_eps_cross = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of dissipation by reynolds term

  std::shared_ptr<std::vector<double>> global_incr_eps_rey;
  global_incr_eps_rey = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of dissipation by continuity stabilisation

  std::shared_ptr<std::vector<double>> global_incr_eps_graddiv;
  global_incr_eps_graddiv = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of dissipation by eddy viscosity

  std::shared_ptr<std::vector<double>> global_incr_eps_eddyvisc;
  global_incr_eps_eddyvisc = std::make_shared<std::vector<double>>(presize, 0.0);

  std::shared_ptr<std::vector<double>> global_incr_eps_avm3;
  global_incr_eps_avm3 = std::make_shared<std::vector<double>>(presize, 0.0);
  // (in plane) averaged values of dissipation by mfs model
  std::shared_ptr<std::vector<double>> global_incr_eps_mfs;
  global_incr_eps_mfs = std::make_shared<std::vector<double>>(presize, 0.0);
  std::shared_ptr<std::vector<double>> global_incr_eps_mfscross;
  global_incr_eps_mfscross = std::make_shared<std::vector<double>>(presize, 0.0);
  std::shared_ptr<std::vector<double>> global_incr_eps_mfsrey;
  global_incr_eps_mfsrey = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of dissipation by viscous forces

  std::shared_ptr<std::vector<double>> global_incr_eps_visc;
  global_incr_eps_visc = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of dissipation/production by convection

  std::shared_ptr<std::vector<double>> global_incr_eps_conv;
  global_incr_eps_conv = std::make_shared<std::vector<double>>(presize, 0.0);

  // (in plane) averaged values of subgrid stresses resulting from supg and cross term

  std::shared_ptr<std::vector<double>> global_incrcrossstress;
  global_incrcrossstress = std::make_shared<std::vector<double>>(stresssize, 0.0);

  // (in plane) averaged values of subgrid stresses resulting from reynolds stresses
  std::shared_ptr<std::vector<double>> global_incrreystress;
  global_incrreystress = std::make_shared<std::vector<double>>(stresssize, 0.0);


  //--------------------------------------------------
  // global sums

  // compute global sum, volume
  discret_->get_comm().SumAll(local_vol->data(), global_vol->data(), presize);

  // compute global sum, element sizes
  discret_->get_comm().SumAll(local_incrhk->data(), global_incrhk->data(), presize);

  // compute global sum, element sizes in viscous regime, Bazilevs parameter
  discret_->get_comm().SumAll(local_incrhbazilevs->data(), global_incrhbazilevs->data(), presize);

  // compute global sum, element sizes
  discret_->get_comm().SumAll(local_incrstrle->data(), global_incrstrle->data(), presize);

  // compute global sum, gradient based element sizes
  discret_->get_comm().SumAll(local_incrgradle->data(), global_incrgradle->data(), presize);

  // compute global sums, stabilisation parameters
  discret_->get_comm().SumAll(local_incrtauM->data(), global_incrtauM->data(), presize);
  discret_->get_comm().SumAll(local_incrtauC->data(), global_incrtauC->data(), presize);

  // compute global sum, mk
  discret_->get_comm().SumAll(local_incrmk->data(), global_incrmk->data(), presize);

  // compute global sums, momentum equation residuals
  discret_->get_comm().SumAll(local_incrres->data(), global_incrres->data(), velsize);
  discret_->get_comm().SumAll(local_incrres_sq->data(), global_incrres_sq->data(), velsize);
  discret_->get_comm().SumAll(local_incrtauinvsvel->data(), global_incrtauinvsvel->data(), velsize);
  discret_->get_comm().SumAll(local_incrabsres->data(), global_incrabsres->data(), presize);

  discret_->get_comm().SumAll(local_incrsvelaf->data(), global_incrsvelaf->data(), velsize);
  discret_->get_comm().SumAll(local_incrsvelaf_sq->data(), global_incrsvelaf_sq->data(), velsize);
  discret_->get_comm().SumAll(local_incrabssvelaf->data(), global_incrabssvelaf->data(), presize);

  // compute global sums, incompressibility residuals
  discret_->get_comm().SumAll(local_incrresC->data(), global_incrresC->data(), presize);
  discret_->get_comm().SumAll(local_incrresC_sq->data(), global_incrresC_sq->data(), presize);

  discret_->get_comm().SumAll(local_incrspressnp->data(), global_incrspressnp->data(), presize);
  discret_->get_comm().SumAll(
      local_incrspressnp_sq->data(), global_incrspressnp_sq->data(), presize);

  // compute global sums, dissipation rates

  discret_->get_comm().SumAll(local_incr_eps_pspg->data(), global_incr_eps_pspg->data(), presize);
  discret_->get_comm().SumAll(local_incr_eps_supg->data(), global_incr_eps_supg->data(), presize);
  discret_->get_comm().SumAll(local_incr_eps_cross->data(), global_incr_eps_cross->data(), presize);
  discret_->get_comm().SumAll(local_incr_eps_rey->data(), global_incr_eps_rey->data(), presize);
  discret_->get_comm().SumAll(
      local_incr_eps_graddiv->data(), global_incr_eps_graddiv->data(), presize);
  discret_->get_comm().SumAll(
      local_incr_eps_eddyvisc->data(), global_incr_eps_eddyvisc->data(), presize);
  discret_->get_comm().SumAll(local_incr_eps_visc->data(), global_incr_eps_visc->data(), presize);
  discret_->get_comm().SumAll(local_incr_eps_conv->data(), global_incr_eps_conv->data(), presize);
  discret_->get_comm().SumAll(local_incr_eps_avm3->data(), global_incr_eps_avm3->data(), presize);
  discret_->get_comm().SumAll(local_incr_eps_mfs->data(), global_incr_eps_mfs->data(), presize);
  discret_->get_comm().SumAll(
      local_incr_eps_mfscross->data(), global_incr_eps_mfscross->data(), presize);
  discret_->get_comm().SumAll(
      local_incr_eps_mfsrey->data(), global_incr_eps_mfsrey->data(), presize);

  // compute global sums, subgrid stresses
  discret_->get_comm().SumAll(
      local_incrcrossstress->data(), global_incrcrossstress->data(), stresssize);
  discret_->get_comm().SumAll(
      local_incrreystress->data(), global_incrreystress->data(), stresssize);


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

    (*sumtau_m_)[rr] += (*global_incrtauM)[rr];
    (*sumtau_c_)[rr] += (*global_incrtauC)[rr];

    (*summk_)[rr] += (*global_incrmk)[rr];

    (*sumres_c_)[rr] += (*global_incrresC)[rr];
    (*sumres_c_sq_)[rr] += (*global_incrresC_sq)[rr];
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
  local_vol = std::make_shared<std::vector<double>>(presize, 0.0);

  local_incrhk = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incrhbazilevs = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incrstrle = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incrgradle = std::make_shared<std::vector<double>>(presize, 0.0);

  local_incrtauC = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incrtauM = std::make_shared<std::vector<double>>(presize, 0.0);

  local_incrmk = std::make_shared<std::vector<double>>(presize, 0.0);

  local_incrres = std::make_shared<std::vector<double>>(velsize, 0.0);
  local_incrres_sq = std::make_shared<std::vector<double>>(velsize, 0.0);
  local_incrsvelaf = std::make_shared<std::vector<double>>(velsize, 0.0);
  local_incrsvelaf_sq = std::make_shared<std::vector<double>>(velsize, 0.0);
  local_incrtauinvsvel = std::make_shared<std::vector<double>>(velsize, 0.0);

  local_incrabsres = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incrabssvelaf = std::make_shared<std::vector<double>>(presize, 0.0);

  local_incrresC = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incrresC_sq = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incrspressnp = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incrspressnp_sq = std::make_shared<std::vector<double>>(presize, 0.0);

  local_incr_eps_pspg = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incr_eps_supg = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incr_eps_cross = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incr_eps_rey = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incr_eps_graddiv = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incr_eps_eddyvisc = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incr_eps_visc = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incr_eps_conv = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incr_eps_avm3 = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incr_eps_mfs = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incr_eps_mfscross = std::make_shared<std::vector<double>>(presize, 0.0);
  local_incr_eps_mfsrey = std::make_shared<std::vector<double>>(presize, 0.0);

  local_incrcrossstress = std::make_shared<std::vector<double>>(stresssize, 0.0);
  local_incrreystress = std::make_shared<std::vector<double>>(stresssize, 0.0);

  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrvol", local_vol);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrhk", local_incrhk);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrhbazilevs", local_incrhbazilevs);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrstrle", local_incrstrle);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrgradle", local_incrgradle);

  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrtauC", local_incrtauC);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrtauM", local_incrtauM);

  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrmk", local_incrmk);

  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrres", local_incrres);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrres_sq", local_incrres_sq);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrabsres", local_incrabsres);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrtauinvsvel", local_incrtauinvsvel);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrsvelaf", local_incrsvelaf);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrsvelaf_sq", local_incrsvelaf_sq);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrabssvelaf", local_incrabssvelaf);

  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrresC", local_incrresC);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrresC_sq", local_incrresC_sq);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrspressnp", local_incrspressnp);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrspressnp_sq", local_incrspressnp_sq);

  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_pspg", local_incr_eps_pspg);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_supg", local_incr_eps_supg);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_cross", local_incr_eps_cross);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_rey", local_incr_eps_rey);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_graddiv", local_incr_eps_graddiv);
  eleparams_.set<std::shared_ptr<std::vector<double>>>(
      "incr_eps_eddyvisc", local_incr_eps_eddyvisc);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_visc", local_incr_eps_visc);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_conv", local_incr_eps_conv);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_avm3", local_incr_eps_avm3);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_mfs", local_incr_eps_mfs);
  eleparams_.set<std::shared_ptr<std::vector<double>>>(
      "incr_eps_mfscross", local_incr_eps_mfscross);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incr_eps_mfsrey", local_incr_eps_mfsrey);

  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrcrossstress", local_incrcrossstress);
  eleparams_.set<std::shared_ptr<std::vector<double>>>("incrreystress", local_incrreystress);

  return;
}  // FLD::TurbulenceStatisticsTgv::EvaluateResiduals



/*----------------------------------------------------------------------*
       Dump the result to file.
  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsTgv::dump_statistics(const int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  std::shared_ptr<std::ofstream> log;
  if (Core::Communication::my_mpi_rank(discret_->get_comm()) == 0)
  {
    std::shared_ptr<std::ofstream> log_res;

    // output of residuals and subscale quantities
    std::string s_res(statistics_outfilename_);
    s_res.append(".dissipation");

    log_res = std::make_shared<std::ofstream>(s_res.c_str(), std::ios::app);

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

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumres_c_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumspressnp_)[rr] / numele_ << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumres_c_sq_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumspressnp_sq_)[rr] / numele_
                 << "  ";

      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumtau_m_)[rr] / numele_ << "  ";
      (*log_res) << std::setw(11) << std::setprecision(4) << (*sumtau_c_)[rr] / numele_ << "  ";

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

}  // TurbulenceStatisticsTvg::time_average_means_and_output_of_statistics


/*----------------------------------------------------------------------*
                  Reset sums and number of samples to 0
  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsTgv::clear_statistics()
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
  for (unsigned rr = 0; rr < sumres_c_->size(); ++rr)
  {
    (*sumabsres_)[rr] = 0.0;
    (*sumabssvelaf_)[rr] = 0.0;

    (*sumhk_)[rr] = 0.0;
    (*sumhbazilevs_)[rr] = 0.0;
    (*sumstrle_)[rr] = 0.0;
    (*sumgradle_)[rr] = 0.0;

    (*sumtau_m_)[rr] = 0.0;
    (*sumtau_c_)[rr] = 0.0;

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

    (*sumres_c_)[rr] = 0.0;
    (*sumspressnp_)[rr] = 0.0;

    (*sumres_c_sq_)[rr] = 0.0;
    (*sumspressnp_sq_)[rr] = 0.0;
  }

  return;
}  // TurbulenceStatisticsTvg::ClearStatistics

FOUR_C_NAMESPACE_CLOSE
