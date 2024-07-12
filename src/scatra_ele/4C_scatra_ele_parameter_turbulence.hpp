/*----------------------------------------------------------------------*/
/*! \file

\brief singleton class holding all static turbulence parameters required for scalar transport
element evaluation

This singleton class holds all static turbulence parameters required for scalar transport element
evaluation. All parameters are usually set only once at the beginning of a simulation, namely during
initialization of the global time integrator, and then never touched again throughout the
simulation. This parameter class needs to coexist with the general parameter class holding all
general static parameters required for scalar transport element evaluation.


\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_PARAMETER_TURBULENCE_HPP
#define FOUR_C_SCATRA_ELE_PARAMETER_TURBULENCE_HPP

#include "4C_config.hpp"

#include "4C_inpar_fluid.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_inpar_turbulence.hpp"
#include "4C_scatra_ele_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleParameterTimInt;

    // class implementation
    class ScaTraEleParameterTurbulence : public ScaTraEleParameterBase
    {
     public:
      //! singleton access method
      static ScaTraEleParameterTurbulence* instance(
          const std::string& disname  //!< name of discretization
      );



      //! set parameters
      void set_parameters(Teuchos::ParameterList& parameters  //!< parameter list
          ) override;

      //! @name return turbulence parameters
      //! @{
      Inpar::FLUID::TurbModelAction turb_model() { return turbmodel_; };
      Inpar::FLUID::ScalarForcing scalar_forcing() { return scalarforcing_; };
      Inpar::ScaTra::FSSUGRDIFF which_fssgd() { return whichfssgd_; };
      bool fssgd() { return fssgd_; };
      double cs() { return cs_; };
      double tpn() { return tpn_; };
      bool cs_av() { return cs_av_; };
      double csgs_sg_vel() { return csgs_sgvel_; };
      double alpha() { return alpha_; }
      bool calc_n() { return calc_n_; };
      double n_vel() { return n_vel_; };
      Inpar::FLUID::RefVelocity ref_vel() { return refvel_; };
      Inpar::FLUID::RefLength ref_length() { return reflength_; };
      double c_nu() { return c_nu_; };
      bool nwl() { return nwl_; };
      bool nwl_scatra() { return nwl_scatra_; };
      bool beta() { return beta_; };
      bool bd_gp() { return bd_gp_; };
      double csgs_sg_phi()
      {
        double tmp = 0.0;
        if (adapt_csgs_phi_ and nwl_)
          tmp = csgs_sgvel_ * mean_cai_;
        else
          tmp = csgs_sgphi_;
        return tmp;
      };
      double c_diff() { return c_diff_; };
      bool mfs_conservative() { return mfs_conservative_; };
      void set_csgs_phi(double meanCai)
      {
        mean_cai_ = meanCai;
        return;
      };
      bool adapt_csgs_phi() { return adapt_csgs_phi_; };
      bool turb_inflow() { return turbinflow_; };
      //! @}

     private:
      //! private constructor for singletons
      ScaTraEleParameterTurbulence(const std::string& disname  //!< name of discretization
      );

      //! @name turbulence parameters
      //! @{
      //! definition of turbulence model
      Inpar::FLUID::TurbModelAction turbmodel_;

      //! define forcing for scalar field
      Inpar::FLUID::ScalarForcing scalarforcing_;

      //! flag to activate AVM3
      bool fssgd_;

      //! type of AVM3
      Inpar::ScaTra::FSSUGRDIFF whichfssgd_;

      //! parameters for subgrid-diffusivity models
      double cs_;
      double tpn_;
      bool cs_av_;

      //! parameters for multifractal subgrid-scale modeling
      double csgs_sgvel_;
      double alpha_;
      bool calc_n_;
      double n_vel_;
      Inpar::FLUID::RefVelocity refvel_;
      Inpar::FLUID::RefLength reflength_;
      double c_nu_;
      bool nwl_;
      bool nwl_scatra_;
      bool beta_;
      bool bd_gp_;
      double csgs_sgphi_;
      double c_diff_;
      bool mfs_conservative_;
      double mean_cai_;
      bool adapt_csgs_phi_;

      //! further parameter
      bool turbinflow_;
      //! @}

      //! parameter class for time integration
      Discret::ELEMENTS::ScaTraEleParameterTimInt* timintparams_;
    };
  }  // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
