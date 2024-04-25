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

namespace DRT
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
      static ScaTraEleParameterTurbulence* Instance(
          const std::string& disname  //!< name of discretization
      );



      //! set parameters
      void SetParameters(Teuchos::ParameterList& parameters  //!< parameter list
          ) override;

      //! @name return turbulence parameters
      //! @{
      INPAR::FLUID::TurbModelAction TurbModel() { return turbmodel_; };
      INPAR::FLUID::ScalarForcing ScalarForcing() { return scalarforcing_; };
      INPAR::SCATRA::FSSUGRDIFF WhichFssgd() { return whichfssgd_; };
      bool FSSGD() { return fssgd_; };
      double Cs() { return cs_; };
      double TPN() { return tpn_; };
      bool CsAv() { return cs_av_; };
      double Csgs_SgVel() { return csgs_sgvel_; };
      double Alpha() { return alpha_; }
      bool Calc_N() { return calc_n_; };
      double N_Vel() { return n_vel_; };
      INPAR::FLUID::RefVelocity RefVel() { return refvel_; };
      INPAR::FLUID::RefLength RefLength() { return reflength_; };
      double C_Nu() { return c_nu_; };
      bool Nwl() { return nwl_; };
      bool Nwl_ScaTra() { return nwl_scatra_; };
      bool Beta() { return beta_; };
      bool BD_Gp() { return bd_gp_; };
      double Csgs_SgPhi()
      {
        double tmp = 0.0;
        if (adapt_csgs_phi_ and nwl_)
          tmp = csgs_sgvel_ * mean_cai_;
        else
          tmp = csgs_sgphi_;
        return tmp;
      };
      double C_Diff() { return c_diff_; };
      bool MfsConservative() { return mfs_conservative_; };
      void SetCsgsPhi(double meanCai)
      {
        mean_cai_ = meanCai;
        return;
      };
      bool AdaptCsgsPhi() { return adapt_csgs_phi_; };
      bool TurbInflow() { return turbinflow_; };
      //! @}

     private:
      //! private constructor for singletons
      ScaTraEleParameterTurbulence(const std::string& disname  //!< name of discretization
      );

      //! @name turbulence parameters
      //! @{
      //! definition of turbulence model
      INPAR::FLUID::TurbModelAction turbmodel_;

      //! define forcing for scalar field
      INPAR::FLUID::ScalarForcing scalarforcing_;

      //! flag to activate AVM3
      bool fssgd_;

      //! type of AVM3
      INPAR::SCATRA::FSSUGRDIFF whichfssgd_;

      //! parameters for subgrid-diffusivity models
      double cs_;
      double tpn_;
      bool cs_av_;

      //! parameters for multifractal subgrid-scale modeling
      double csgs_sgvel_;
      double alpha_;
      bool calc_n_;
      double n_vel_;
      INPAR::FLUID::RefVelocity refvel_;
      INPAR::FLUID::RefLength reflength_;
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
      DRT::ELEMENTS::ScaTraEleParameterTimInt* timintparams_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
