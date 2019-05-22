/*----------------------------------------------------------------------*/
/*!

\brief singleton class holding all static turbulence parameters required for scalar transport
element evaluation

This singleton class holds all static turbulence parameters required for scalar transport element
evaluation. All parameters are usually set only once at the beginning of a simulation, namely during
initialization of the global time integrator, and then never touched again throughout the
simulation. This parameter class needs to coexist with the general parameter class holding all
general static parameters required for scalar transport element evaluation.

\maintainer Anh-Tu Vuong

\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_PARAMETER_TURBULENCE_H
#define SCATRA_ELE_PARAMETER_TURBULENCE_H

#include "scatra_ele_parameter_base.H"

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_turbulence.H"

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
          const std::string& disname,                           //!< name of discretization
          const ScaTraEleParameterTurbulence* delete_me = NULL  //!< creation/destruction indication
      );

      //! singleton destruction
      virtual void Done();

      //! set parameters
      void SetParameters(Teuchos::ParameterList& parameters  //!< parameter list
      );

      //! @name return turbulence parameters
      //! @{
      INPAR::FLUID::TurbModelAction TurbModel() { return turbmodel_; };
      INPAR::FLUID::ScalarForcing ScalarForcing() { return scalarforcing_; };
      INPAR::SCATRA::FSSUGRDIFF WhichFssgd() { return whichfssgd_; };
      bool FSSGD() { return fssgd_; };
      double Cs() { return Cs_; };
      double TPN() { return tpn_; };
      bool CsAv() { return Cs_av_; };
      double Csgs_SgVel() { return Csgs_sgvel_; };
      double Alpha() { return alpha_; }
      bool Calc_N() { return calc_N_; };
      double N_Vel() { return N_vel_; };
      INPAR::FLUID::RefVelocity RefVel() { return refvel_; };
      INPAR::FLUID::RefLength RefLength() { return reflength_; };
      double C_Nu() { return c_nu_; };
      bool Nwl() { return nwl_; };
      bool Nwl_ScaTra() { return nwl_scatra_; };
      bool Beta() { return beta_; };
      bool BD_Gp() { return BD_gp_; };
      double Csgs_SgPhi()
      {
        double tmp = 0.0;
        if (adapt_Csgs_phi_ and nwl_)
          tmp = Csgs_sgvel_ * meanCai_;
        else
          tmp = Csgs_sgphi_;
        return tmp;
      };
      double C_Diff() { return c_diff_; };
      bool MfsConservative() { return mfs_conservative_; };
      void SetCsgsPhi(double meanCai)
      {
        meanCai_ = meanCai;
        return;
      };
      bool AdaptCsgsPhi() { return adapt_Csgs_phi_; };
      bool TurbInflow() { return turbinflow_; };
      //! @}

     private:
      //! private constructor for singletons
      ScaTraEleParameterTurbulence(const std::string& disname  //!< name of discretization
      );

      //! private destructor for singletons
      virtual ~ScaTraEleParameterTurbulence() { return; };

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
      double Cs_;
      double tpn_;
      bool Cs_av_;

      //! parameters for multifractal subgrid-scale modeling
      double Csgs_sgvel_;
      double alpha_;
      bool calc_N_;
      double N_vel_;
      INPAR::FLUID::RefVelocity refvel_;
      INPAR::FLUID::RefLength reflength_;
      double c_nu_;
      bool nwl_;
      bool nwl_scatra_;
      bool beta_;
      bool BD_gp_;
      double Csgs_sgphi_;
      double c_diff_;
      bool mfs_conservative_;
      double meanCai_;
      bool adapt_Csgs_phi_;

      //! further parameter
      bool turbinflow_;
      //! @}

      //! parameter class for time integration
      DRT::ELEMENTS::ScaTraEleParameterTimInt* timintparams_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
