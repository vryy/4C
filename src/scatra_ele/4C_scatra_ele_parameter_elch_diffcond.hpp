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
#ifndef FOUR_C_SCATRA_ELE_PARAMETER_ELCH_DIFFCOND_HPP
#define FOUR_C_SCATRA_ELE_PARAMETER_ELCH_DIFFCOND_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_parameter_base.hpp"

#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleParameterStd;

    // class implementation
    class ScaTraEleParameterElchDiffCond : public ScaTraEleParameterBase
    {
     public:
      //! singleton access method
      static ScaTraEleParameterElchDiffCond* instance(
          const std::string& disname  //!< name of discretization
      );

      //! set parameters
      void set_parameters(Teuchos::ParameterList& parameters  //!< parameter list
          ) override;

      //! return flag for current as solution variable
      bool cur_sol_var() const { return cursolvar_; }

      //! return flag for diffusion potential
      bool diffusion_coeff_based() const { return diffusioncoefbased_; }

      //! return Newman constants
      double newman_constdata() const { return newmanconsta_; }
      double newman_const_b() const { return newmanconstb_; }
      double newman_const_c() const { return newmanconstc_; }
      double epsilon0() const { return epsilon_0_; }

     private:
      //! private constructor for singletons
      ScaTraEleParameterElchDiffCond(const std::string& disname  //!< name of discretization
      );

      //! flag if current is used as a solution variable
      bool cursolvar_;

      //! mat_diffcond: flag if diffusion potential is based on diffusion coefficients or
      //! transference number
      bool diffusioncoefbased_;

      //! switch for dilute and concentrated solution theory (diffusion potential in current
      //! equation):
      //!    A          B
      //!   |--|  |----------|
      //!   z_1 + (z_2 - z_1) t_1
      //! ------------------------ (RT/F kappa 1/c_k grad c_k)
      //!      z_1 z_2
      //!     |________|
      //!         C
      double newmanconsta_;
      double newmanconstb_;
      double newmanconstc_;
      double epsilon_0_;
    };
  }  // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
