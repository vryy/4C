// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_PARAMETER_LSREINIT_HPP
#define FOUR_C_SCATRA_ELE_PARAMETER_LSREINIT_HPP

#include "4C_config.hpp"

#include "4C_inpar_levelset.hpp"
#include "4C_scatra_ele_parameter_base.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    class ScaTraEleParameterLsReinit : public ScaTraEleParameterBase
    {
     protected:
      /// (private) protected constructor, since we are a Singleton.
      /// this constructor is called from a derived class
      /// -> therefore, it has to be protected instead of private
      ScaTraEleParameterLsReinit(const std::string& disname  //!< name of discretization
      );

     public:
      /// Singleton access method
      static ScaTraEleParameterLsReinit* instance(
          const std::string& disname  //!< name of discretization
      );

      /*========================================================================*/
      //! @name set-routines
      /*========================================================================*/

      //! set parameters
      void set_parameters(Teuchos::ParameterList& parameters  //!< parameter list
          ) override;

      /*========================================================================*/
      //! @name access-routines
      /*========================================================================*/

      Inpar::ScaTra::ReInitialAction reinit_type() const { return reinittype_; };
      Inpar::ScaTra::SmoothedSignType sign_type() const { return signtype_; };
      Inpar::ScaTra::CharEleLengthReinit char_ele_length_reinit() const
      {
        return charelelengthreinit_;
      };
      double interface_thickness_fac() const { return interfacethicknessfac_; };
      bool use_projected_vel() const { return useprojectedreinitvel_; };
      Inpar::ScaTra::LinReinit lin_form() const { return linform_; };
      Inpar::ScaTra::ArtDiff art_diff() const { return artdiff_; };

      /// access the penalty parameter for the elliptical reinitialization
      double penalty_para() const
      {
        FOUR_C_ASSERT(alphapen_ > 0.0, "The penalty parameter have to be larger than zero!");
        return alphapen_;
      };

      bool project() const { return project_; };
      bool lumping() const { return lumping_; };
      double project_diff() const { return projectdiff_; };
      Inpar::ScaTra::DiffFunc diff_fct() const { return difffct_; };

     private:
      // reinit type
      Inpar::ScaTra::ReInitialAction reinittype_;
      // sign function for phi
      Inpar::ScaTra::SmoothedSignType signtype_;
      // element length for smoothing of sign function
      Inpar::ScaTra::CharEleLengthReinit charelelengthreinit_;
      // interface thickness factor (multiple of characteristic element length)
      double interfacethicknessfac_;
      // from of velocity evaluation
      bool useprojectedreinitvel_;
      // form of linearization of nonlinear terms
      Inpar::ScaTra::LinReinit linform_;
      // form of artificial diffusion
      Inpar::ScaTra::ArtDiff artdiff_;
      // penalty parameter of elliptic reinitialization
      double alphapen_;
      // use L2 projection
      bool project_;
      // diffusion for L2 projection
      double projectdiff_;
      // mass lumping for L2 projection
      bool lumping_;
      // function for diffusivity
      Inpar::ScaTra::DiffFunc difffct_;
    };
  }  // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
