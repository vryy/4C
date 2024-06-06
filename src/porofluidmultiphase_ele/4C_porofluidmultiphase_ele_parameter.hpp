/*----------------------------------------------------------------------*/
/*! \file
 \brief container class holding parameters for element evaluation (singleton)

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_ELE_PARAMETER_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_ELE_PARAMETER_HPP


#include "4C_config.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    /// Evaluation of general parameters (constant over time)
    class PoroFluidMultiPhaseEleParameter
    {
     public:
      //! singleton access method
      static PoroFluidMultiPhaseEleParameter* Instance(
          const std::string& disname  //!< name of discretization
      );

      //! set parameters
      void set_time_step_parameters(Teuchos::ParameterList& parameters  //!< parameter list
      );

      //! set parameters
      void set_general_parameters(Teuchos::ParameterList& parameters  //!< parameter list
      );

      //! @name access methods
      double Time() const { return time_; };
      bool IsGenAlpha() const { return is_genalpha_; };
      bool IsStationary() const { return is_stationary_; };
      double Dt() const { return dt_; };
      double TimeFac() const { return timefac_; };
      double TimeFacRhs() const { return timefacrhs_; };
      double TimeFacRhsTau() const { return timefacrhstau_; };
      double AlphaF() const { return alpha_f_; };
      bool IsAle() const { return is_ale_; };
      bool BiotStab() const { return stab_biot_; };
      int NdsDisp() const { return nds_disp_; };
      int NdsVel() const { return nds_vel_; };
      int NdsSolidPressure() const { return nds_solidpressure_; };
      int NdsScalar() const { return nds_scalar_; };
      bool HasScalar() const { return nds_scalar_ > -1; };
      int num_domain_int_functions() const { return domainint_funct_.size(); };
      std::vector<int> DomainIntFunctions() const { return domainint_funct_; };
      //@}

     private:
      //! private constructor for singletons
      PoroFluidMultiPhaseEleParameter(const std::string& disname  //!< name of discretization
      );

      //! @name parameters potentially changing every time step

      //! current total time
      double time_;
      //! current time step
      double dt_;
      //! time integration factor for left hand side
      double timefac_;
      //! time integration factor for right hand side
      double timefacrhs_;
      //! (time integration factor for right hand side (* (stabilization parameter)
      double timefacrhstau_;
      //! alpha_f parameter from generalized alpha time integration
      double alpha_f_;

      //@}

      //! @name (almost) constant parameters over simulation time

      //! generalized-alpha flag
      bool is_genalpha_;
      //! instationary flag
      bool is_stationary_;
      //! ALE flag
      bool is_ale_;
      //! flag for biot stabilization
      bool stab_biot_;
      //! number of dof set related to mesh displacements
      int nds_disp_;
      //! number of dof set related to mesh velocities
      int nds_vel_;
      //! number of dof set related to solid pressure
      int nds_solidpressure_;
      //! number of dof set related to scalar field
      int nds_scalar_;
      //! setup flag
      bool isset_generalparams_;
      //! domain integral functions
      std::vector<int> domainint_funct_;
      //@}

    };  // class PoroFluidMultiPhaseEleParameter
  }     // namespace ELEMENTS
}  // namespace Discret



FOUR_C_NAMESPACE_CLOSE

#endif
