/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluation of general fluid parameter for porofluid

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_PARAMETER_PORO_HPP
#define FOUR_C_FLUID_ELE_PARAMETER_PORO_HPP

#include "4C_config.hpp"

#include "4C_fluid_ele_parameter.hpp"
#include "4C_inpar_poroelast.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    class FluidEleParameterPoro : public FluidEleParameter
    {
     public:
      /// Singleton access method
      static FluidEleParameterPoro* instance(
          Core::UTILS::SingletonAction action = Core::UTILS::SingletonAction::create);

      void set_element_poro_parameter(Teuchos::ParameterList& params, int myrank);

      //! print parameter to screen
      void print_fluid_parameter_poro() const;

      //! partial integration of porosity gradient in continuity equation
      bool poro_conti_part_int() const { return poro_conti_partint_; };

      //! biot stabilization
      bool stab_biot() const { return stab_biot_; };

      //! add convective term
      bool convective_term() const { return poro_convective_term_; };

      //! scaling factor for biot stabilization
      double stab_biot_scaling() const { return stab_biot_scaling_; };

      //! flag for inclusion of transient terms in continuity equation
      bool is_stationary_conti() const
      {
        return (not(transient_terms_ == Inpar::PoroElast::transient_all or
                    transient_terms_ == Inpar::PoroElast::transient_continuity_only));
      };

      //! flag for inclusion of transient terms in momentum equation
      bool is_stationary_momentum() const
      {
        return (not(transient_terms_ == Inpar::PoroElast::transient_all or
                    transient_terms_ == Inpar::PoroElast::transient_momentum_only));
      };


     private:
      //! Flag SetGeneralParameter was called
      bool set_fluid_parameter_poro_;

      //! partial integration of porosity gradient in continuity equation
      bool poro_conti_partint_;

      //! Flag for biot stabilization
      bool stab_biot_;

      //! scaling factor for biot stabilization
      double stab_biot_scaling_;

      //! additional convective term
      bool poro_convective_term_;

      //! type of handling transient terms
      Inpar::PoroElast::TransientEquationsOfPoroFluid transient_terms_;

      /// private Constructor since we are a Singleton.
      FluidEleParameterPoro();
    };

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
