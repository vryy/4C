/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluation of general fluid parameter for porofluid

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_FLUID_ELE_PARAMETER_PORO_H
#define BACI_FLUID_ELE_PARAMETER_PORO_H

#include "baci_fluid_ele_parameter.H"
#include "baci_inpar_poroelast.H"
#include "baci_utils_singleton_owner.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidEleParameterPoro : public FluidEleParameter
    {
     public:
      ~FluidEleParameterPoro() override {}

      /// Singleton access method
      static FluidEleParameterPoro* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      void SetElementPoroParameter(Teuchos::ParameterList& params, int myrank);

      //! print parameter to screen
      void PrintFluidParameterPoro() const;

      //! partial integration of porosity gradient in continuity equation
      bool PoroContiPartInt() const { return poro_conti_partint_; };

      //! biot stabilization
      bool StabBiot() const { return stab_biot_; };

      //! add convective term
      bool ConvectiveTerm() const { return poro_convective_term_; };

      //! scaling factor for biot stabilization
      double StabBiotScaling() const { return stab_biot_scaling_; };

      //! flag for inclusion of transient terms in continuity equation
      bool IsStationaryConti() const
      {
        return (not(transient_terms_ == INPAR::POROELAST::transient_all or
                    transient_terms_ == INPAR::POROELAST::transient_continuity_only));
      };

      //! flag for inclusion of transient terms in momentum equation
      bool IsStationaryMomentum() const
      {
        return (not(transient_terms_ == INPAR::POROELAST::transient_all or
                    transient_terms_ == INPAR::POROELAST::transient_momentum_only));
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
      INPAR::POROELAST::TransientEquationsOfPoroFluid transient_terms_;

      /// private Constructor since we are a Singleton.
      FluidEleParameterPoro();
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
