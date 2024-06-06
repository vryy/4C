/*--------------------------------------------------------------------------*/
/*! \file

\brief singleton class holding all static parameters required for Lubrication element evaluation

This singleton class holds all static parameters required for Lubrication element evaluation. All
parameters are usually set only once at the beginning of a simulation, namely during initialization
of the global time integrator, and then never touched again throughout the simulation. This
parameter class needs to coexist with the general parameter class holding all general static
parameters required for Lubrication element evaluation.

\level 3


*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_LUBRICATION_ELE_PARAMETER_HPP
#define FOUR_C_LUBRICATION_ELE_PARAMETER_HPP

#include "4C_config.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    /// Evaluation of general parameters (constant over time)
    class LubricationEleParameter
    {
     public:
      //! singleton access method
      static LubricationEleParameter* Instance(
          const std::string& disname  //!< name of discretization
      );

      //! set parameters
      void SetTimeParameters(Teuchos::ParameterList& parameters  //!< parameter list
      );

      //! set parameters
      void set_general_parameters(Teuchos::ParameterList& parameters  //!< parameter list
      );

      double Time() const { return time_; };

      /// return function for the modified reynolds equation flag
      bool ModifiedReynolds() const { return modified_reynolds_; };
      /// return function for the Add squeeze term in reynolds equ. flag
      bool AddSqz() const { return addsqz_; };
      /// return function for the pure lubrication reynolds equ. flag
      bool PureLub() const { return purelub_; };
      /// return function for surface roughness value in modified reynolds equ.
      double RoughnessDeviation() const { return roughness_deviation_; };

     private:
      //! private constructor for singletons
      LubricationEleParameter(const std::string& disname  //!< name of discretization
      );

      double time_;
      /// modified reynolds equation flag
      bool modified_reynolds_;
      /// Add squeeze term in reynolds equ. flag
      bool addsqz_;
      /// pure lubrication reynolds equ. flag
      bool purelub_;
      /// surface roughness STD value in modified reynolds equ.
      double roughness_deviation_;

    };  // class LubricationEleParameter
  }     // namespace ELEMENTS
}  // namespace Discret



FOUR_C_NAMESPACE_CLOSE

#endif
