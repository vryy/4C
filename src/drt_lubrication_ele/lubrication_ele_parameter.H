/*--------------------------------------------------------------------------*/
/*!
\file lubrication_ele_parameter.H

\brief singleton class holding all static parameters required for Lubrication element evaluation

This singleton class holds all static parameters required for Lubrication element evaluation. All
parameters are usually set only once at the beginning of a simulation, namely during initialization
of the global time integrator, and then never touched again throughout the simulation. This
parameter class needs to coexist with the general parameter class holding all general static
parameters required for Lubrication element evaluation.

\level 3

\maintainer Mostafa Faraji

*/
/*--------------------------------------------------------------------------*/

#ifndef SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_PARAMETER_H_
#define SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_PARAMETER_H_

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_singletondestruction.H"

namespace DRT
{
  namespace ELEMENTS
  {
    /// Evaluation of general parameters (constant over time)
    class LubricationEleParameter : public DRT::SingletonDestruction
    {
     public:
      //! singleton access method
      static LubricationEleParameter* Instance(
          const std::string& disname,                      //!< name of discretization
          const LubricationEleParameter* delete_me = NULL  //!< creation/destruction indication
      );

      //! Singleton destruction method called from global problem
      void Done();

      //! dtor
      ~LubricationEleParameter() {}

      //! set parameters
      void SetTimeParameters(Teuchos::ParameterList& parameters  //!< parameter list
      );

      //! set parameters
      void SetGeneralParameters(Teuchos::ParameterList& parameters  //!< parameter list
      );

      double Time() const { return time_; };

      /// return function for the modified reynolds equation flag
      bool ModifiedReynolds() const { return modified_reynolds_; };
      /// return function for the Add squeeze term in reynolds equ. flag
      bool AddSqz() const { return addsqz_; };
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
      /// surface roughness STD value in modified reynolds equ.
      double roughness_deviation_;

    };  // class LubricationEleParameter
  }     // namespace ELEMENTS
}  // namespace DRT



#endif /* SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_PARAMETER_H_ */
