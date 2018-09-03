/*--------------------------------------------------------------------------*/
/*!
\file lubrication_ele_parameter.H

\brief singleton class holding all static parameters required for Lubrication element evaluation

This singleton class holds all static parameters required for Lubrication element evaluation. All
parameters are usually set only once at the beginning of a simulation, namely during initialization
of the global time integrator, and then never touched again throughout the simulation. This
parameter class needs to coexist with the general parameter class holding all general static
parameters required for Lubrication element evaluation.

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
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

     private:
      //! private constructor for singletons
      LubricationEleParameter(const std::string& disname  //!< name of discretization
      );

      double time_;

    };  // class LubricationEleParameter
  }     // namespace ELEMENTS
}  // namespace DRT



#endif /* SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_PARAMETER_H_ */
