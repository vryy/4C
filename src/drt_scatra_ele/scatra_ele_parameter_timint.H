/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_timint.H

\brief singleton class holding all static time integration parameters required for scalar transport
element evaluation

This singleton class holds all static time integration parameters required for scalar transport
element evaluation. All parameters are usually set only once at the beginning of a simulation,
namely during initialization of the global time integrator, and then never touched again throughout
the simulation. This parameter class needs to coexist with the general parameter class holding all
general static parameters required for scalar transport element evaluation.

<pre>
\maintainer Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15237

\level 1
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_PARAMETER_TIMINT_H
#define SCATRA_ELE_PARAMETER_TIMINT_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "scatra_ele_parameter_base.H"

namespace DRT
{
  namespace ELEMENTS
  {
    /// Evaluation of general parameters (constant over time)
    class ScaTraEleParameterTimInt : public ScaTraEleParameterBase
    {
     public:
      //! singleton access method
      static ScaTraEleParameterTimInt* Instance(
          const std::string& disname,                       //!< name of discretization
          const ScaTraEleParameterTimInt* delete_me = NULL  //!< creation/destruction indication
      );

      //! Singleton destruction method called from global problem
      void Done();

      //! dtor
      ~ScaTraEleParameterTimInt() {}

      //! set parameters
      void SetParameters(Teuchos::ParameterList& parameters  //!< parameter list
      );

      //! print parameters to screen
      void PrintScaTraTimeParameter();

      bool IsGenAlpha() const { return is_genalpha_; };
      bool IsStationary() const { return is_stationary_; };
      bool IsIncremental() const { return is_incremental_; };
      double Time() const { return time_; };
      double Dt() const { return dt_; };
      double TimeFac() const { return timefac_; };
      double TimeFacRhs() const { return timefacrhs_; };
      double TimeFacRhsTau() const { return timefacrhstau_; };
      double AlphaF() const { return alphaF_; };

     private:
      //! private constructor for singletons
      ScaTraEleParameterTimInt();

      bool is_genalpha_;
      bool is_stationary_;
      bool is_incremental_;
      double time_;
      double dt_;
      double timefac_;
      double timefacrhs_;
      double timefacrhstau_;
      double alphaF_;
    };  // class ScaTraEleParameterTimInt
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
