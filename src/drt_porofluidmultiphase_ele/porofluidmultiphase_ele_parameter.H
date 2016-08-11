/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_ele_parameter.H

 \brief container class holding parameters for element evaluation (singleton)

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_PARAMETER_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_PARAMETER_H_


#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_singletondestruction.H"

namespace DRT
{
  namespace ELEMENTS
  {
    /// Evaluation of general parameters (constant over time)
    class PoroFluidMultiPhaseEleParameter : public DRT::SingletonDestruction
    {
      public:

        //! singleton access method
        static PoroFluidMultiPhaseEleParameter* Instance(
            const std::string&           disname,   //!< name of discretization
            const PoroFluidMultiPhaseEleParameter* delete_me = NULL //!< creation/destruction indication
            );

        //! Singleton destruction method called from global problem
        void Done();

        //! dtor
        ~PoroFluidMultiPhaseEleParameter() {}

        //! set parameters
        void SetTimeStepParameters(
            Teuchos::ParameterList& parameters   //!< parameter list
            );

        //! set parameters
        void SetGeneralParameters(
            Teuchos::ParameterList& parameters   //!< parameter list
            );

        //! @name access methods
        double Time() const {return time_;};
        bool IsGenAlpha() const {return is_genalpha_;};
        bool IsStationary() const {return is_stationary_;};
        double Dt() const {return dt_;};
        double TimeFac() const {return timefac_;};
        double TimeFacRhs() const {return timefacrhs_;};
        double TimeFacRhsTau() const {return timefacrhstau_;};
        double AlphaF() const {return alphaF_;};
        bool IsAle() const {return is_ale_;};
        int NdsDisp() const {return nds_disp_;};
        int NdsVel() const {return nds_vel_;};
        int NdsSolidPressure() const {return nds_solidpressure_;};
        //@}

      private:

        //! private constructor for singletons
        PoroFluidMultiPhaseEleParameter(
            const std::string& disname   //!< name of discretization
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
        double alphaF_;

        //@}

        //! @name (almost) constant parameters over simulation time

        //! generalized-alpha flag
        bool is_genalpha_;
        //! instationary flag
        bool is_stationary_;
        //! ALE flag
        bool is_ale_;
        //! number of dof set related to mesh displacements
        int nds_disp_;
        //! number of dof set related to mesh velocities
        int nds_vel_;
        //! number of dof set related to solid pressure
        int nds_solidpressure_;
        //! setup flag
        bool isset_generalparams_;
        //@}

    }; // class PoroFluidMultiPhaseEleParameter
  } // namespace ELEMENTS
} // namespace DRT



#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_PARAMETER_H_ */
