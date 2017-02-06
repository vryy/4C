/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_lsreinit.H

\brief singleton class holding all static levelset reinitialization parameters required for element evaluation

This singleton class holds all static levelset reinitialization parameters required for element evaluation. All
parameters are usually set only once at the beginning of a simulation, namely during initialization of the global
time integrator, and then never touched again throughout the simulation. This parameter class needs to coexist with
the general parameter class holding all general static parameters required for scalar transport element evaluation.

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>

\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_PARAMETER_LSREINIT_H
#define SCATRA_ELE_PARAMETER_LSREINIT_H

#include "scatra_ele_parameter_base.H"

#include "../drt_inpar/inpar_levelset.H"
#include "../drt_lib/drt_dserror.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class ScaTraEleParameterLsReinit : public ScaTraEleParameterBase
    {
      protected:
        /// (private) protected constructor, since we are a Singleton.
        /// this constructor is called from a derived class
        /// -> therefore, it has to be protected instead of private
        ScaTraEleParameterLsReinit(
            const std::string& disname   //!< name of discretization
            );

      public:

        virtual ~ScaTraEleParameterLsReinit() {}

        /// Singleton access method
        static ScaTraEleParameterLsReinit* Instance(
            const std::string&   disname,                       //!< name of discretization
            const ScaTraEleParameterLsReinit* delete_me = NULL  //!< creation/destruction indication
            );

        /// called upon destruction
        virtual void Done();

        /*========================================================================*/
        //! @name set-routines
        /*========================================================================*/

        //! set parameters
        void SetParameters(
            Teuchos::ParameterList& parameters   //!< parameter list
            );

        /*========================================================================*/
        //! @name access-routines
        /*========================================================================*/

        INPAR::SCATRA::ReInitialAction ReinitType() const {return reinittype_;};
        INPAR::SCATRA::SmoothedSignType SignType() const {return signtype_;};
        INPAR::SCATRA::CharEleLengthReinit CharEleLengthReinit() const {return charelelengthreinit_;};
        double InterfaceThicknessFac() const {return interfacethicknessfac_;};
        bool UseProjectedVel() const {return useprojectedreinitvel_;};
        INPAR::SCATRA::LinReinit LinForm() const {return linform_;};
        INPAR::SCATRA::ArtDiff ArtDiff() const {return artdiff_;};

        /// access the penalty parameter for the elliptical reinitialization
        double PenaltyPara() const
        {
          dsassert(alphapen_>0.0, "The penalty parameter have to be larger than zero!");
          return alphapen_;
        };

        bool Project() const {return project_;};
        bool Lumping() const {return lumping_;};
        double ProjectDiff() const {return projectdiff_;};
        INPAR::SCATRA::DiffFunc DiffFct() const {return difffct_;};

      private:

        // reinit type
        INPAR::SCATRA::ReInitialAction reinittype_;
        // sign function for phi
        INPAR::SCATRA::SmoothedSignType signtype_;
        // element length for smoothing of sign function
        INPAR::SCATRA::CharEleLengthReinit charelelengthreinit_;
        // interface thickness factor (multiple of characteristic element length)
        double interfacethicknessfac_;
        // from of velocity evaluation
        bool useprojectedreinitvel_;
        // form of linearization of nonlinear terms
        INPAR::SCATRA::LinReinit linform_;
        // form of artificial diffusion
        INPAR::SCATRA::ArtDiff artdiff_;
        // penalty parameter of elliptic reinitialization
        double alphapen_;
        // use L2 projection
        bool project_;
        // diffusion for L2 projection
        double projectdiff_;
        // mass lumping for L2 projection
        bool lumping_;
        // function for diffusivity
        INPAR::SCATRA::DiffFunc difffct_;
    };
  }
}
#endif
