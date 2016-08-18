/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_elch.H

\brief singleton class holding all static electrochemistry parameters required for element evaluation

This singleton class holds all static electrochemistry parameters required for element evaluation. All parameters are
usually set only once at the beginning of a simulation, namely during initialization of the global time integrator,
and then never touched again throughout the simulation. This parameter class needs to coexist with the general parameter
class holding all general static parameters required for scalar transport element evaluation.

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_PARAMETER_ELCH_H
#define SCATRA_ELE_PARAMETER_ELCH_H

#include "../drt_inpar/inpar_elch.H"

#include "scatra_ele_parameter_base.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    class ScaTraEleParameterElch : public ScaTraEleParameterBase
    {
      public:

        //! destructor
        virtual ~ScaTraEleParameterElch(){return;};

        //! singleton access method
        static ScaTraEleParameterElch* Instance(
            const std::string&   disname,                   //!< name of discretization
            const ScaTraEleParameterElch* delete_me = NULL  //!< creation/destruction indication
            );

        //! singleton destruction
        virtual void Done();

        //! return flag for coupling of lithium-ion flux density and electric current density at Dirichlet and Neumann boundaries
        bool BoundaryFluxCoupling() const {return boundaryfluxcoupling_;};

        //! set parameters
        void SetParameters(
            Teuchos::ParameterList& parameters   //!< parameter list
            );

        //! return type of closing equation for electric potential
        INPAR::ELCH::EquPot EquPot() const {return equpot_;};

        //! return Faraday constant
        double Faraday() const {return faraday_;};

        //! return dielectric constant
        double Epsilon() const {return epsilon_;};

        //! return constant F/RT
        double FRT() const {return frt_;};

      private:

        //! private constructor for singletons
        ScaTraEleParameterElch(
            const std::string& disname   //!< name of discretization
            );

        //! flag for coupling of lithium-ion flux density and electric current density at Dirichlet and Neumann boundaries
        bool boundaryfluxcoupling_;

        //! equation used for closing of the elch-system
        enum INPAR::ELCH::EquPot equpot_;

        //! Faraday constant
        const double faraday_;

        //! dielectric constant
        const double epsilon_;

        //! pre-calculation of regularly used constant F/RT
        //! (a division is much more expensive than a multiplication)
        double frt_;
    };
  } // namespace ELEMENTS
} // namespace DRT
#endif
