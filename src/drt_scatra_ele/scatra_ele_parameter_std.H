/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_std.H

\brief singleton class holding all static parameters required for the evaluation of a standard scalar transport element

This singleton class holds all static parameters required for the evaluation of a standard scalar transport element,
e.g., stabilization parameters and finite difference check parameters. All parameters are usually set only once at
the beginning of a simulation, namely during initialization of the global time integrator, and then never touched again
throughout the simulation. Enhanced scalar transport problems, such as electrochemistry and levelset problems, instantiate
additional, problem specific singleton classes holding additional static parameters required for element evaluation. These
additional singleton classes are not meant to be derived from, but rather to coexist with this general class.

<pre>
\level 1

\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_PARAMETER_STD_H
#define SCATRA_ELE_PARAMETER_STD_H

#include "scatra_ele_parameter_base.H"

#include "../drt_inpar/inpar_scatra.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleParameterTimInt;

    // class implementation
    class ScaTraEleParameterStd : public ScaTraEleParameterBase
    {
      public:

        //! singleton access method
        static ScaTraEleParameterStd* Instance(
            const std::string&           disname,            //!< name of discretization
            const ScaTraEleParameterStd* delete_me = NULL  //!< creation/destruction indication
            );

        //! singleton destruction
        virtual void Done();

        //! set parameters
        void SetParameters(
            Teuchos::ParameterList& parameters   //!< parameter list
            );

        //! @name return general parameters
        //! @{
        bool IsAle() const {return is_ale_;};
        bool IsConservative() const {return is_conservative_;};
        bool SphericalCoords() const {return sphericalcoords_;};
        INPAR::SCATRA::FluxType CalcFluxDomain() const {return calcflux_domain_;};
        const Teuchos::RCP<std::vector<int> > WriteFluxIds() const {return writefluxids_;};
        INPAR::SCATRA::FDCheck FDCheck() const {return fdcheck_;};
        double FDCheckEps() const {return fdcheckeps_;};
        double FDCheckTol() const {return fdchecktol_;};
        int ProbNum() const {return probnum_;};
        //! @}

        //! @name return stabilization parameters
        //! @{
        INPAR::SCATRA::StabType StabType() const {return stabtype_;};
        INPAR::SCATRA::TauType TauDef() const {return whichtau_;};
        INPAR::SCATRA::CharEleLength CharEleLength() const {return charelelength_;};
        double USFEMGLSFac() const {return diffreastafac_;};
        bool RBSubGrVel() const {return sgvel_;};
        bool ASSGD() const {return assgd_;};
        INPAR::SCATRA::AssgdType ASSGDType() const {return whichassgd_;};
        bool TauGP() const {return tau_gp_;};
        bool MatGP() const {return mat_gp_;};
        double TauValue() const {return tau_value_;};
        //! @}

      private:

        //! private constructor for singletons
        ScaTraEleParameterStd(
            const std::string& disname   //!< name of discretization
            );

        //! private destructor for singletons
        virtual ~ScaTraEleParameterStd(){return;};

        //! @name general parameters
        //! @{
        //! flag for ALE
        bool is_ale_;

        //! flag for conservative form
        bool is_conservative_;

        //! flag for use of spherical coordinates
        bool sphericalcoords_;

        //! flag for writing the flux vector fields
        INPAR::SCATRA::FluxType calcflux_domain_;

        //! ids of scalars for which flux vectors are written (starting with 1)
        Teuchos::RCP<std::vector<int> > writefluxids_;

        //! flag for finite difference check
        INPAR::SCATRA::FDCheck fdcheck_;

        //! perturbation magnitude for finite difference check
        double fdcheckeps_;

        //! relative tolerance for finite difference check
        double fdchecktol_;
        //! @}

        //! @name stabilization parameters
        //! @{
        //! type of stabilization
        INPAR::SCATRA::StabType stabtype_;

        //! definition of stabilization parameter
        INPAR::SCATRA::TauType whichtau_;

        //! definition of characteristic element length
        INPAR::SCATRA::CharEleLength charelelength_;

        //! parameter to switch between SUPG, GLS and USFEM
        double diffreastafac_;

        //! flag to include residual-based subgrid-scale velocity
        bool sgvel_;

        //! flag to active artificial diffusion
        bool assgd_;

        //! definition of artificial diffusion
        INPAR::SCATRA::AssgdType whichassgd_;

        //! flag for evaluation of tau at Gauss point
        bool tau_gp_;

        //! flag for evaluation of material at Gauss point
        bool mat_gp_;

        //! numerical value for stabilization parameter
        double tau_value_;

        //! problem number
        int probnum_;
        //! @}

        //! parameter class for time integration
        DRT::ELEMENTS::ScaTraEleParameterTimInt* scatraparatimint_;
    };
  }
}
#endif
